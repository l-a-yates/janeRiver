library(tidyverse)
library(spOccupancy)
library(ompr)
library(spNNGP)
library(coda)

rm(list=ls())
select <- dplyr::select
theme_set(theme_classic())

#------------------------------------------------------------------
# Data preparation for spOccupancy
#------------------------------------------------------------------

# import occupancy data
jr <- read_csv("data/jr_daily_occ_2022_03_24.csv")

# create 3D occupancy array
y <- jr %>% 
  arrange(date) %>% 
  select(date) %>% 
  distinct %>% 
  mutate(visit = 1:n()) %>% 
  right_join(jr, by = "date") %>% 
  select(species = common, site = cam, visit, occ) %>% 
  xtabs(occ ~ species + site + visit, data = ., na.action = na.pass)

# data dimension names
sites <- dimnames(y)$site; sites
species <- dimnames(y)$species
visits <- dimnames(y)$visit

# check that reshaped occupancy values match by melting array back to long format: YES
identical(reshape2::melt(y, value.name = "occ") %>%
            as_tibble() %>%
            arrange(site,visit,species) %>% pull(occ),
            jr %>% arrange(cam,date,common)  %>% pull(occ))

# import and clean: detection covariates, occupancy covariates, and spatial data
jr_fov_veg <- readxl::read_xlsx("data/jr-vegdata.xlsx","JR-FOV_veg") %>% # field of view (fov)
  select(-Notes) %>%
  mutate(camera = factor(camera, levels = paste0("JR-C",1:25), 
                         labels = paste0("JR-C",str_pad(1:25, 2, pad = 0)))) %>% 
  filter(camera %in% sites) %>% 
  rename(c.b.grass = `carpet/button.grass`,
         c.g.sedge = cut.grass.sedge,
         site = camera) %>%
  mutate(across(any_of(c("h.cam","h.mid","l.fov","l.trk")),as.numeric)) %>% 
  mutate(across(where(is.character), factor))

jr_ba_veg <- readxl::read_xlsx("data/jr-vegdata.xlsx","JR-BROAD_veg") %>%  # broad area (ba)
  select(-Notes) %>% 
  mutate(camera = factor(camera, levels = paste0("JR-C",1:25), 
                         labels = paste0("JR-C",str_pad(1:25, 2, pad = 0)))) %>% 
  filter(camera %in% sites) %>% 
  rename(c.b.grass = carpet.grass,
         c.g.sedge = cut.grass.sedge,
         site = camera) %>% 
  mutate(across(where(is.character), factor))

# replace broad area NAs with field of view values for JR-C16 
jr_ba_veg[16,] <- jr_fov_veg[16,names(jr_ba_veg[16,])]

jr_spatial <- read_csv("data/jr_spatial_2022_03_24.csv") %>% # lat/lon plus projected coords x/y
  filter(cam %in% sites) %>% 
  select(X = x, Y = y)

# check that order of site names match: TRUE
identical(dimnames(y)$site,as.character(jr_fov_veg$site))
identical(dimnames(y)$site,as.character(jr_ba_veg$site))


# collate data into list for spOccupancy
jr_data <- list(y = y,
                occ.covs = jr_ba_veg %>% select(-site),
                det.covs = jr_fov_veg %>% select(-site) %>% as.list,
                coords = jr_spatial)

str(jr_data)

#------------------------------------------------------------------
# Model fitting
#------------------------------------------------------------------

occ.ms.formula <- ~ . 
det.ms.formula <- ~ . 

species_sub <- c("Tasmanian Devil", "Cat", "Bennett's Wallaby", "Tasmanian Pademelon", "Spotted-tail Quoll")  

run_date <- "2022_03_30"
run_id <- ""
save_dir <- paste0("results/m_",run_date,"_",run_id)
if(!dir.exists(save_dir)) dir.create(save_dir)

jr_spOcc <- jr_data
jr_spOcc$y <- jr_spOcc$y[species_sub,,]

ms.inits <- list(alpha.comm = 0, 
                 beta.comm = 0, 
                 beta = 0, 
                 alpha = 0,
                 tau.sq.beta = 1, 
                 tau.sq.alpha = 1, 
                 z = apply(jr_spOcc$y, c(1, 2), max, na.rm = TRUE))

ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72), 
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

fit <- msPGOcc(occ.formula = occ.ms.formula, 
                  det.formula = det.ms.formula, 
                  data = jr_spOcc, 
                  inits = ms.inits, 
                  n.samples = 30000, 
                  priors = ms.priors, 
                  n.omp.threads = 3, 
                  verbose = TRUE, 
                  n.report = 1500, 
                  n.burn = 10000,
                  n.thin = 50, 
                  n.chains = 3)

saveRDS(fit,paste0(save_dir,"/fit_",run_date,".rds"))


### BRMS
if(F){
library(brms)

jr_long <- 
  jr %>% 
  arrange(date) %>% 
  select(date) %>% 
  distinct %>% 
  mutate(visit = 1:n()) %>% 
  right_join(jr, by = "date") %>% 
  select(species = common, site = cam, visit, occ) %>% 
  left_join(jr_fov_veg %>% select(site, fov.canopy = canopy, fov.c.b.grass = c.b.grass), by = "site") %>% 
  left_join(jr_ba_veg %>% select(site, ba.canopy = canopy, ba.water = water, ba.rock = rock), by = "site")

species_sub <- c("Tasmanian Devil", "Bennett's Wallaby", "Tasmanian Pademelon", "Spotted-tail Quoll")  

jr_brms <- 
  jr_long %>% 
  filter(species %in% species_sub) %>% 
  pivot_wider(values_from = occ, names_from = species) %>% 
  rename(y_D = "Tasmanian Devil", 
         y_W = "Bennett's Wallaby", 
         y_P = "Tasmanian Pademelon", 
         y_Q = "Spotted-tail Quoll") %>% 
  mutate(obs = 1:n())

bf(y_D ~ p*z, 
   p ~ fov.canopy + fov.c.b.grass, 
   z ~ ba.canopy + ba.water + ba.rock + 1|a|obs, nl = TRUE)

form <- bf(mvbind(y_D,y_W,y_P, y_Q) ~ inv_logit(p)*inv_logit(z), 
            p ~ 1 + fov.canopy + fov.c.b.grass, 
            z ~ 1 + ba.canopy + ba.water + ba.rock, nl = TRUE) + 
  bernoulli() +
  set_rescor(FALSE)

form

for(nl in c("p","z")) for(resp in c("yD","yW","yP","yQ")){
  if(nl=="p" & resp == "yD"){
    prior_1 <- prior_string("normal(0, 5)", nlpar = "p", resp = "yD")
  } else {
    prior_1 = prior_1 + prior_string("normal(0, 5)", nlpar = nl, resp = resp)}
}
prior_1

make_stancode(form, jr_brms, prior = prior_1)

fit.brms.1 <- brm(form, data = jr_brms, chains = 4, cores = 4, prior = prior_1)
fit.brms.1
}
