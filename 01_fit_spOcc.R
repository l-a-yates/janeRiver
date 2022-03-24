library(tidyverse)
library(spOccupancy)

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
  mutate(camera = factor(camera, levels = paste0("JR-C",1:25), labels = paste0("JR-C",str_pad(1:25, 2, pad = 0)))) %>% 
  filter(camera %in% sites) %>% 
  rename(c.b.grass = `carpet/button.grass`,
         c.g.sedge = cut.grass.sedge,
         site = camera) %>%
  mutate(across(any_of(c("h.cam","h.mid","l.fov","l.trk")),as.numeric)) %>% 
  mutate(across(where(is.character), factor))

jr_ba_veg <- readxl::read_xlsx("data/jr-vegdata.xlsx","JR-BROAD_veg") %>%  # broad area (ba)
  select(-Notes) %>% 
  mutate(camera = factor(camera, levels = paste0("JR-C",1:25), labels = paste0("JR-C",str_pad(1:25, 2, pad = 0)))) %>% 
  filter(camera %in% sites) %>% 
  rename(c.b.grass = carpet.grass,
         c.g.sedge = cut.grass.sedge,
         site = camera) %>% 
  mutate(across(where(is.character), factor))

jr_spatial <- read_csv("data/jr_spatial_2022_03_24.csv") %>% # lat/lon plus projected coords x/y
  filter(cam %in% sites) %>% 
  select(X = x, Y = y)

# collate data into list for spOccupancy
jr_data <- list(y = y,
                occ.covs = jr_ba_veg,
                det.covs = jr_fov_veg,
                coords = jr_spatial)

str(jr_data)


#------------------------------------------------------------------
# Model fitting
#------------------------------------------------------------------

occ.ms.formula <- ~ canopy + shrub + rock + water + l.litter
det.ms.formula <- ~ canopy + shrub + rock + h.mid

ms.inits <- list(alpha.comm = 0, 
                 beta.comm = 0, 
                 beta = 0, 
                 alpha = 0,
                 tau.sq.beta = 1, 
                 tau.sq.alpha = 1, 
                 z = apply(jr_data$y, c(1, 2), max, na.rm = TRUE))

ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72), 
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

out.ms <- msPGOcc(occ.formula = occ.ms.formula, 
                  det.formula = det.ms.formula, 
                  data = jr_data, 
                  inits = ms.inits, 
                  n.samples = 30000, 
                  priors = ms.priors, 
                  n.omp.threads = 3, 
                  verbose = TRUE, 
                  n.report = 6000, 
                  n.burn = 10000,
                  n.thin = 50, 
                  n.chains = 3)

