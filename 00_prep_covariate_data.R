
library(tidyverse)
library(spOccupancy)
library(forcats)


rm(list=ls())
select <- dplyr::select
theme_set(theme_classic())

#------------------------------------------------------------------
# Load occupancy data and shape into 3D array
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

#----------------------------------
# Load field-of-view covariate data
#----------------------------------

# FOV
# st.type: collapse road and w.track to track
# canopy/sub.canopy/shrub/ground sum values and use <=4 as threshold for low vs high cover 
# b.ground: delete
# l.litter: delete
# w.litter: delete
# rock: collapse (1,2,3) together
# c.b.grass: delete
# c.g.sedge: delete
# tree.fern: delete
# fern: delete
# moss: collapse (0-1) and (2-3)
# water: delete
# slope: keep
# treefalls: collapse (minor + major)  (none)
# scratch: delete
# gt.no: delete

jr_fov_veg_raw <- readxl::read_xlsx("data/jr-vegdata.xlsx","JR-FOV_veg") 

rem_vars <- c("b.ground","l.litter","w.litter","carpet/button.grass", 
              "cut.grass.sedge", "fern","water","scratch","gt.no")
cov_vars <- c("canopy","sub.canopy","shrub", "ground")

jr_fov_veg <- 
  jr_fov_veg_raw %>% 
  mutate(site = factor(camera, levels = paste0("JR-C",1:25), 
                         labels = paste0("JR-C",str_pad(1:25, 2, pad = 0)))) %>% 
  filter(site %in% sites) %>% 
  select(-Notes,-camera, -all_of(rem_vars)) %>% 
  mutate(cover = {rowSums(across(all_of(cov_vars),as.numeric))>=4} %>% 
           factor(levels = c("FALSE","TRUE"),labels = c("0","1"))) %>% # low/high
  select(-all_of(cov_vars)) %>% 
  mutate(across(any_of(c("h.cam","h.mid","l.fov","l.trk","gt.no")),as.numeric)) %>% 
  mutate(across(where(is.character), factor)) %>% 
  mutate(track = fct_collapse(st.type, `0` = c("arena"),`1` = c("road","g.trail","w.track"))) %>% 
  select(-st.type) %>% 
  mutate(rock = fct_collapse(rock, `1` = c("1","2","3"))) %>% 
  mutate(moss = fct_collapse(moss, `0` = c("0","1"), `1` = c("2","3"))) %>% #low/high
  mutate(treefalls = fct_collapse(treefalls, `0` = c("none"), `1` = c("minor","major"))) %>% 
  mutate(slope = fct_collapse(slope, `0` = c("flat"), `1` = c("low"))) # flat/low




#----------------------------------
# Load broad area covariate data
#----------------------------------

# BA
# canopy: delete
# sub.canopy: collapse (0,1) (2-)
# shrub: collapse (0,1) (2-)
# ground: collapse (0,1) (2-)
# b.ground: collapse (0,1) (2-)
# l.litter: delete
# w.litter:  collapse (0,1) (2-)
# rock: delete
# c.b.grass: (0) (1-)
# c.g.grass: (0) (1-)
# tree.fern: delete
# fern: (0) (1-)
# moss: (0,1) (2,3)
# water: (none) (still, flowing, muddy)
# slope: (flat, low) (mod,steep)
# treefalls: (major, minor) (none)
# scratch: delete

jr_ba_veg_raw <- readxl::read_xlsx("data/jr-vegdata.xlsx","JR-BROAD_veg")

jr_ba_veg_ <- jr_ba_veg_raw %>% 
  mutate(camera = factor(camera, levels = paste0("JR-C",1:25), 
                         labels = paste0("JR-C",str_pad(1:25, 2, pad = 0)))) %>% 
  select(-Notes, -canopy,-l.litter,-rock, - tree.fern, -scratch)

# replace broad area NAs with field of view values for unmeasured site JR-C16 
jr_ba_veg_[16,] <- rename(jr_fov_veg_raw, `carpet.grass` = `carpet/button.grass`)[16,names(jr_ba_veg_[16,])]

jr_ba_veg <- 
  jr_ba_veg_ %>% 
  filter(camera %in% sites) %>% 
  mutate(sub.canopy = fct_collapse(sub.canopy, `0` = c("0","1"), other_level = "1"),
         shrub = fct_collapse(shrub, `0` = c("0","1"), other_level = "1"),
         ground = fct_collapse(ground, `0` = c("0","1"), other_level = "1"),
         b.ground = fct_collapse(b.ground, `0` = c("0","1"), other_level = "1"),
         w.litter = fct_collapse(w.litter, `0` = c("0","1"), other_level = "1"),
         carpet.grass = fct_collapse(carpet.grass, `0` = c("0"), other_level = "1"), 
         cut.grass.sedge = fct_collapse(cut.grass.sedge, `0` = c("0"), other_level = "1"), 
         fern = fct_collapse(fern, `0` = c("0"), other_level = "1"),
         moss = fct_collapse(moss, `0` = c("0","1"), other_level = "1"), 
         water = fct_collapse(water, "0" = c("none"),other_level = "1"), 
         slope = fct_collapse(slope, "0" = c("flat","low"), other_level = "1"), 
         treefalls = fct_collapse(treefalls, `0` = c("none"), other_level = "1")) %>% 
  rename(site = camera, grass = carpet.grass, sedge = cut.grass.sedge, cover = sub.canopy)
  
# Load spatial data
jr_spatial <- read_csv("data/jr_spatial_2022_03_24.csv") %>% # lat/lon plus projected coords x/y
  filter(cam %in% sites) %>% 
  select(X = x, Y = y)

# check that order of site names match: TRUE
identical(dimnames(y)$site,as.character(jr_fov_veg$site))
identical(dimnames(y)$site,as.character(jr_ba_veg$site))



#--------------------------------------------
# collate data into list for spOccupancy etc.
#--------------------------------------------

jr_occ_data <- list(y = y,
                occ.covs = jr_ba_veg %>% select(-site),
                det.covs = jr_fov_veg %>% select(-site) %>% as.list %>% map(as.vector),
                coords = jr_spatial)

str(jr_occ_data)

#saveRDS(jr_occ_data, "data/jr_occ_data_2022_04_07.rds")
