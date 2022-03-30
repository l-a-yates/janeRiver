library(tidyverse)
library(spOccupancy)
library(ompr)
library(spNNGP)
library(coda)

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

# import and clean: detection covariates, occupancy covariates, and spatial data
jr_fov_veg <- readxl::read_xlsx("data/jr-vegdata.xlsx","JR-FOV_veg") %>% # field of view (fov)
  select(-Notes) %>%
  mutate(camera = factor(camera, levels = paste0("JR-C",1:25), 
                         labels = paste0("JR-C",str_pad(1:25, 2, pad = 0)))) %>% 
  filter(camera %in% sites) %>% 
  rename(c.b.grass = `carpet/button.grass`,
         c.g.sedge = cut.grass.sedge,
         site = camera) %>%
  mutate(across(any_of(c("h.cam","h.mid","l.fov","l.trk","gt.no")),as.numeric)) %>% 
  mutate(across(where(is.character), factor))


# canopy/sub.canopy/shrub/ground combine to 0 (low cover) 1 (high cover) 
jr_fov_veg %>% 
  select(canopy,sub.canopy,shrub, ground) %>% 
  mutate(across(everything(), as.numeric))  %>% 
  mutate(across(everything(), ~ .x -1)) %>% 
  mutate(cover = (canopy + sub.canopy + shrub + ground)>=4)

# FOV
# st.type: collapse road and w.track to track
# canopy/sub.canopy/shrub/ground combine to 0 (low cover) 1 (high cover) 
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



jr_fov_veg %>% 
  select(where(is.factor), -site) %>% 
  map(table)




#----------------------------------
# Load broad area covariate data
#----------------------------------



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

jr_ba_veg %>% 
  select(-site) %>% 
  mutate(across(everything(), ~ .x %>% as.character() %>% as.factor)) %>% 
  map(table)


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

# Load spatial data
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
