library(tidyverse)
library(spOccupancy)
#library(coda)
library(bayesplot)
library(ggpubr)

rm(list=ls())
select <- dplyr::select
theme_set(theme_classic())

jr_data <- readRDS("data/jr_occ_data_2022_04_07.rds")

y <- jr_data$y

apply(y,c(1,2),sum, na.rm = T) %>% 
  knitr::kable()

# data dimension names
sites <- dimnames(y)$site; sites
species <- dimnames(y)$species; species
visits <- dimnames(y)$visit; visits

# subset data
species_sub <- c("Tasmanian Devil", 
                 "Cat",
                 "Bennett's Wallaby",
                 "Tasmanian Pademelon", 
                 "Spotted-tail Quoll")  

# look at each species across all sites
y %>% 
  reshape2::melt(value.name = "occ") %>%
  as_tibble() %>% 
  filter(species %in% species_sub) %>% 
  mutate(occ = factor(occ)) %>% 
  ggplot() +
  geom_raster(aes(x = visit, y = site, fill = occ)) +
  scale_fill_manual(values = c(`0` = "grey10",`1` = "blue")) +
  facet_wrap(~ species) +
  theme(legend.position = "none")

#ggsave("plots/dataVis/siteOcc_by_species.pdf")

# look at each site acros species
y %>% 
  reshape2::melt(value.name = "occ") %>%
  as_tibble() %>% 
  #filter(site == "JR-C01") %>% 
  #filter(site %in% sites[1:6]) %>% 
  filter(species %in% species_sub) %>% 
  mutate(occ = factor(occ)) %>% 
  ggplot() +
  geom_raster(aes(x = visit, y = species, fill = occ)) +
  scale_fill_manual(values = c(`0` = "grey10",`1` = "blue")) +
  facet_wrap(~ site) +
  theme(legend.position = "none")

#ggsave("plots/dataVis/speciesOcc_by_site.pdf")

# occupancy covariate values (binary)
jr_data$occ.covs %>% 
  mutate(site = sites) %>% 
  pivot_longer(-site, names_to = 'predictor', values_to = "value") %>% 
  ggplot() +
  geom_raster(aes(x = predictor,y = site, fill = value)) +
  scale_fill_manual(values = c(`0` = "grey10",`1` = "blue")) +
  theme(legend.position = "none") +
  labs(title = "Occupancy covariates")

#ggsave("plots/dataVis/occ_covariates.pdf")

# detection covariate values (binary)
jr_data$det.covs %>% 
  bind_cols() %>% 
  select(where(is_character)) %>% 
  mutate(site = sites) %>% 
  pivot_longer(-site, names_to = 'predictor', values_to = "value") %>% 
  ggplot() +
  geom_raster(aes(x = predictor,y = site, fill = value)) +
  scale_fill_manual(values = c(`0` = "grey10",`1` = "blue")) +
  theme(legend.position = "none") +
  labs(title = "Detection covariates")

#ggsave("plots/dataVis/det_covariates.pdf")

