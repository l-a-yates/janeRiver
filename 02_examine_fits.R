library(tidyverse)
library(spOccupancy)
library(coda)
library(bayesplot)

rm(list=ls())
select <- dplyr::select
theme_set(theme_classic())

#-------------
# Examine fits
#-------------


# MODEL 1: 2022_03_28 --------------------------------------------------
#occ.ms.formula <- ~ . - site 
#det.ms.formula <- ~ . - site 
#species_sub1 <- c("Tasmanian Devil", "Bennett's Wallaby", 
#                  "Tasmanian Pademelon", "Spotted-tail Quoll")  

run_date <- "2022_03_28"
run_id <- ""
save_dir <- paste0("results/m_",run_date,"_",run_id);save_dir
fit <- readRDS(paste0(save_dir,"/fit_",run_date,".rds"))


# community-level estimates
fit$beta.comm.samples %>% mcmc_areas # occupancy 
fit$alpha.comm.samples %>% mcmc_areas()  # detection

# species-level estimates
fit$alpha.samples %>% mcmc_areas
fit$beta.samples %>% colnames

c("Devil","Wallaby","Pademelon","Quoll") %>% 
  walk(~ {fit$beta.samples %>% as.data.frame %>% 
         select(contains(.x)) %>% 
         mcmc_areas() + labs(title = .x)} %>% 
         ggsave(paste0(save_dir,"/beta_",.x,".pdf"), plot = .))


