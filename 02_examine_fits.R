library(tidyverse)
library(spOccupancy)
#library(coda)
library(bayesplot)
library(ggpubr)

rm(list=ls())
select <- dplyr::select
theme_set(theme_classic())

#-------------
# Examine fits
#-------------


# MODEL 3: 2022_04_11 --------------------------------------------------
#occ.ms.formula <- ~ . 
#det.ms.formula <- ~ . 
#species_sub1 <- c("Tasmanian Devil", "cat", "Bennett's Wallaby", 
#                  "Tasmanian Pademelon", "Spotted-tail Quoll")  

run_date <- "2022_04_11"
run_id <- ""
save_dir <- paste0("results/m_",run_date,"_",run_id);save_dir
fit <- readRDS(paste0(save_dir,"/fit_",run_date,".rds"))

summary(fit)

# select model: alpha (detection) or beta (occupancy) model
var_type <- "beta"

c("Devil","Cat", "Wallaby","Pademelon","Quoll") %>% 
  map(~ fit[[paste0(var_type,".samples")]] %>% 
        as.data.frame %>% 
        select(contains(.x)) %>% 
        rename_with(~ .x %>% str_remove("\\-.*")) %>% 
        mcmc_areas() + 
        xlim(-20,20) +
        labs(title = str_remove(.x,".* "))) %>% 
  ggarrange(plotlist = .) %>% 
  annotate_figure(top = paste(if_else(var_type=="alpha", "Detection","Occupancy"), "model") %>% 
                    text_grob(size = 14, face = "bold"))

#ggsave(paste0(save_dir,"/plot_",if_else(var_type=="alpha", "Detection","Occupancy"),".jpg"))







# MODEL 2: 2022_04_07 --------------------------------------------------
#occ.ms.formula <- ~ . - site 
#det.ms.formula <- ~ . - site 
#species_sub1 <- c("Tasmanian Devil", "cat", "Bennett's Wallaby", 
#                  "Tasmanian Pademelon", "Spotted-tail Quoll")  

run_date <- "2022_04_07"
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

fit$beta.samples %>% as.data.frame %>% colnames()
  rename_with(~ .x %>% str_remove("-Tasmanian Devil"))

plots <- c("Devil","Cat", "Wallaby","Pademelon","Quoll") %>% 
  map(~ fit$beta.samples %>% as.data.frame %>% 
         select(contains(.x)) %>% 
         rename_with(~ .x %>% str_remove("\\-.*")) %>% 
         mcmc_areas() + 
         xlim(-20,20) +
         labs(title = str_remove(.x,".* ")))  
         #ggsave(paste0(save_dir,"/beta_",.x,".pdf"), plot = .))

ggpubr::ggarrange(plotlist = plots)
fit$X.p
