library(tidyverse)
library(spOccupancy)
library(ompr)
library(spNNGP)
library(coda)

rm(list=ls())
select <- dplyr::select
theme_set(theme_classic())

# set run details
run_date <- "2022_04_11"
run_id <- ""

# create save directory
save_dir <- paste0("results/m_",run_date,"_",run_id)
if(!dir.exists(save_dir)) dir.create(save_dir)

jr_data <- readRDS("data/jr_occ_data_2022_04_07.rds")

occ.ms.formula <- as.formula(paste0("~ ",paste0(colnames(jr_data$occ.covs), collapse = " + "))) 
det.ms.formula <- as.formula(paste0("~ ",paste0(names(jr_data$det.covs), collapse = " + "))) 


# subset data
species_sub <- c("Tasmanian Devil", 
                 "Cat", 
                 "Bennett's Wallaby", 
                 "Tasmanian Pademelon", 
                 "Spotted-tail Quoll")  

jr_spOcc <- jr_data
jr_spOcc$y <- jr_spOcc$y[species_sub,,]


# model summary
message("occupancy model:")
message(occ.ms.formula)
message("detection model:")
message(det.ms.formula)
message("Species:")
message(paste(species_sub, collapse = ", "))


# set initial values
ms.inits <- list(alpha.comm = 0, 
                 beta.comm = 0, 
                 beta = 0, 
                 alpha = 0,
                 tau.sq.beta = 1, 
                 tau.sq.alpha = 1, 
                 z = apply(jr_spOcc$y, c(1, 2), max, na.rm = TRUE))

# set priors
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                  alpha.comm.normal = list(mean = 0, var = 2.72), 
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

# fit model
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

