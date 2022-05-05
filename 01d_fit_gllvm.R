library(tidyverse)
library(gllvm)

rm(list=ls())
select <- dplyr::select
theme_set(theme_classic())

# set run details
run_date <- "2022_05_05"
run_id <- ""

# create save directory
save_dir <- paste0("results/m_",run_date,"_",run_id)
if(!dir.exists(save_dir)) dir.create(save_dir)

jr_data <- readRDS("data/jr_occ_data_2022_04_07.rds")


