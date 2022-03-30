#----------------------------------------------------------------------------
# 
# Import survey and camera-operation data for Jane River and
# create csv for daily occupancy by site and species.
# csv is long form: i.e, one row per (site, day, species) combination
#
# Author: Luke Yates
# Date created: 24th March, 2022
#
#----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(ggspatial)

rm(list=ls())

select <- dplyr::select
theme_set(theme_classic())

#source("camelot_analysis_func.r")
source("99_functions.R")



jr_raw <- read_csv("data/survey-export_2022-02-25_1551.csv")
jr_obs <- pre_process(jr_raw)

# load and summarise camera operation times
opt_time <- read_csv("data/jr-opt-time.csv")

jr_op_days <- process_op_days(opt_time)
jr_op_days %>% 
  group_by(cam) %>% 
  summarise(n = sum(op,na.rm = T)) %>% pull(n)

# merge daily occupancy detection with operational days
# creates one row for each (camera, day, species) combination
# occ: 0 = species not-detected; 1 = species detected; NA = camera not operational
jr <- jr_obs %>% 
  transmute(cam, date = date(date.time), common, occ = 1) %>% 
  distinct() %>% 
  group_by(common) %>% 
  group_modify(~ .x %>% 
        full_join(jr_op_days, by = c("date","cam"))) %>% 
  mutate(occ = case_when(op==1 ~ replace_na(occ,0)), op = NULL) %>% 
  relocate(cam,date) %>% 
  filter(cam != "JR-C17") %>% # remove faulty camera with only 1 day of operation
  arrange(cam,date) %>% 
  ungroup()

#write_csv(jr,"data/jr_daily_occ_2022_03_24.csv")




jr_spatial <- jr_obs %>% 
  distinct(cam, lat, lon) %>% 
  add_geometry()
jr_spatial
#write_csv(jr_spatial,"data/jr_spatial_2022_03_24.csv")


# take a look at the data
jr_spatial %>% 
  ggplot() +
  geom_sf(aes(col = cam)) + 
  geom_text(aes(x = 421000, y = y, label = cam), size = 3, hjust = 0, nudge_x = 300) +
  xlim(c(415000, 424500)) +
  labs(title = "Camera locations", col = "Camera", x = "", y = "") +
  theme_void() +
  theme(legend.position = "none", 
        panel.border = element_rect(linetype = 1, size = 1, fill = NA)) 

jr %>% 
  group_by(cam,date) %>% 
  summarise(n = min(sum(occ,na.rm = T),1)) %>% 
  group_by(cam) %>% 
  summarise(n = sum(n)) %>% 
  right_join(jr_obs %>% distinct(cam, lat, lon)) %>% 
  add_geometry() %>% 
  ggplot() +
  geom_sf(aes(size = n)) + 
  geom_text(aes(x = 421000, y = y, label = cam), size = 3, hjust = 0, nudge_x = 300) +
  xlim(c(415000, 424500)) +
  labs(title = "Camera locations", size = "# species-days-detected", x = "", y = "") +
  theme_void() +
  theme(legend.position = "right", 
        panel.border = element_rect(linetype = 1, size = 1, fill = NA)) 

#ggsave("plots/jr_sites_activity.pdf", width =4, height = 7)
