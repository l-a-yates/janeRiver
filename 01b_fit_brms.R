
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
