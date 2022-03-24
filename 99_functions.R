#---------------------------------------------------
#
# Luke's functions
#
#---------------------------------------------------

pre_process <- function(data) data %>% 
  select(cam = `Trap Station Name`,
         lon  = `Camelot GPS Longitude`,
         lat = `Camelot GPS Latitude`,
         date.time = `Date/Time`,
         flash = Flash,
         class = Class,
         order = Order,
         common = `Species Common Name`,
         scientific = Species...52,
         count = `Sighting Quantity`,
         life.stage = `Life stage`,
         sex = Sex,
         iso = `ISO Speed Ratings`,
         image = `Media ID`) %>% 
  mutate(flash = flash %>% as.factor %>% as.numeric %>% {.-1}, # 1 = flash fired (night shot), 0 = day shot as a numeric
         life.stage = replace_na(life.stage, "Adult"), # blank life-stage entries get filled as Adult
         sex = replace_na(sex, "Unknown"), # blank sex entries get filled as Unknown
         iso = case_when(iso<=50 ~ iso), # Exclude images that were not pre-classified
         cam = factor(cam, levels = paste0("JR-C",1:25)), # order camera levels correctly
         date.time = ymd_hms(date.time, tz = "Australia/Queensland")) %>% # use QLD tz because no daylight saving
  arrange(cam, date.time)



add_geometry <- function(x) x %>% 
  sf::st_as_sf(coords = c("lon","lat"), crs = 4326) %>%  # 4326 is the EPS code for WS84
  sf::st_transform(crs = 32755) %>% # 32755 is the EPS code for UTM zone 55S (South) 
  mutate(X= sf::st_coordinates(geometry)[,"X"],
         Y= sf::st_coordinates(geometry)[,"Y"])


process_op_days <- function(op_data) op_data %>% 
  transmute(cam = camera, 
            start_date = dmy(s0.date),
            end_date = start_date + days(s1.rd)) %>% 
  drop_na %>% # remove cameras with no run days
  group_split(cam) %>% 
  imap_dfr(~ tibble(cam = .x$cam[1], # generate date sequence of active operation days
                    date = seq(ymd(.x$start_date[1]),ymd(.x$end_date[1]),by = "days"),
                    op = 1)) %>% 
  pivot_wider(names_from = date, values_from = op) %>% # generate entry for each day for each cam
  pivot_longer(-cam, names_to = "date", values_to = "op") %>% 
  mutate(date = ymd(date))







# Barry's original function
create_op_mat <- function(op) {
  n.cam <- length(op$camera) # number of cameras
  op <- op[,apply(op, 2, function(x) !all(is.na(x)))] # trim service columns where all cameras have NA
  n.service <- (ncol(op)-2)/2 # number of services
  date.cols <- seq(2,ncol(op),by=2) # columns with service dates
  
  for(c in date.cols) op[,c] <- as.Date(op[,c],format="%d/%m/%Y",tz="Australia/Queensland") # QLD tz b/c no daylight savings
  
  start.date <- as.Date(min(apply(op[,date.cols],2, min,na.rm=T),na.rm=T)) # Earliest date a camera was deployed
  end.date <- as.Date(max(apply(op[,date.cols],2,max, na.rm=T),na.rm=T)) # Most recent date a camera was serviced
  duration <- as.numeric(difftime(end.date,start.date,units='days'))+1 # Total days covered by op.df
  
  odf <- as.data.frame(matrix(0,nrow=duration,ncol=length(op$camera)+1)) # initialise df for operating days
  colnames(odf) <- c('date',as.character(op$camera)) # cameras as column names
  odf$date <- seq.Date(start.date,end.date,by="days") # create a sequence of dates (days) from first image to last
  
  for(c in 1:n.cam){
    for(s in 1:n.service) {
      sd <- which(odf$date==op[c,date.cols[s]]) # match service date to appropriate row in op.df
      rd <- op[c,date.cols[s]+1] # run days since last service
      if(!is.na(rd)) odf[sd:(sd+rd),which(colnames(odf)==op[c,1])] <- 1 # fill in service days unless it was the last service
    }
  } # Fill odf with 1 if cam operating on a given day, 0 if not. Used to calculate effort per camera.
  return(odf)
} # create matrix of camera operation by day




