library(tidyverse)
library(unmarked)

rm(list=ls())
select <- dplyr::select
theme_set(theme_classic())

# set run details
run_date <- "2022_04_11"
run_id <- "um"

# create save directory
save_dir <- paste0("results/m_",run_date,"_",run_id)
if(!dir.exists(save_dir)) dir.create(save_dir)

# load and organise data
jr_data <- readRDS("data/jr_occ_data_2022_04_07.rds")

n_visits <- dim(jr_data$y)[3] # 366
n_sites <- dim(jr_data$y)[2] # 24


# subset data
species_sub <- c("Tasmanian Devil", 
                 #"Cat",
                 #"Bennett's Wallaby",
                 #"Tasmanian Pademelon", 
                 "Spotted-tail Quoll")  

n_species <- length(species_sub)

jr_data_occMulti <- list() 
jr_data_occMulti$y <- jr_data$y[species_sub,,] %>% 
  reshape2::melt(value.name = "occ") %>% 
  split(unique(.$species))  %>% 
  map(~ .x %>%   
              pivot_wider(-1, names_from = visit, values_from = occ) %>% 
              select(-site) %>% 
              as.matrix)

#jr_data_occMulti$det.covs <-  jr_data$det.covs %>% map(~ .x %>% rep(n_visits) %>% matrix(n_sites,n_visits))
jr_data_occMulti$occ.covs <- jr_data$occ.covs

data_jr <- unmarkedFrameOccuMulti(y=jr_data_occMulti$y,
                                 siteCovs=jr_data_occMulti$occ.covs,
                                 obsCovs=NULL,
                                 maxOrder = 2)
data_jr %>% plot
summary(data_jr)

# Length should match number/order of columns in fDesign
#occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov3','~1','~1','~1','~1')
occFormulas <- rep("0", ncol(data_jr@fDesign)); occFormulas
occFormulas[1:(n_species)] <- "~1"

#occFormulas[1:(n_species*(n_species+1)*0.5)] <- "~1"


colnames(data_jr@fDesign)

#Length should match number/order of species in data@ylist
detFormulas <- rep("~1", length(data_jr@ylist));detFormulas #c('~1','~1','~1')

fit <- occuMulti(detformulas = detFormulas, 
                 stateformulas = occFormulas,
                 data = data_jr, 
                 se = F, 
                 maxOrder = 2, 
                 engine = "R",
                 penalty = 2)
 
fit_opt <- optimizePenalty(fit, k = 10, boot = 100)

fit_opt

gd <- getDes(data_jr, detformulas = detFormulas, stateformulas = occFormulas, maxOrder = 2)
gd$dmF


## Not run: 
#Simulate 3 species data
N <- 1000
nspecies <- 3
J <- 5

occ_covs <- as.data.frame(matrix(rnorm(N * 10),ncol=10))
names(occ_covs) <- paste('occ_cov',1:10,sep='')

det_covs <- list()
for (i in 1:nspecies){
  det_covs[[i]] <- matrix(rnorm(N*J),nrow=N)
}
names(det_covs) <- paste('det_cov',1:nspecies,sep='')

#True vals
beta <- c(0.5,0.2,0.4,0.5,-0.1,-0.3,0.2,0.1,-1,0.1)
f1 <- beta[1] + beta[2]*occ_covs$occ_cov1
f2 <- beta[3] + beta[4]*occ_covs$occ_cov2
f3 <- beta[5] + beta[6]*occ_covs$occ_cov3
f4 <- beta[7]
f5 <- beta[8]
f6 <- beta[9]
f7 <- beta[10]
f <- cbind(f1,f2,f3,f4,f5,f6,f7)
z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
colnames(z) <- paste('sp',1:nspecies,sep='')
dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)

psi <- exp(f %*% t(dm))
psi <- psi/rowSums(psi)

#True state
ztruth <- matrix(NA,nrow=N,ncol=nspecies)
for (i in 1:N){
  ztruth[i,] <- as.matrix(z[sample(8,1,prob=psi[i,]),])
}

p_true <- c(0.6,0.7,0.5)

# fake y data
y <- list()

for (i in 1:nspecies){
  y[[i]] <- matrix(NA,N,J)
  for (j in 1:N){
    for (k in 1:J){
      y[[i]][j,k] <- rbinom(1,1,ztruth[j,i]*p_true[i])
    }
  }
}
names(y) <- c('coyote','tiger','bear')

#Create the unmarked data object
data = unmarkedFrameOccuMulti(y=y,siteCovs=occ_covs,obsCovs=det_covs)

#Summary of data object
summary(data)
plot(data)

# Look at f parameter design matrix
data@fDesign

# Formulas for state and detection processes

# Length should match number/order of columns in fDesign
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov3','~1','~1','~1','~1')

#Length should match number/order of species in data@ylist
detFormulas <- c('~1','~1','~1')

fit <- occuMulti(detFormulas,occFormulas,data)

#Look at output
fit

plot(fit)

#Compare with known values
cbind(c(beta,log(p_true/(1-p_true))),fit@opt$par)

#predict method
lapply(predict(fit,'state'),head)
lapply(predict(fit,'det'),head)

#marginal occupancy
head(predict(fit,'state',species=2))
head(predict(fit,'state',species='bear'))
head(predict(fit,'det',species='coyote'))

#probability of co-occurrence of two or more species
(predict(fit, 'state', species=c('coyote','tiger')))

#conditional occupancy
head(predict(fit,'state',species=2,cond=3)) #tiger | bear present
head(predict(fit,'state',species='tiger',cond='bear')) #tiger | bear present
head(predict(fit,'state',species='tiger',cond='-bear')) #bear absent
head(predict(fit,'state',species='tiger',cond=c('coyote','-bear')))

#residuals (by species)
lapply(residuals(fit),head)

#ranef (by species)
ranef(fit, species='coyote')

#parametric bootstrap
bt <- parboot(fit,nsim=30)

#update model
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov2+occ_cov3','~1','~1','~1','~1')
fit2 <- update(fit,stateformulas=occFormulas)

#List of fitted models
fl <- fitList(fit,fit2)
coef(fl)

#Model selection
modSel(fl)

#Fit model while forcing some natural parameters to be 0
#For example: fit model with no species interactions
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov2+occ_cov3','0','0','0','0')
fit3 <- occuMulti(detFormulas,occFormulas,data)

#Alternatively, you can force all interaction parameters above a certain
#order to be zero with maxOrder. This will be faster.
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov2+occ_cov3')
fit4 <- occuMulti(detFormulas,occFormulas,data,maxOrder=1)

## End(Not run)



getDes <- function(umf, detformulas, stateformulas, maxOrder, na.rm=TRUE, warn=FALSE,
         newdata=NULL, type="state")
{
  
  #Format formulas
  #Workaround for parameters fixed at 0
  fixed0 <- stateformulas %in% c("~0","0")
  stateformulas[fixed0] <- "~1"
  
  stateformulas <- lapply(stateformulas,as.formula)
  detformulas <- lapply(detformulas,as.formula)
  
  #Generate some indices
  S <- length(umf@ylist) # of species
  if(missing(maxOrder)){
    maxOrder <- S
  }
  z <- expand.grid(rep(list(1:0),S))[,S:1] # z matrix
  colnames(z) <- names(umf@ylist)
  M <- nrow(z) # of possible z states
  
  # f design matrix
  if(maxOrder == 1){
    dmF <- as.matrix(z)
  } else {
    dmF <- model.matrix(as.formula(paste0("~.^",maxOrder,"-1")),z)
  }
  nF <- ncol(dmF) # of f parameters
  
  J <- ncol(umf@ylist[[1]]) # max # of samples at a site
  N <- nrow(umf@ylist[[1]]) # of sites
  
  #Check formulas
  if(length(stateformulas) != nF)
    stop(paste(nF,"formulas are required in stateformulas list"))
  if(length(detformulas) != S)
    stop(paste(S,"formulas are required in detformulas list"))
  
  if(is.null(siteCovs(umf))) {
    site_covs <- data.frame(placeHolderSite = rep(1, N))
  } else {
    site_covs <- siteCovs(umf)
  }
  
  if(is.null(obsCovs(umf))) {
    obs_covs <- data.frame(placeHolderObs = rep(1, J*N))
  } else {
    obs_covs <- obsCovs(umf)
  }
  
  #Add site covs to obs covs if we aren't predicting with newdata
  # Record future column names for obsCovs
  col_names <- c(colnames(obs_covs), colnames(site_covs))
  
  # add site covariates at observation-level
  obs_covs <- cbind(obs_covs, site_covs[rep(1:N, each = J),])
  colnames(obs_covs) <- col_names
  
  #Re-format ylist
  index <- 1
  ylong <- lapply(umf@ylist, function(x) {
    colnames(x) <- 1:J
    x <- cbind(x,site=1:N,species=index)
    index <<- index+1
    x
  })
  ylong <- as.data.frame(do.call(rbind,ylong))
  ylong <- reshape(ylong, idvar=c("site", "species"), varying=list(1:J),
                   v.names="value", direction="long")
  ylong <- reshape(ylong, idvar=c("site","time"), v.names="value",
                   timevar="species", direction="wide")
  ylong <- ylong[order(ylong$site, ylong$time), ]
  
  #Remove missing values
  if(na.rm){
    naSiteCovs <- which(apply(site_covs, 1, function(x) any(is.na(x))))
    if(length(naSiteCovs>0)){
      stop(paste("Missing site covariates at sites:",
                 paste(naSiteCovs,collapse=", ")))
    }
    
    naY <- apply(ylong, 1, function(x) any(is.na(x)))
    naCov <- apply(obs_covs, 1, function(x) any(is.na(x)))
    navec <- naY | naCov
    
    sites_with_missingY <- unique(ylong$site[naY])
    sites_with_missingCov <- unique(ylong$site[naCov])
    
    ylong <- ylong[!navec,,drop=FALSE]
    obs_covs <- obs_covs[!navec,,drop=FALSE]
    
    no_data_sites <- which(! 1:N %in% ylong$site)
    if(length(no_data_sites>0)){
      stop(paste("All detections and/or detection covariates are missing at sites:",
                 paste(no_data_sites,collapse=", ")))
    }
    
    if(sum(naY)>0&warn){
      warning(paste("Missing detections at sites:",
                    paste(sites_with_missingY,collapse=", ")))
    }
    if(sum(naCov)>0&warn){
      warning(paste("Missing detection covariate values at sites:",
                    paste(sites_with_missingCov,collapse=", ")))
    }
    
  }
  
  #Start-stop indices for sites
  yStart <- c(1,1+which(diff(ylong$site)!=0))
  yStop <- c(yStart[2:length(yStart)]-1,nrow(ylong))
  
  y <- as.matrix(ylong[,3:ncol(ylong)])
  
  #Indicator matrix for no detections at a site
  Iy0 <- do.call(cbind, lapply(umf@ylist,
                               function(x) as.numeric(rowSums(x, na.rm=T)==0)))
  
  #Save formatted covariate frames for use in model frames
  #For predicting with formulas etc
  site_ref <- site_covs
  obs_ref <- obs_covs
  
  #Assign newdata as the covariate frame if it is provided
  if(!is.null(newdata)){
    if(type == "state"){
      site_covs <- newdata
    } else if(type == "det"){
      obs_covs <- newdata
    }
  }
  
  #Design matrices + parameter counts
  #For f/occupancy
  fInd <- c()
  sf_no0 <- stateformulas[!fixed0]
  var_names <- colnames(dmF)[!fixed0]
  dmOcc <- lapply(seq_along(sf_no0),function(i){
    fac_col <- site_ref[, sapply(site_ref, is.factor), drop=FALSE]
    mf <- model.frame(sf_no0[[i]], site_ref)
    xlevs <- lapply(fac_col, levels)
    xlevs <- xlevs[names(xlevs) %in% names(mf)]
    out <- model.matrix(sf_no0[[i]],
                        model.frame(stats::terms(mf), site_covs, na.action=stats::na.pass, xlev=xlevs))
    colnames(out) <- paste('[',var_names[i],'] ',
                           colnames(out), sep='')
    fInd <<- c(fInd,rep(i,ncol(out)))
    out
  })
  fStart <- c(1,1+which(diff(fInd)!=0))
  fStop <- c(fStart[2:length(fStart)]-1,length(fInd))
  occParams <- unlist(lapply(dmOcc,colnames))
  nOP <- length(occParams)
  
  #For detection
  dInd <- c()
  dmDet <- lapply(seq_along(detformulas),function(i){
    fac_col <- obs_ref[, sapply(obs_ref, is.factor), drop=FALSE]
    mf <- model.frame(detformulas[[i]], obs_ref)
    xlevs <- lapply(fac_col, levels)
    xlevs <- xlevs[names(xlevs) %in% names(mf)]
    out <- model.matrix(detformulas[[i]],
                        model.frame(stats::terms(mf), obs_covs, na.action=stats::na.pass, xlev=xlevs))
    colnames(out) <- paste('[',names(umf@ylist)[i],'] ',
                           colnames(out),sep='')
    dInd <<- c(dInd,rep(i,ncol(out)))
    out
  })
  dStart <- c(1,1+which(diff(dInd)!=0)) + nOP
  dStop <- c(dStart[2:length(dStart)]-1,length(dInd)+nOP)
  detParams <- unlist(lapply(dmDet,colnames))
  #nD <- length(detParams)
  
  #Combined
  paramNames <- c(occParams,detParams)
  nP <- length(paramNames)
  
  mget(c("N","S","J","M","nF","fStart","fStop","fixed0","dmF","dmOcc","dmDet",
         "dStart","dStop","y","yStart","yStop","Iy0","z","nOP","nP","paramNames"))
}
