library(psych)
library(polycor)
library(dplyr)
library(geosphere)
library(Matrix)
library(spdep)
library(stargazer)
library(memisc)
library(spatialreg)
#####################################################################################
##     To calculate peer network averages, add network fixed effects,             ###
##      and export to datasets used in regression analysis                        ###
##     Due restricted information (farm coordinates), raw data are not published  ###
#####################################################################################
rm(list=ls())


setwd("")

## raw data
dat <- read.csv('rawdat.csv', stringsAsFactors = F)


na2zero <- function(dat, var){
  dat[,var][is.na(dat[,var]==T)] <- 0
  return(dat[,var])
}


## add calculated vars

add_vars <- function(dat){
  ## binary participation
  dat$part_bin <- 0
  dat$part_bin[dat$part_gen==1] <- 1
  
  # machine - 1: available; 0: not available but possible
  dat$machine <- 0
  dat$machine[dat$av_machinery %in% c(2,3,4,5)] <- 1

  dat$machine_prob <- 0
  dat$machine_prob[dat$av_machinery==11] <- 1
  
  # own personal experience
  dat$exp <- 0
  dat$exp[dat$exp_herbfree %in% c(1,2)] <- 1
  
  # own or others' experience
  dat$exp_herbfree2 <- 0
  dat$exp_herbfree2[dat$exp_herbfree %in% c(1,2, 3, 4)] <- 1
  
  ## wheat land ownership: % of leased land for wheat
  dat$lease <- dat$leasedland
  
  ## farm succession - 1: not established
  dat$succN <- 0
  dat$succN[dat$succession==3] <- 1
  
  ## college degree
  dat$college <- 0
  dat$college[dat$education %in% c(6, 7)] <- 1
  
  dat$soil_conservation <- 0
  dat$soil_conservation[dat$directfed_directseed==1 | dat$directfed_mulchseed==1] <- 1
  
  dat$directcanton_BEpartyears <- na2zero(dat, 'directcanton_BEpartyears')
  dat$directcanton_BEtotyears <- na2zero(dat, 'directcanton_BEtotyears')
  
  dat$suitability_slope <- na2zero(dat, 'suitability_slope')
  
  return(dat)
}


dat <- add_vars(dat)

##########################
## Duplicate coords    ###
##  add random number  ###
########################## 
length(unique(dat$lat))
length(unique(dat$lon))

dat %>% group_by(lat) %>% filter(n()>1) %>%  summarize(n=n())

## find non-unique lat values

duplatdf <- dat %>% group_by(lat) %>% filter(n()>1)
duplat <- duplatdf$lat

dat$lat2 <- dat$lat
for (i in 1:nrow(dat)) {
  if(dat$lat2[i] %in% duplat){
    dat$lat2[i] <- dat$lat2[i]+runif(1, -0.000005, 0.000005)
  #  print(c(dat$lat[i], dat$lat2[i]))
  } 
}


obsbyKt <- table(dat$Kt)
datsub <- dat[dat$Kt %in% names(obsbyKt[obsbyKt>20]), ]


#########################################################################################################################
## remove farms with incorrect lon/lat - e.g. reported canton in Thurgau but long/lat is in Uri, and isolated points  ###
#########################################################################################################################

dmat <- distm(datsub[, c('lon','lat2')])
ptcoords <- SpatialPoints(datsub[, c('lon','lat2')], proj4string=CRS("+init=epsg:4326"))
knb <- knearneigh(ptcoords, 1)
nbb <- knn2nb(knb)
dsts<-nbdists(nbb, as.matrix(datsub[, c('lon','lat2')]), longlat = T)
distvec <- unlist(dsts)

datsub$min_dist2nb <- distvec

## further remove those with canton misreported
mislabeled <- datsub[(datsub$lat2<47.3 & datsub$lon>9) | datsub$lon > 9.45, ]

## remove those with nearest neighbor 10km away and mislabeled
## data for baseline model
datsub2 <- datsub[datsub$min_dist2nb < 10 & !datsub$X %in% mislabeled$X, ]


#############################################################################################
###  Prepare raw data for regression - individual and peer average dependent variable and ###
###   covariates                                                                          ###
###  Two datasets will be exported, respectively include and exclude information on the   ###
###    tendency to consult peers, via (un)commenting respective lines and variables       ###
###    (indicated in function datnw_by_kt below)                                          ###
#############################################################################################


## for robustness check that restricts to those who know about the program 
#datsub2 <-  datsub2[datsub2$part_gen!=0,]


datnw_by_kt <- function(ktdat, nnb){
  dist_kt <- distm(ktdat[, c('lon','lat2')])
  fcoords <- SpatialPoints(ktdat[, c('lon','lat2')], proj4string=CRS("+init=epsg:4326"))
  knb <- knearneigh(fcoords, nnb)
  nbb <- knn2nb(knb)
  dsts<-nbdists(nbb, fcoords)
  idw<-lapply(dsts, function(x)1/(x))
  W1<-nb2listw(nbb, glist=idw, style = "W", zero.policy = TRUE)  
  W <- as(W1, "CsparseMatrix")
  all(W == t(W))
  
  ###################################################################
  nbmat0 <- nb2mat(nbb, style = 'W') # row standardized
  nbmat1 <- nb2mat(nbb, style = 'B') # binary
  distnb <- dist_kt*nbmat1
  distnb_res <- distnb
  distnb_res[distnb_res>10000] <- 0
  nbmat_res <- distnb_res
  nbmat_res[nbmat_res>0] <- 1   
  nn_res <- rowSums(nbmat_res)
  nn_res[nn_res==0] <- 99
  
  wtmat_unw <- nbmat_res/nn_res  ## row standardize on restricted, but not inverse-weighted by distance
  
  y_nb <- as.data.frame(wtmat_unw %*% ktdat$part_bin)
  ###################################################################
  
  wtmat <- as.matrix(W)
  wtmat_info <- ktdat$statement_agdesc_consult * wtmat_unw / 3   ## if not to weight by dist in infomat, use wtmat_unw, otherwise use wtmat
  
  ################
  
  wtmat2 <- wtmat %*% wtmat
  wtmat3 <- wtmat %*% wtmat %*% wtmat
  
  wtmat_info2 <- wtmat_info %*% wtmat_info
  wtmat_info3 <- wtmat_info %*% wtmat_info %*% wtmat_info
  
  # matrix to demean from whole canton
  Ktones <- rep(1, nrow(ktdat)) %*% t(rep(1, nrow(ktdat)))/nrow(ktdat)
  
  # dependent variable and individual independent variables
  yXi <- ktdat[, c('part_bin','shareincome_arable', 'leasedland', 'succN', 'shareincome_offfarm')]
  
  # for models without information exchange between farmers, comment out 'statement_agdesc_consult',
  Xc <- ktdat[, c('statement_agdesc_consult',
                  'age', 'college', 'language', 
                  'agland','share_wheat', 'workforce',  
                  'mean_ertrag_dt_ha_plz','soil_conservation', 
                  'weed_problematic_sum_rel','suitability_grains', 'suitability_slope', 'Tabs_mean', 'Precsum_mean', 
                  'machine','BioFrac')]
  
  ## contextual factors that may have effect via info exchange - should not distinguish
  Xc2 <- ktdat[, c('exp_herbfree2','expyield_decr','expyieldrisk', 'riskpref_plantprotect',  
                   'statement_pos_environ', 'statement_openinnovat' )]
  
  XcW <- wtmat %*% as.matrix(Xc)
  XcW2 <- wtmat2 %*% as.matrix(Xc)
  XcW3 <- wtmat3 %*% as.matrix(Xc)
  
  colnames(XcW) <- paste0(colnames(XcW),'W')
  colnames(XcW2) <- paste0(colnames(XcW),'W2')
  colnames(XcW3) <- paste0(colnames(XcW),'W3')
  XcW <- as.data.frame(XcW)
  
  Xc2W <- wtmat %*% as.matrix(Xc2)
  Xc2W2 <- wtmat2 %*% as.matrix(Xc2)
  Xc2W3 <- wtmat3 %*% as.matrix(Xc2)
  
  colnames(Xc2W) <- paste0(colnames(Xc2W),'W') 
  colnames(Xc2W2) <- paste0(colnames(Xc2W),'W2') 
  colnames(Xc2W3) <- paste0(colnames(Xc2W),'W3')
  Xc2W <- as.data.frame(Xc2W)
  
  ## alternative weights for y - wtmat for models without information exchange between farmers; wtmat_info otherwise
  yW <- wtmat %*% ktdat$part_bin
  
  yW <- wtmat_info %*% ktdat$part_bin
  
  ## add canton for probit
  kt <- ktdat$Kt
  dativ <- cbind(yXi, kt, yW, Xc, XcW, XcW2, XcW3, Xc2, Xc2W, Xc2W2, Xc2W3)
  ## also return wtmat
  rtls <- list(dativ, wtmat, Ktones, y_nb)
  
  return(rtls)
}


## list of df by canton - use datsub2
ktls <- list()
kt <- unique(datsub2$Kt)
for (i in kt) {
  ktdat <- datsub2[datsub2$Kt==i, ]
  ktls <- c(ktls, list(ktdat))
}

## apply fn to each in the list
ktdat1ls <- lapply(ktls, datnw_by_kt, nnb=10)

## separate the returned lists of two dfs into two lists - dativ and wtmat
dativls <- list()
wtmatls <- list()
hmatls <- list()
y_nbls <- list()
for(i in 1:length(ktdat1ls)){
  dativls <- c(dativls, list(ktdat1ls[[i]][[1]]))
  wtmatls <- c(wtmatls, list(ktdat1ls[[i]][[2]]))
  hmatls <- c(hmatls, list(ktdat1ls[[i]][[3]]))
  y_nbls <- c(y_nbls, list(ktdat1ls[[i]][[4]]))
}



## stack in to one df
dat1 <- bind_rows(dativls, .id = "column_label") 

kt_dmean <- bdiag(hmatls)
eye <- diag(nrow = nrow(dat1))
dmeanmat2 <- eye-as.matrix(kt_dmean) ## demean by canton level

ll <- 7  # length of vec yXi + kt plus one, i.e., column of yW, change with how many columns are in Xi
dmeanedx <- dmeanmat2 %*% as.matrix(dat1[, c(2:(ll-1), (ll+1):ncol(dat1))])  ## demean y, Wy, X, WX - exclude kt
dat_demeaned <- as.data.frame(dmeanedx)

## save compiled demeaned dataset with peer averages for upload to research collection

## dataset with social ties defined by only spatial proximity (see notes for function datnw_by_kt)
write.csv(dat_demeaned, "dat_nb_stack_sp.csv", row.names = F)

## dataset with social ties defined by spatial proximity and tendency to consult peers (see notes for function datnw_by_kt)
write.csv(dat_demeaned, "dat_nb_stack_info.csv", row.names = F)
