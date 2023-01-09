#############################################################
##  Models to estimate spillover effects with and without  ##
##   accounting for farmers' tendency to consult peers on  ##
##   agricultural decisions                                ##
#############################################################

## read compiled, demeaned datasets with peer averages from research collection (raw data not published due to restricted information (i.e., spatial coordinates of farms))

## dataset with social ties defined by only spatial proximity - for models (1) and (3)
dat_RC <- read.csv("dat_nb_stack_sp.csv", stringsAsFactors = F)

## dataset with social ties defined by spatial proximity and tendency to consult peers - for models (2) and (4)
dat_RC <- read.csv("dat_nb_stack_info.csv", stringsAsFactors = F)


# first stage regression
ivm <- lm(yW~., data = dat_RC[, c(2:ncol(dat_RC))]) ## exclude kt for de-meaned data
summary(ivm)

yWhat <- ivm$fitted.values  ## predicted demeaned Wy
y <- dat_RC$part_bin # demeaned y

datreg <- cbind(y, yWhat, dat_RC)

# second stage regressions
mm01 <- lm(y ~ yWhat + leasedland + succN + shareincome_offfarm +
             age + college + agland + share_wheat + workforce + mean_ertrag_dt_ha_plz + language +
             soil_conservation + machine + exp_herbfree2 + 
             ageW + collegeW + aglandW + share_wheatW + workforceW + mean_ertrag_dt_ha_plzW + languageW+
             soil_conservationW + machineW +  exp_herbfree2W +
             BioFracW+weed_problematic_sum_relW + suitability_grainsW + suitability_slopeW + Tabs_meanW + Precsum_meanW,  
           data = datreg)
summary(mm01)



mm02 <- lm(y ~ yWhat + leasedland + succN + shareincome_offfarm +
             age + college + agland + share_wheat + workforce + mean_ertrag_dt_ha_plz + language +
             soil_conservation + machine +  exp_herbfree2 +
             expyield_decr + expyieldrisk + riskpref_plantprotect +  statement_pos_environ  +  statement_openinnovat +
             ageW + collegeW + aglandW + share_wheatW + workforceW + mean_ertrag_dt_ha_plzW + languageW+
             soil_conservationW + machineW + exp_herbfree2W +
             expyield_decrW + expyieldriskW + riskpref_plantprotectW + statement_pos_environW + statement_openinnovatW + 
             BioFracW+weed_problematic_sum_relW + suitability_grainsW + suitability_slopeW + Tabs_meanW + Precsum_meanW,  
           data = datreg)
summary(mm02)


# Moran's I test for residuals
dist_all <- distm(datsub2[, c('lon','lat2')])
fcoords_all <- SpatialPoints(datsub2[, c('lon','lat2')], proj4string=CRS("+init=epsg:4326"))
knb_all <- knearneigh(fcoords_all, 10)
nbb_all <- knn2nb(knb_all)
dsts_all <- nbdists(nbb_all, fcoords_all) 

idw_all <- lapply(dsts_all, function(x)1/(x))
W_all <- nb2listw(nbb_all, glist=idw_all, style = "W", zero.policy = TRUE)

lm.morantest(mm02, listw=W_all)

### models reported in Table 2
model1 <- mm01 #
model2 <- mm01 # info 
model3 <- mm02 #
model4 <- mm02 # info

## first stage
ivm1 <- ivm
ivm2 <- ivm  # info

## models removing openness to innovation, or alternative behavioral chars 
model_ex_openinno1 <- mm02
model_ex_openinno2 <- mm02 # info

model_alt_behav1 <- mm02
model_alt_behav2 <- mm02 # info

