#title: Global Distribution, Diversity, and Ecological Niche of Picozoa, a widespread and enigmatic marine protist lineage

############################################################################################################## 

##### Script: Picozoa_distribution_and_diversity.R
##### Created: 11/02/2024
##### Authors: Daniele De Angelis daniele.deangelis_at_uniroma1.it
##### Objective: to compute SMDs of Picozoas in the sunlit Ocean
##### NOTE: for data processing please refer to the manuscript: Global Distribution, Diversity, and Ecological Niche of Picozoa, a widespread and enigmatic marine protist lineage. This script is designed such that each code chunk can be executed independently. T

############################################################################################################## 

#### SDM FOR LIGHT OCEAN ####
#1 load counts per sampling station per OTU
#2 convert count data to presence/absence data ["We used phytoplankton presence data, rather than abundance data, as the former are less sensitive to differences in sampling methods and are more widely available." Righetti et al. 2019]
#3 select mean Bio-oracle variables for surface layer [possibly variables available for both present and future]
#4 for each OTU, use all other sampling stations as "absences", similar to target-group approach
#5 run ensemble modeling 
#6 evaluate models
#7 binarize models based on TSS
#8 create species (OTU) richness maps

### LOAD LIBRARIES ####
library(raster)
library(sdmpredictors)
# library(terra)
library(ggplot2)
library(dplyr)

##########################################################################################
####---- LOAD RASTER VARIABLES AND CREATE RASTER STACK SEPARATED BY DEPTH LEVEL ----######


MLD<-raster::raster("Ecological Niche/env_variables/M02.tif")
names(MLD)<-c("mld")

# Open environmental variables (surface layer)
(names <- list.files(path = pathpath, 
                     pattern = ".tif", full.names = TRUE))

index <- grep("EPI",ignore.case = T,names) #select only mean variables
stacked_epi <- stack(names[index])
names(stacked_epi) <- gsub(".tif","", gsub(pathpath,"", names[index]))
names(stacked_epi) <- gsub("_withNA","", names(stacked_epi) ) 
names(stacked_epi)
dim(stacked_epi)

index <- grep("MESO",ignore.case = T,names) #select only mean variables
stacked_meso <- stack(names[index])
names(stacked_meso) <- gsub(".tif","", gsub(pathpath,"", names[index]))
names(stacked_meso) <- gsub("_withNA","", names(stacked_meso) ) 
names(stacked_meso)
dim(stacked_meso)

index <- grep("BATHY",ignore.case = T,names) #select only mean variables
stacked_bathy <- stack(names[index])
names(stacked_bathy) <- gsub(".tif","", gsub(pathpath,"", names[index]))
names(stacked_bathy) <- gsub("_withNA","", names(stacked_bathy)) 
names(stacked_bathy)
dim(stacked_bathy)
# 
# plot(stacked_epi)
# plot(stacked_meso)
# plot(stacked_bathy)

Nstar_epi<-(stacked_epi[["n_epi"]]-16*stacked_epi[["p_epi"]])
Nstar_meso<-(stacked_meso[["n_meso"]]-16*stacked_meso[["p_meso"]])
Nstar_bathy<-(stacked_bathy[["n_bathy"]]-16*stacked_bathy[["p_bathy"]])

names(Nstar_epi)<-"Nstar_epi"
names(Nstar_meso)<-"Nstar_meso"
names(Nstar_bathy)<-"Nstar_bathy"

stacked_epi<-raster::stack(MLD,stacked_epi, Nstar_epi)
stacked_meso<-raster::stack(MLD, stacked_meso, Nstar_meso)
stacked_bathy<-raster::stack(MLD, stacked_bathy, Nstar_bathy)

usdm::vif(data.frame(as.matrix(stacked_epi)))
usdm::vif(data.frame(as.matrix(stacked_meso)))
usdm::vif(data.frame(as.matrix(stacked_bathy)))
# 
# > usdm::vif(stacked_epi)
# Variables         VIF
# 1        mld    2.359181
# 2      A_epi  238.063864
# 3      C_epi 4481.436112
# 4    chl_epi    1.670151
# 5      I_epi  286.897494
# 6      n_epi   53.291392
# 7      o_epi  435.866493
# 8  O2sat_epi  159.693295
# 9      p_epi   53.248729
# 10     s_epi  218.810633
# 11    si_epi    9.702234
# 12     t_epi 4100.985783
# 13   ugo_epi    1.148185
# 14   vgo_epi    1.010233
# > usdm::vif(stacked_meso)
# Variables         VIF
# 1         mld    1.911089
# 2      A_meso 1367.805581
# 3      C_meso 3130.951190
# 4    chl_meso    3.817895
# 5      I_meso   24.120536
# 6      n_meso   28.152927
# 7      o_meso  968.051816
# 8  O2sat_meso  539.343940
# 9      p_meso   44.710589
# 10     s_meso  112.155303
# 11    si_meso   10.260671
# 12     t_meso 2907.929870
# 13   ugo_meso    1.257206
# 14   vgo_meso    1.018514
# > usdm::vif(stacked_bathy)
# Variables         VIF
# 1          mld    1.602082
# 2      A_bathy 5235.747066
# 3      C_bathy 3672.092848
# 4    chl_bathy    2.612742
# 5      I_bathy   52.738340
# 6      n_bathy   81.231316
# 7      o_bathy 2752.671188
# 8  O2sat_bathy 2051.807500
# 9      p_bathy  111.261465
# 10     s_bathy  223.835306
# 11    si_bathy   35.338344
# 12     t_bathy 3310.197968
# 13   ugo_bathy    1.081424
# 14   vgo_bathy    1.018837

usdm::vif(data.frame(as.matrix(stacked_epi[[which(!names(stacked_epi) %in% c("n_epi","p_epi", "A_epi", "C_epi", "I_epi", "O2sat_epi","t_epi"))]])))
usdm::vif(data.frame(as.matrix(stacked_meso[[which(!names(stacked_meso) %in% c("n_meso","p_meso", "A_meso", "C_meso", "I_meso", "O2sat_meso", "t_meso"))]])))
usdm::vif(data.frame(as.matrix(stacked_bathy[[which(!names(stacked_bathy) %in% c("n_bathy","p_bathy", "A_bathy", "C_bathy", "I_bathy", "O2sat_bathy", "t_bathy"))]])))


##########################################################################################
####---- OPEN ALL SAMPLING STATION DATA WITH OTU COUNTS AND ABSENCES ---------------######

ABS <- read.csv("picozoa_absences_SDM.csv") #sampling stations with 0 picozoa
length(which(is.na(ABS$depth)))/nrow(ABS)*100


MM <- read.csv("picozoa_dataSDM.csv")

range(MM$depth, na.rm=T)
hist(MM$depth, xlab="Depth (m)", main="ASV frequency along depth levels")
length(which(is.na(MM$water_colum)))/nrow(MM)*100



names(ABS)
names(MM)

library(plyr)

DF <- rbind.fill(MM,ABS) 
dim(DF)
tail(DF,10)


any(is.na(DF$water_colum))

DF$water_colum[DF$water_colum=="SURF"] <- "EPI"


png("Sampling_Stations_PerDepth.png")
barplot(table(DF$water_colum), main="N of sampling stations per depth layer")
dev.off()

write.csv(DF, "MASTER_DF_picozoi_01APR24.csv")

# convert ASV counts to presence-absence
ASV_MATRIX <- DF[,c(3:176)]
ASV_MATRIX[is.na(ASV_MATRIX)] <- 0
ASV_MATRIX[ASV_MATRIX>0] <- 1
table(unlist(c(ASV_MATRIX)))

# convert ASV counts to presence-absence using a cut-off >= 2
# ASV_MATRIX <- DF[,c(3:176)]
# ASV_MATRIX[ASV_MATRIX<=2] <- 0
# ASV_MATRIX[ASV_MATRIX>1] <- 1
# table(unlist(c(ASV_MATRIX)))
# #data loss
# (1-(11402/20749))*100

# create spatial vector
DF_bin<- cbind(DF[,c(1:2)], ASV_MATRIX, DF[,c(177:ncol(DF))])
names(DF_bin)
DF_bin <- DF_bin[!is.na(DF_bin$long),] #6 sampling stations have no coords
DF_bin <- DF_bin[!is.na(DF_bin$water_colum),]  #80 sampling stations have no indication of the depth
DF_bin <- SpatialPointsDataFrame(data.frame(DF_bin$long,DF_bin$lat), data=DF_bin)


# extract variables
(DF_epi <- DF_bin[DF_bin@data$water_colum == "EPI",])
(DF_meso <- DF_bin[DF_bin@data$water_colum == "MESO",])
(DF_bathy <- DF_bin[DF_bin@data$water_colum == "BATHY",])
ex.epi <- cbind(DF_epi@data, terra::extract(stacked_epi, DF_epi))
ex.meso <- cbind(DF_meso@data, terra::extract(stacked_meso, DF_meso))
ex.bathy <- cbind(DF_bathy@data, terra::extract(stacked_bathy, DF_bathy))


# replace names and Rbind data.frames
(names(ex.epi) <- gsub("epi","",names(ex.epi)))
(names(ex.meso) <- gsub("meso","",names(ex.meso)))
(names(ex.bathy) <- gsub("bathy","",names(ex.bathy)))

rb1 <- rbind(ex.epi, ex.meso, ex.bathy) 

#MOST ABUNDANT
info <- xlsx::read.xlsx("pOTU_tax_information_Jan24.xlsx",1)
info <- info[,c("Final_Name","Lat_distribution")]
selected_otus <- info[info$Lat_distribution!="LA","Final_Name"]

df.<-rb1[,selected_otus]
df.<-lapply(df.[], as.character)
df.<-lapply(df.[], as.numeric)
df.<-data.frame(df., "water_colum"=rb1$water_colum)

(nepi<-(colSums(df.[df.$water_colum=="EPI",-ncol(df.)])))
(nmeso<-(colSums(df.[df.$water_colum=="MESO",-ncol(df.)])))
(nbathy<-(colSums(df.[df.$water_colum=="BATHY",-ncol(df.)])))

N_ASV_PER_DEPTH<-data.frame(nepi, nmeso, nbathy)
# write.csv(N_ASV_PER_DEPTH,"N_ASV_PER_DEPTH.csv")


###----------------TRY MODELS -------------####
#NO MaxEnt, as we have presence/absence

library(MASS)

write.csv(rb1, "Picozoa_DF_with_WOA_Vars_01APR24.csv")

rb1$water_colum <- factor(rb1$water_colum, levels = c("EPI", "MESO", "BATHY"))

#CHECK NUMBER OF OBSERVATIONS PER ASV
(num_obs<-rb1 %>%
    group_by(water_colum) %>%
    dplyr::summarise(across(all_of(selected_otus), sum))%>%data.frame())


(summary1<-t(data.frame(lapply(num_obs[2:ncol(num_obs)], as.numeric))))
length(which(sort(rowSums(summary1))>20)) #52 OTUs have >20 detections across depth levels
summary2<-data.frame(summary1)
names(summary2)<-c("Epi", "Meso", "Bathy")
summary2$ASV<-rownames(summary2)

head(dd<-reshape::melt(summary2))
dd$ASV<-as.factor(dd$ASV)
names(dd)<-c("ASV", "Depth", "N_obs")

#######################################################################################
#EXCLUDE REPEATED SAMPLES: take only one row for each combination of depth layer ######
#######################################################################################
set.seed(999)
ix<-paste(rb1$long,rb1$lat,rb1$water_colum)
head(sort(table(ix), decreasing=T))

rb1$tempID<-ix
rb2_ <- rb1 %>% group_by(tempID) %>% slice_sample(n=1)
rb2_$water_colum <- factor(rb2_$water_colum, levels = c("EPI", "MESO", "BATHY"))

################################################################################################
#CHECK NUMBER OF OBSERVATIONS PER ASV ################################################
################################################################################################
(num_obs<-rb2_ %>%
   group_by(water_colum) %>%
   dplyr::summarise(across(all_of(selected_otus), sum))%>%data.frame())


(summary1<-t(data.frame(lapply(num_obs[2:ncol(num_obs)], as.numeric))))
length(which(sort(rowSums(summary1))>20)) #61 OTUs have >20 detections across depth levels


summary2<-data.frame(summary1)
names(summary2)<-c("Epi", "Meso", "Bathy")
summary2$ASV<-rownames(summary2)
summary3<-summary2[summary2$ASV %in% names(which(sort(rowSums(summary1))>20)),]

head(dd<-reshape::melt(summary3))
dd$ASV<-as.factor(dd$ASV)
names(dd)<-c("pOTU", "Depth", "N_obs")

(prevalence<-round(summary1/(nrow(rb2_)-summary1),2))
prevalence<-prevalence[which(row.names(prevalence) %in% selected_otus),]
length(which(prevalence[,1]>1))
length(which(prevalence[,1]<0.01))

################################################
############# PLOT ASV NUMBER ######################
gridExtra::grid.arrange(gg1,gg2, ncol=1)
################################################


################################################
###INCLUDE ONLY OTU with MORE THAN 20 OBSERVATIONS###
selected_otus<-names(which(sort(rowSums(summary1))>20))
################################################



################################################
################################################
#------------- SDM ---------------##############
#################################################
################################################
################################################


library(h2o)
h2o.init()
rb2_[3:176] <- lapply(rb2_[3:176], factor)
rb2_$water_colum <- factor(rb2_$water_colum, levels=c("EPI","MESO","BATHY"))
h2o_df<- as.h2o(rb2_)
table(as.vector(h2o_df$water_colum))

# susbet variables to keep based on previous step (VIF check)

vars_epi<-(stacked_epi[[which(!names(stacked_epi) %in% c("n_epi","p_epi", "A_epi", "C_epi", "O2sat_epi"))]])
vars_meso<-(stacked_meso[[which(!names(stacked_meso) %in% c("n_meso","p_meso", "A_meso", "C_meso", "O2sat_meso"))]])
vars_bathy<-(stacked_bathy[[which(!names(stacked_bathy) %in% c("n_bathy","p_bathy", "A_bathy", "C_bathy", "O2sat_bathy"))]])


usdm::vif(data.frame(rb2_[,c("mld", "o_", "s_", "n_","p_", "A_", "C_", "I_", "O2sat_","t_", "Nstar_")]))
usdm::vif(data.frame(rb2_[,c("mld", "o_", "s_", "A_", "C_", "I_", "O2sat_","t_", "Nstar_")]))
usdm::vif(data.frame(rb2_[,c("mld", "o_", "s_", "A_", "I_","t_", "Nstar_")]))
usdm::vif(data.frame(rb2_[,c("mld", "o_", "s_", "I_","t_", "Nstar_")]))


library(GGally)
ggcorr(data.frame(rb2_[,c("mld", "o_", "s_", "I_", "t_", "Nstar_")]),
       digits = 2, label =T )



#######################################################################################
# MESS ANALYSIS ######################################################################## 
#######################################################################################
library(dismo)
mess.epi<-mess(vars_epi, rb1[rb1$water_colum == "EPI", gsub("epi","",names(vars_epi))], full=FALSE, filename="Mess_epi_woa_01APR24.tif", overwrite=T)
mess.meso<-mess(vars_meso, rb1[rb1$water_colum == "MESO", gsub("meso","",names(vars_meso))], full=FALSE, filename="Mess_meso_woa_01APR24.tif",  overwrite=T)
mess.bathy<-mess(vars_bathy, rb1[rb1$water_colum == "BATHY", gsub("bathy","",names(vars_bathy))], full=FALSE, filename="Mess_bathy_woa_01APR24.tif",  overwrite=T)

par(mfrow=c(3,2))
plot(mess.epi, main="MESS")
hist(values(mess.epi),1000, main="Histogram of mess values")
plot(mess.meso)
hist(values(mess.meso),1000, main=" ")
plot(mess.bathy)
hist(values(mess.bathy),1000, main=" ")
# 




dddfff <- data.frame(as.data.frame(vars_epi), water_colum="EPI")
(names(dddfff)<-gsub("epi","",names(dddfff)))
h2o.nd.epi <- as.h2o(dddfff) #newdata1

dddfff2 <- data.frame(as.data.frame(vars_meso), water_colum="MESO")
(names(dddfff2)<-gsub("meso","",names(dddfff2)))
h2o.nd.meso <- as.h2o(dddfff2) #newdata2

dddfff3 <- data.frame(as.data.frame(vars_bathy), water_colum="BATHY")
(names(dddfff3)<-gsub("bathy","",names(dddfff3)))
h2o.nd.bathy <- as.h2o(dddfff3) #newdata3


PREDICTOR_NAMES <-   c(gsub("epi","",names(vars_epi)), "water_colum")

TSS_EPI<-c()
TSS_MESO<-c()
TSS_BATHY<-c()
AUC_EPI<-c()
AUC_MESO<-c()
AUC_BATHY<-c()
AUCC_list<-list()
start<-Sys.time()
AUC_df_final<-c()

set.seed(999)
setwd("C:/workD/ddeangelis/Picozoa/results_EM_21apr24_varimp")
for(otu in 1:length(selected_otus)){
  
  glm_model<-NULL
  rf_model<-NULL
  ann_model<-NULL
  gam_model<-NULL
  
  j <- selected_otus[otu]
  print(j)
  
  ### --- RUN GLM MODEL --- ####
  glm_model <- h2o.glm(x = PREDICTOR_NAMES,
                       y = j,
                       training_frame = h2o_df,
                       model_id = "GLM_Model",
                       family="binomial",
                       nfolds = 5,
                       standardize = F)
  
  
  
  ### --- RUN RANDOM FOREST MODEL --- ####
  rf_model <- h2o.randomForest(y = j,
                               x = PREDICTOR_NAMES,
                               training_frame = h2o_df,
                               ntrees = 750,
                               nfolds = 5,
                               model_id = "RF_model",
                               categorical_encoding = "OneHotExplicit")
  
  
  ### --- RUN ANN MODEL --- ####
  ann_model <- h2o.deeplearning(y = j,
                                x = PREDICTOR_NAMES,
                                training_frame = h2o_df,
                                nfolds = 5,
                                model_id = "ANN_model",
                                standardize = F)
  
  ### --- RUN BOOSTED REGRESSION TREE --- ####
  gbm_model <- h2o.gbm(y = j,
                       x = PREDICTOR_NAMES,
                       training_frame = h2o_df,
                       nfolds = 5,
                       model_id = "GBM_model",
                       categorical_encoding = "OneHotExplicit")
  
  
  MODELS <- list("GLM"=glm_model, "RF"=rf_model, "ANN"=ann_model, "GBM"=gbm_model)
  MPERFORM <- lapply(MODELS, FUN = h2o.performance)
  saveRDS(MPERFORM, paste0("ModelPerformance_",j,".rds"))
  saveRDS(MODELS, paste0("Models_",j,".rds"))
  
  
  ### --- SUMMARIZE MODEL EVALUATION METRICS --- ####
  AUCC <- lapply(MODELS, FUN = h2o.auc)
  AUC_thr <- 0.7
  
  AUCC_list[[j]]<-AUCC
  saveRDS(AUCC_list, "AUC_list2.rds")
  
  
  Modeval_Summary <- data.frame(ALGO = names(MODELS), 
                                AUC = unlist(AUCC))
  
  
  models_above_auc_threshold_T <- Modeval_Summary[Modeval_Summary$AUC>=AUC_thr,]
  models_above_auc_threshold <- models_above_auc_threshold_T$ALGO
  n_good_models <- length(models_above_auc_threshold)
  GOOD_MODELS <- MODELS[names(MODELS) %in% models_above_auc_threshold]
  
  if(n_good_models>0){
    
    # ### --- VARIABLE IMPORTANCE --- #####

    varIMP <- list("GLM"= h2o.varimp(glm_model)[,c(1,4)],
                   "RF"=h2o.varimp(rf_model)[,c(1,4)],
                   "ANN"=h2o.varimp(ann_model)[,c(1,4)],
                   "GBM"=h2o.varimp(gbm_model)[,c(1,4)])
    
    
    
    variable_importance<-varIMP[which(names(varIMP) %in% names(GOOD_MODELS))]
    saveRDS(variable_importance, paste0("VarImp_",j,".rds"))
    
    if(!n_good_models>0){
      paste0("Ouch! All models got an AUC score <", AUC_thr, " - the analysis stops here")}else if(n_good_models==1){
        paste0("Humm... Only ", models_above_auc_threshold, " model got an AUC score >",AUC_thr, " - the final model corresponds to this model. Check 'Model Results' window.")}else{
          paste0(n_good_models, " models have an AUC score >", AUC_thr, " and were used in the ensemble model. Check 'Model Results' window.")}
    
    
    ### --- PREDICT RASTERS -----#####
    
    ##---- Epipelagic -----###
    predictions <- lapply(GOOD_MODELS, FUN = h2o.predict, h2o.nd.epi) #predict models
    
    m_layer <- vars_epi[[1]]
    l2 <- lapply(predictions, "[[", 3) # take third column (with predicted probability of presence)
    vals_p1 <- lapply(l2, function (x) ifelse(is.na(values(m_layer)), NA, as.numeric(as.vector(x)))) # for each algorithm above threshold create vectors with either NA (if raster has NA in that cel, or predicted probability of presence)
    range01 <- function(x){(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))}
    rs <- lapply(vals_p1, FUN = range01) #scale values between 0 and 1
    
    RStack_epi <- stack()
    for(th in 1:length(models_above_auc_threshold)){ #create raster stack
      temp <- m_layer
      values(temp) <- as.matrix(rs[[models_above_auc_threshold[th]]])
      RStack_epi <- stack(RStack_epi, temp)
    }
    names(RStack_epi) <- models_above_auc_threshold
    RStack_epi
    ##plot(RStack_epi)
    
    
    ##---- Mesopelagic -----###
    predictions <- lapply(GOOD_MODELS, FUN = h2o.predict, h2o.nd.meso) #predict models
    
    m_layer <- vars_meso[[1]]
    l2 <- lapply(predictions, "[[", 3) # take third column (with predicted probability of presence)
    vals_p1 <- lapply(l2, function (x) ifelse(is.na(values(m_layer)), NA, as.numeric(as.vector(x)))) # for each algorithm above threshold create vectors with either NA (if raster has NA in that cel, or predicted probability of presence)
    rs <- lapply(vals_p1, FUN = range01) #scale values between 0 and 1
    
    RStack_meso <- stack()
    for(th2 in 1:length(models_above_auc_threshold)){ #create raster stack
      temp <- m_layer
      values(temp) <- as.matrix(rs[[models_above_auc_threshold[th2]]])
      RStack_meso <- stack(RStack_meso, temp)
    }
    names(RStack_meso) <- models_above_auc_threshold
    RStack_meso
    #plot(RStack_meso)  
    
    
    ##---- Bathypelagic -----###
    predictions <- lapply(GOOD_MODELS, FUN = h2o.predict, h2o.nd.bathy) #predict models
    
    m_layer <- vars_bathy[[1]]
    l2 <- lapply(predictions, "[[", 3) # take third column (with predicted probability of presence)
    vals_p1 <- lapply(l2, function (x) ifelse(is.na(values(m_layer)), NA, as.numeric(as.vector(x)))) # for each algorithm above threshold create vectors with either NA (if raster has NA in that cel, or predicted probability of presence)
    rs <- lapply(vals_p1, FUN = range01) #scale values between 0 and 1
    
    RStack_bathy <- stack()
    for(th3 in 1:length(models_above_auc_threshold)){ #create raster stack
      temp <- m_layer
      values(temp) <- as.matrix(rs[[models_above_auc_threshold[th3]]])
      RStack_bathy <- stack(RStack_bathy, temp)
    }
    names(RStack_bathy) <- models_above_auc_threshold
    RStack_bathy
   
    
    
    PLOT_LIST <- lapply(GOOD_MODELS, FUN=h2o.partialPlot, data = h2o_df, #set this to MODELS to show response curves for all models independently to AUC values
                        cols = PREDICTOR_NAMES,
                        nbins = 50,
                        plot = F)
    
    saveRDS(PLOT_LIST, paste0(j, "_PPlot", "_epi.rds"))
    
    meanEPI = mean(RStack_epi, na.rm=T)
    meanMESO = mean(RStack_meso, na.rm=T)
    meanBATHY = mean(RStack_bathy, na.rm=T)
    
    writeRaster(meanEPI, paste0(j, "_meanModel_", "epi.tif"), overwrite=T)
    # writeRaster(meanMESO, paste0(j, "_meanModel_", "meso.tif"), overwrite=T)
    # writeRaster(meanBATHY, paste0(j, "_meanModel_", "bathy.tif"), overwrite=T)
    
    if(n_good_models>1){
      #min
      minEPI = min(RStack_epi, na.rm=F)
      minMESO = min(RStack_meso, na.rm=F)
      minBATHY = min(RStack_bathy, na.rm=F)
      
      writeRaster(minEPI, paste0(j, "_minModel_", "epi.tif"), overwrite=T)
      # writeRaster(minMESO, paste0(j, "_minModel_", "meso.tif"), overwrite=T)
      # writeRaster(minBATHY, paste0(j, "_minModel_", "bathy.tif"), overwrite=T)
      
      #max
      maxEPI = max(RStack_epi, na.rm=F)
      # maxMESO = max(RStack_meso, na.rm=F)
      # maxBATHY = max(RStack_bathy, na.rm=F)
      
      writeRaster(maxEPI, paste0(j, "_maxModel_", "epi.tif"), overwrite=T)
      # writeRaster(maxMESO, paste0(j, "_maxModel_", "meso.tif"), overwrite=T)
      # writeRaster(maxBATHY, paste0(j, "_maxModel_", "bathy.tif"), overwrite=T)
       
      #sd
      sdEPI = calc(RStack_epi, sd, na.rm=F)
      # sdMESO = calc(RStack_meso, sd, na.rm=F)
      # sdBATHY = calc(RStack_bathy, sd, na.rm=F)
       
      writeRaster(sdEPI, paste0(j, "_sdModel_", "epi.tif"), overwrite=T)
      # writeRaster(sdMESO, paste0(j, "_sdModel_", "meso.tif"), overwrite=T)
      # writeRaster(sdBATHY, paste0(j, "_sdModel_", "bathy.tif"), overwrite=T)
      
      #cv
      cvEPI = raster::cv(RStack_epi, na.rm=F)
      # cvMESO = raster::cv(RStack_meso, na.rm=F)
      # cvBATHY = raster::cv(RStack_bathy, na.rm=F)
      
      writeRaster(cvEPI, paste0(j, "_CV_", "epi.tif"), overwrite=T)
      # writeRaster(cvMESO, paste0(j, "_CV_", "meso.tif"), overwrite=T)
      # writeRaster(cvBATHY, paste0(j, "_CV_", "bathy.tif"), overwrite=T)
    }
    
    # png(paste0(j,".png"), width = 500, height = 800, units = "px", pointsize = 12)
    # par(mfrow=c(3,1))
    # image(meanEPI, main=paste(j), ylab="Epipelagic")
    # image(meanMESO, main=paste(j), ylab="Mesopelagic")
    # image(meanBATHY, main=paste(j), ylab="Bathypelagic")
    # dev.off()
    
    
    # BINARIZE MODELS ######
    epi
    p_EPI=terra::extract(meanEPI, DF_bin[DF_bin@data[DF_bin$water_colum=="EPI",j]==1,]@coords)
    a_EPI=terra::extract(meanEPI, DF_bin[DF_bin@data[DF_bin$water_colum=="EPI",j]==0,]@coords)

    if((!is.null(p_EPI))&(length(na.omit(p_EPI)))>0){
      eval_EPI<- dismo::evaluate(p_EPI,a_EPI)
      TSS_EPI[j]<-dismo::threshold(eval_EPI, "spec_sens")
      meanEPIbin<-meanEPI
      values(meanEPIbin)<-ifelse(values(meanEPIbin)>=TSS_EPI[j], 1, 0)
      writeRaster(meanEPIbin, paste0(j, "_TSS_BinModel_", "epi.tif"), overwrite=T)
    }else{
      TSS_EPI[j]<-NA
      print(paste0("No epipelagic presences for otu ", j))}
     
    # p_MESO=terra::extract(meanMESO, DF_bin[DF_bin@data[DF_bin$water_colum=="MESO",j]==1,]@coords)
    # a_MESO=terra::extract(meanMESO, DF_bin[DF_bin@data[DF_bin$water_colum=="MESO",j]==0,]@coords)
    # 
    # if((!is.null(p_MESO))&(length(na.omit(p_MESO)))>0){
    #   eval_MESO<- dismo::evaluate(p_MESO,a_MESO)
    #   TSS_MESO[j]<-dismo::threshold(eval_MESO, "spec_sens")
    #   meanMESObin<-meanMESO
    #   values(meanMESObin)<-ifelse(values(meanMESObin)>=TSS_MESO[j], 1, 0)
    #   writeRaster(meanMESObin, paste0(j, "_TSS_BinModel_", "meso.tif"), overwrite=T)
    # }else{
    #   TSS_MESO[j]<-NA
    #   print(paste0("No mesopelagic presences for otu ", j))}
    # 
    # #bathy
    # p_BATHY=terra::extract(meanBATHY, DF_bin[DF_bin@data[DF_bin$water_colum=="BATHY",j]==1,]@coords)
    # a_BATHY=terra::extract(meanBATHY, DF_bin[DF_bin@data[DF_bin$water_colum=="BATHY",j]==0,]@coords)
    # 
    # if((!is.null(p_BATHY))&(length(na.omit(p_BATHY)))>0){
    #   eval_BATHY<- dismo::evaluate(p_BATHY,a_BATHY)
    #   TSS_BATHY[j]<-dismo::threshold(eval_BATHY, "spec_sens")
    #   meanBATHYbin<-meanBATHY
    #   values(meanBATHYbin)<-ifelse(values(meanBATHYbin)>=TSS_BATHY[j], 1, 0)
    #   writeRaster(meanBATHYbin, paste0(j, "_TSS_BinModel_", "bathy.tif"), overwrite=T)
    # }else{
    #   TSS_BATHY[j]<-NA
    #   print(paste0("No bathypelagic presences for otu ", j))}
    
  }else{print(paste0("The following OTU did not give models with AUC>threshold: ", j))}
  
  # AUC_df<-data.frame(AUCC)
  # rownames(AUC_df)<-paste(j)
  # AUC_df_final<-rbind(AUC_df_final,AUC_df)
  # write.csv(AUC_df_final,"AUC_3foldCV.csv")
  
  end<-Sys.time() 
}
end-start 
tss_tab<-data.frame(TSS_EPI, TSS_MESO, TSS_BATHY)
write.csv(tss_tab, "TSS_TABLE.csv")




