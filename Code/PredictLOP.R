# General information ------------------------------------------

# Use MCMC samples (fit) to generate
## prediction (Deviance)
## in principle the same as calculating Deviance, just for test data in all folds
## difference: individual posterior parameter samples are generated from group-level
## distributions.
## We parallized across folds in the routines below since the approximation by
## sampling can be prohibitively slow for some data sets, models and model types

# for
## models without items effects (functions in LOPThetaInd.R)
## models without item effects, random crit (criteria are given for individuals as deviations from grand mean) (functions in LOPThetaInd_gam.R)
## models with item effects, random crit (functions in LOPTheta.R)

# (LOPTheta files have some redundancies, but it was the easiest way to keep it clean for modeltypes and avoid crossover)

# in the routine below, the correct procedure for model type is indicated by using or changing column: itemInd in data file
# original hardcoding has itemInd = "Item" for all data sets that have item information and "Ind" otherwise
## itemInd == "Ind" : No item effects
## itemInd == "Indgam" : No item effects, random crit
## itemInd == "Item" : Item effects (only for data sets that actually contain item information)


# Note: due to damage on the servers, we do not have the posterior samples / fits for the
# training fits for any KFCV or LOP fits

#library("rstan")
#library("tidybayes")
#library("brms")
library("cmdstanr")
library("Matrix")
library("tidyverse")
#library("bayestestR")
library("foreach")
library("loo")

options(contrasts=c('contr.orthonorm', 'contr.poly'))


rds_file <- list.files("../ProcessedData/", full.names = TRUE)
rds_file <- rds_file[grepl(".rds",rds_file)]

# put all data in a single object
alldata <- purrr::map(c(rds_file), readRDS)

# list data/experiment names
listexp <- alldata %>% purrr::map(. %>% ungroup() %>% distinct(exp)) %>% bind_rows() %>% .$exp
listexpitem <- alldata %>% purrr::map(. %>% filter(itemInd=="Item") %>% ungroup() %>% distinct(exp)) %>% bind_rows() %>% .$exp

source("MakeStanDataModel.R")

# No item effects -------------------------------------------------

# prediction functions for this kind of model in here
source("LOPThetaInd.R")

# for the DP-SDT model with response mapping for recollection (DPRMSDT/DPRM2SDT),
# we fitted two versions. We initially fitted DPRMSDT where recollection is spread
# across all 'old' confidence ratings. For some data sets this resulted in a
# high percentage of divergent transitions. For those data sets, we additionally
# fitted a reduced DPRMSDT model where recollection was only spread across the
# top two 'old' confidence ratings. In the manuscript, both versions are represented
# by DPRMSDT

# dprmexp <- DPRMSDT fitted
# dprm2exp <- DPRM2SDT fitted

dprm2exp <- c("AP2007_e1","AP2007_e3","BKSR2013_e1","NFR2013_e1","NFR2013_e3","GKH1999_e1",
              "JCD2012_e1a","JCM2019_e1","KAWY2013_e4","LBA2019_e1","LM2020_e1","TR2017_e1",
              "WKS2020_e1","KL2012_e2","KL2012_e3", "BTL2013_e1","CSR2015_e1","DR2012_e1b",
              "DW2020_e1","KK2015_e1","KUO2017_e2","MKG2018_e2","OZH2010_e1","PCM2006_e1",
              "PRM2010_e1","RSP2012_e1","SB2020_e1","SBT2018_e1","SCR2019_e1","JW2019_e1",
              "SD2014_e2","SHJ2005_e1","TMP2014_e1","WHD2018_e1", "LP2013_e2")

dprmexp <- c("AP2007_e2","APH2016_e1","BSG2014_e1","D2007_e1a","D2007_e1b","FBH2013_e1",
             "FGR2019_e1","FO2016_e2a", "FO2016_e2b","FO2016_e3","GKH1999_e2","GKH1999_e3",
             "GKH1999_e4","HDM2006_e1","HDM2006_e2","HUW2015_e1","JCD2012_e1b","JWH2009_e1",
             "KAWY2013_e2","KAWY2013_e3","KFH2013_e1","KL2012_e1","KL2012_e4","KY2010_e1",
             "KY2011_e1","KY2016_e1","LBA2019_e2", "LP2006_e2","LP2013_e1","NFR2013_e2",
             "OBD2017_e1", "QGM2021_e1","RS2009_e1","RS2009_e2","RSM2009_e1","SB2020_e2",
             "SB2020_e3","SD2004_e2","SD2014_e1", "USS2015_e1","WKH2020_e1",
             "WMS2012_e1","ZMD2011_e1","ZOL2021_e1")



for(i in listexp){

  #i <- "FGR2019_e1"
  model <- "EVSDT"

  #fitfiles <- list.files(paste0("../Fits/NoItemEffects/KFCV/",i),full.names = T)

  exp_pred <- alldata[[which(listexp == i)]] %>%
    mutate(rating = as.numeric(as.character(rating))) %>%
    mutate(itemInd = "Ind")

  #print(unique(exp_pred$exp))

  mcoptions <- list(preschedule=TRUE, set.seed=FALSE)
  doParallel::registerDoParallel(cores=10)
  print(expstring)
  allfolds <- foreach(fold = 1:10

                      ,.options.multicore=mcoptions) %dopar% {

                        prepmodel_pred <- selectModel(exp_pred,model)

                        testf <- makeFreqDataInd(exp_pred,model) %>% .$LOPFolds == fold

                        tmpdat_pred <- prepLOP_Ind(prepmodel_pred$tmpdat,testf,model)

                        # predict out-of-sample probabilities with samples from group-level
                        pred <- predictionLOP_Ind(tmpdat_pred,fit_f,model)

                        # call the fit file for the training data for the specific fold, model, experiment, modeltype
                        fit_f <- readRDS(paste0("../Fits/NoItemEffects/LOP/",expstring,"/fit_",model,"_LOP_",unique(exp_pred$exp),"_LOPFold",fold,".rds"))


                        # Determine diagnostics
                        divergent <- sum(fit_f$sampler_diagnostics()[,,2])/length(fit_f$sampler_diagnostics()[,,2])

                        fullres <-  data.frame(matrix(pred %>% unlist()*-2,nrow=6000,byrow=T)) %>%
                          mutate(DivChainsPercent = divergent) %>%
                          mutate(model = model,
                                 type = "LOP",
                                 Fold = fold,
                                 experiment = expstring,
                                 manipulation = unique(exp_pred$manipulation),
                                 itemInd = unique(exp_pred$itemInd))

                        return(fullres)

                      }

  saveRDS(allfolds, paste0("../Predictions/NoItemEffects/LOP/",model,"/DevT1_", model,"_LOP_",unique(exp_pred$exp),".rds"))


}

# No Item Effects, Random Crit -------------------------------------------------

# prediction functions for this kind of model in here
source("LOPThetaInd_gam.R")

selectdata <- c("APH2016_e1", "D2007_e1b","DW2020_e1", "FGR2019_e1", "FO2016_e2a",
                "KFH2013_e1", "KY2011_e1", "SB2020_e2", "SD2014_e1", "ZOL2021_e1")


for(i in selectdata){

  #i <- "FGR2019_e1"
  model <- "EVgamSDT"

  #fitfiles <- list.files(paste0("../Fits/NoItemEffects/KFCV/",i),full.names = T)

  exp_pred <- alldata[[which(listexp == i)]] %>%
    mutate(rating = as.numeric(as.character(rating))) %>%
    mutate(itemInd = "Indgam")

  #print(unique(exp_pred$exp))

  mcoptions <- list(preschedule=TRUE, set.seed=FALSE)
  doParallel::registerDoParallel(cores=10)
  print(expstring)
  allfolds <- foreach(fold = 1:10

                      ,.options.multicore=mcoptions) %dopar% {

                        prepmodel_pred <- selectModel(exp_pred,model)

                        testf <- makeFreqDataInd(exp_pred,model) %>% .$LOPFolds == fold

                        tmpdat_pred <- prepLOP_Ind(prepmodel_pred$tmpdat,testf,model)

                        # predict out-of-sample probabilities with samples from group-level
                        pred <- predictionLOP_Ind(tmpdat_pred,fit_f,model)

                        # call the fit file for the training data for the specific fold, model, experiment, modeltype
                        fit_f <- readRDS(paste0("../Fits/NoItemEffectsRandomCrit/LOP/",expstring,"/fit_",model,"_LOP_",unique(exp_pred$exp),"_LOPFold",fold,".rds"))


                        # Determine diagnostics
                        divergent <- sum(fit_f$sampler_diagnostics()[,,2])/length(fit_f$sampler_diagnostics()[,,2])

                        fullres <-  data.frame(matrix(pred %>% unlist()*-2,nrow=6000,byrow=T)) %>%
                          mutate(DivChainsPercent = divergent) %>%
                          mutate(model = model,
                                 type = "LOP",
                                 Fold = fold,
                                 experiment = expstring,
                                 manipulation = unique(exp_pred$manipulation),
                                 itemInd = unique(exp_pred$itemInd))

                        return(fullres)

                      }

  saveRDS(allfolds, paste0("../Predictions/NoItemEffectsRandomCrit/LOP/",model,"/DevT1_", model,"_LOP_",unique(exp_pred$exp),".rds"))


}

# Item effects, random crit ----------------------------------------------------

# prediction functions for this kind of model in here
source("LOPTheta.R")

selectdata <- c("APH2016_e1", "D2007_e1b","DW2020_e1", "FGR2019_e1", "FO2016_e2a",
                "KFH2013_e1", "KY2011_e1", "SB2020_e2", "SD2014_e1", "ZOL2021_e1")

for(i in selectdata){

  #i <- "FGR2019_e1"
  model <- "EVgamSDT"

  #fitfiles <- list.files(paste0("../Fits/NoItemEffects/KFCV/",i),full.names = T)

  exp_pred <- alldata[[which(listexp == i)]] %>%
    mutate(rating = as.numeric(as.character(rating)))

  #print(unique(exp_pred$exp))

  mcoptions <- list(preschedule=TRUE, set.seed=FALSE)
  doParallel::registerDoParallel(cores=10)
  print(expstring)
  allfolds <- foreach(fold = 1:10

                      ,.options.multicore=mcoptions) %dopar% {

                        prepmodel_pred <- selectModel(exp_pred,model)

                        tmpdat_pred<- makeLOtest_item(exp_pred,prepmodel_pred$tmpdat,fold,model,"LOP","test")

                        # predict out-of-sample probabilities with samples from group-level
                        pred <- predictionLOP(tmpdat_pred,fit_f,model)

                        # call the fit file for the training data for the specific fold, model, experiment, modeltype
                        fit_f <- readRDS(paste0("../Fits/ItemEffects/LOP/",expstring,"/fit_",model,"_LOP_",unique(exp_pred$exp),"_LOPFold",fold,".rds"))


                        # Determine diagnostics
                        divergent <- sum(fit_f$sampler_diagnostics()[,,2])/length(fit_f$sampler_diagnostics()[,,2])

                        fullres <-  data.frame(matrix(pred %>% unlist()*-2,nrow=6000,byrow=T)) %>%
                          mutate(DivChainsPercent = divergent) %>%
                          mutate(model = model,
                                 type = "LOP",
                                 Fold = fold,
                                 experiment = expstring,
                                 manipulation = unique(exp_pred$manipulation),
                                 itemInd = unique(exp_pred$itemInd))

                        return(fullres)

                      }

  saveRDS(allfolds, paste0("../Predictions/ItemEffects/LOP/",model,"/DevT1_", model,"_LOP_",unique(exp_pred$exp),".rds"))


}