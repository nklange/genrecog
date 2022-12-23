# General information ------------------------------------------

# Use MCMC samples (fit) to generate
## prediction (Deviance)
## posterior p-values
## LOOIC (we generated LOOIC by running predictions a second time
# -> separate files in Predictions/MODELTYPE/Full/LOOIC folder
# cleaned the code here to make it a single step output)

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

#library("rstan")
#library("tidybayes")
#library("brms")
library("cmdstanr")
library("Matrix")
library("tidyverse")
#library("bayestestR")
#library("foreach")
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

  fitfiles <- list.files(paste0("../Fits/NoItemEffects/Full/",model),full.names = T)

  exp_pred <- alldata[[which(listexp == i)]] %>%
    mutate(rating = as.numeric(as.character(rating))) %>%
    mutate(itemInd = "Ind")

  #print(unique(exp_pred$exp))

  prepmodel_pred <- selectModel(exp_pred,model)
  tmpdat_pred <- prepmodel_pred$tmpdat

  fit_f <- readRDS(fitfiles[grepl(unique(exp_pred$exp),fitfiles)])

  # predicted probabilities, when having fit full data set / for training data
  pred <- predictionFULL_Ind(tmpdat_pred,fit_f,model)

  # make Deviance value
  Dev <- getDeviancefromTheta_Ind(tmpdat_pred$Y,pred)
  # ALTERNATIVE, for deviance summed by participant rather than summed across participants:
  #  ids <- makeFreqDataInd(exp_pred,model)
  #  Dev <- getDeviancefromTheta_Ind_individual(tmpdat_pred$Y,pred,ids$id)

  # Determine posterior p-values
  T1s <- parallel::mclapply(pred,function(x) T1variouslevel(x,exp_pred,model), mc.cores = 18, mc.preschedule = TRUE)

  # Determine diagnostics
  divergent <- sum(fit_f$sampler_diagnostics()[,,2])/length(fit_f$sampler_diagnostics()[,,2])

  # make LOOIC
  ids <- makeFreqDataInd(exp_pred,model)
  pwLL <- getpwLL_Ind(ids,pred) # point-wise log-probs
  matrixpwLL<-t(array(unlist(pwLL), dim = c(dim(exp_pred)[1],6000)))

  loop <- loo::loo(matrixpwLL, reloo = T,r_eff = NA,cores = 1)



  information <-
    tibble(model = modelpred,
           type = "Full",
           experiment = unique(exp_pred$exp),
           manipulation = unique(exp_pred$manipulation),
           itemInd = unique(exp_pred$itemInd))
  deviances <-  tibble(Dev = Dev) %>%
    mutate(DivChainsPercent = divergent)

  saveRDS(list(information = information,
               deviances = deviances,
               LOOIC = loop,
               T1s = lapply(T1s, function(x) x$T1) ,
               predictedfrequency = lapply(T1s, function(x) x$avpredfreq)),
          paste0("../Predictions/NoItemEffects/Full/",model,"/DevT1_", model,"_full_",unique(exp_pred$exp),".rds"))

  # when saving by-individual deviance
  # paste0("../Predictions/NoItemEffects/Full/IndividualFull/DevT1_", model,"_full_",unique(exp_pred$exp),"_individual.rds"))


}

# No Item Effects, Random Crit -------------------------------------------------

# prediction functions for this kind of model in here
source("LOPThetaInd_gam.R")
fitfiles <- list.files(paste0("../Fits/NoItemEffectsRandomCrit/Full/"),full.names = T)

selectdata <- c("APH2016_e1", "D2007_e1b","DW2020_e1", "FGR2019_e1", "FO2016_e2a",
                "KFH2013_e1", "KY2011_e1", "SB2020_e2", "SD2014_e1", "ZOL2021_e1")


for(i in selectdata){


  model <- "EVgamSDT"

  exp_pred <- alldata[[which(listexp == i)]] %>%
    mutate(rating = as.numeric(as.character(rating))) %>%
    mutate(itemInd = "Indgam")

  #print(unique(exp_pred$exp))

  prepmodel_pred <- selectModel(exp_pred,model)
  tmpdat_pred <- prepmodel_pred$tmpdat

  fittmp <- fitfiles[grepl(unique(exp_pred$exp),fitfiles)]
  fit_f <- readRDS(fittmp[grepl(model,fittmp)])

  # predicted probabilities, when having fit full data set / for training data
  pred <- predictionFULL_Ind(tmpdat_pred,fit_f,model)

  # make Deviance value
  Dev <- getDeviancefromTheta_Ind(tmpdat_pred$Y,pred)
  # ALTERNATIVE, for deviance summed by participant rather than summed across participants:
  #  ids <- makeFreqDataInd(exp_pred,model)
  #  Dev <- getDeviancefromTheta_Ind_individual(tmpdat_pred$Y,pred,ids$id)

  # Determine posterior p-values
  T1s <- parallel::mclapply(pred,function(x) T1variouslevel(x,exp_pred,model), mc.cores = 18, mc.preschedule = TRUE)

  # Determine diagnostics
  divergent <- sum(fit_f$sampler_diagnostics()[,,2])/length(fit_f$sampler_diagnostics()[,,2])

  # make LOOIC
  ids <- makeFreqDataInd(exp_pred,model)
  pwLL <- getpwLL_Ind(ids,pred) # point-wise log-probs
  matrixpwLL<-t(array(unlist(pwLL), dim = c(dim(exp_pred)[1],6000)))

  loop <- loo::loo(matrixpwLL, reloo = T,r_eff = NA,cores = 1)

  information <-
    tibble(model = modelpred,
           type = "Full",
           experiment = unique(exp_pred$exp),
           manipulation = unique(exp_pred$manipulation),
           itemInd = unique(exp_pred$itemInd))
  deviances <-  tibble(Dev = Dev) %>%
    mutate(DivChainsPercent = divergent)

  saveRDS(list(information = information,
               deviances = deviances,
               LOOIC = loop,
               T1s = lapply(T1s, function(x) x$T1) ,
               predictedfrequency = lapply(T1s, function(x) x$avpredfreq)),
          paste0("../Predictions/NoItemEffectsRandomCrit/Full/",model,"/DevT1_", model,"_full_",unique(exp_pred$exp),".rds"))

  # when saving by-individual deviance
  # paste0("../Predictions/NoItemEffectsRandomCrit/Full/IndividualFull/DevT1_", model,"_full_",unique(exp_pred$exp),"_individual.rds"))


}

# Item effects, random crit ----------------------------------------------------

# prediction functions for this kind of model in here
source("LOPTheta.R")
fitfiles <- list.files(paste0("../Fits/ItemEffectsRandomCrit/Full/"),full.names = T)

selectdata <- c("APH2016_e1", "D2007_e1b","DW2020_e1", "FGR2019_e1", "FO2016_e2a",
                "KFH2013_e1", "KY2011_e1", "SB2020_e2", "SD2014_e1", "ZOL2021_e1")

for(i in selectdata){

  model <- "EVgamSDT"

  exp_pred <- alldata[[which(listexp == i)]] %>%
    mutate(rating = as.numeric(as.character(rating)))

  #print(unique(exp_pred$exp))

  prepmodel_pred <- selectModel(exp_pred,model)
  tmpdat_pred <- prepmodel_pred$tmpdat

  fittmp <- fitfiles[grepl(unique(exp_pred$exp),fitfiles)]
  fit_f <- readRDS(fittmp[grepl(model,fittmp)])

  # predicted probabilities, when having fit full data set / for training data
  pred <- predictionFULL(tmpdat_pred,fit_f,model)

  # make Deviance value
  Dev <- getDeviancefromTheta(tmpdat_pred$Y,pred)
  # ALTERNATIVE, for deviance summed by participant rather than summed across participants:
  #  Dev <- getDeviancefromTheta_individual(tmpdat_pred$Y,pred,exp_pred$id)

  # Determine posterior p-values
  T1s <- parallel::mclapply(pred,function(x) T1variouslevel(x,exp_pred,model), mc.cores = 18, mc.preschedule = TRUE)

  # Determine diagnostics
  divergent <- sum(fit_f$sampler_diagnostics()[,,2])/length(fit_f$sampler_diagnostics()[,,2])

  # make LOOIC

  pwLL <- getpwLL(tmpdat_pred$Y,pred)
  matrixpwLL<-t(array(unlist(pwLL), dim = c(length(tmpdat_pred$Y),6000)))
  loop <- loo::loo(matrixpwLL, reloo = T,r_eff = NA,cores = 1)

  information <-
    tibble(model = modelpred,
           type = "Full",
           experiment = unique(exp_pred$exp),
           manipulation = unique(exp_pred$manipulation),
           itemInd = unique(exp_pred$itemInd))
  deviances <-  tibble(Dev = Dev) %>%
    mutate(DivChainsPercent = divergent)

  saveRDS(list(information = information,
               deviances = deviances,
               LOOIC = loop,
               T1s = lapply(T1s, function(x) x$T1) ,
               predictedfrequency = lapply(T1s, function(x) x$avpredfreq)),
          paste0("../Predictions/ItemEffectsRandomCrit/Full/",model,"/DevT1_", model,"_full_",unique(exp_pred$exp),".rds"))

  # when saving by-individual deviance
  # paste0("../Predictions/ItemEffectsRandomCrit/Full/IndividualFull/DevT1_", model,"_full_",unique(exp_pred$exp),"_individual.rds"))


}