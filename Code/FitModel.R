# General information ------------------------------------------

# Fit full data sets for comparison by Deviance and LOOIC
# Fit training data for K-Fold crossvalidation (KFCV) where K=10
# Fit training data Leave-Participant-Out crossvalidation (LOP)

# fit all for
## models without items effects
## models without item effects, random crit (criteria are given for individuals as deviations from grand mean)
## models with item effects, random crit

# in the fitting routine below, model type is indicated by using or changing column: itemInd in data file
# original hardcoding has itemInd = "Item" for all data sets that have item information and "Ind" otherwise
## itemInd == "Ind" : No item effects
## itemInd == "Indgam" : No item effects, random crit
## itemInd == "Item" : Item effects (only for data sets that actually contain item information)

# Folds for KFCV and LOP were hardcoded into the data files during preprocessing as
## KFCV: Folds column
## LOP: LOPFolds column

# Note: the data prepared for stan (named tmpdat, below) can contain variables not relevant for the model
# Stan only uses the variables actually called in the stanmodel
# in KFCV, LOP only variables relevant to the model are changed

# For more details about models, see the manuscript


#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

#library("rstan")
#library("tidybayes")
#library("brms")
library("cmdstanr") # to run mcmc stan
library("Matrix") # tidyverse seems to need Matrix to work
library("tidyverse")
library("bayestestR") #for model matrices
options(contrasts=c('contr.orthonorm', 'contr.poly'))


rds_file <- list.files("../ProcessedData/", full.names = TRUE)
rds_file <- rds_file[grepl(".rds",rds_file)]

# put all data in a single object
alldata <- purrr::map(c(rds_file), readRDS)

# list data/experiment names
listexp <- alldata %>% purrr::map(. %>% ungroup() %>% distinct(exp)) %>% bind_rows() %>% .$exp
listexpitem <- alldata %>% purrr::map(. %>% filter(itemInd=="Item") %>% ungroup() %>% distinct(exp)) %>% bind_rows() %>% .$exp
#print(listexp)
#print(listexpitem)

# all functions called in this file are in:
source("MakeStanDataModel.R")

# Fit No Item Effects, Full/Deviance fit ---------------------------

for(u in listexp){

  # u <- "FGR2019_e1"

  expdata<- alldata[[which(listexp == u)]] %>%
    mutate(rating = as.integer(rating)) %>% # for some data sets it's as character
    mutate(itemInd = "Ind")

  #models: EVSDT, UVSDT, Gumbel, DPSDT, DPRMSDT, DPRM2SDT, M0SDT, MASDT, 2HTM
  model <- "EVSDT"

  print(unique(expdata$exp))
  print(model)

  tmp<-selectModel(expdata,model)
  # create list with two objects:
  ## tmp$tmpdat = data for stan,
  ## tmp$model = stanmodel file

  tmpdat<- tmp$tmpdat
  #str(tmpdat)

  # make cmdstanr model from file
  stanmodel <- tmp$model %>%
    write_stan_file() %>%
    cmdstan_model() ->
    cmdstanr_model

  tmpdat$exp <- 0 # otherwise crashes

  fit<- stanmodel$sample(
    data = tmpdat,
    adapt_delta = 0.99,
    max_treedepth = 20,
    seed = 667669,
    init = 0.1,
    chains = 6, parallel_chains = 6)

  fit$save_object(file = paste0("../Fits/NoItemEffects/Full/",unique(expdata$exp),"/fit_",model,"_Full_",unique(expdata$exp),".rds"))

}


# Fit No Item Effects, KFCV fit ---------------------------

for(u in listexp){

  expdata<- alldata[[which(listexp == u)]] %>%
    mutate(rating = as.integer(rating)) %>% # for some data sets it's as character
    mutate(itemInd = "Ind")

  #models: EVSDT, UVSDT, Gumbel, DPSDT, DPRMSDT, DPRM2SDT, M0SDT, MASDT, 2HTM
  model <- "EVSDT"

  print(unique(expdata$exp))
  print(model)

  for(fold in c(1:10)){

    #fold <- 1

    tmp<-selectModel(expdata %>% filter(Folds != fold),model) # all but 1 fold

    # create list with two objects:
    ## tmp$tmpdat = data for stan,
    ## tmp$model = stanmodel file

    tmpdat<- tmp$tmpdat
    #str(tmpdat)

    # make cmdstanr model from file
    stanmodel <- tmp$model %>%
      write_stan_file() %>%
      cmdstan_model() ->
      cmdstanr_model

    tmpdat$exp <- 0 # otherwise crashes

    fit<- stanmodel$sample(
      data = tmpdat,
      adapt_delta = 0.99,
      max_treedepth = 20,
      seed = 667669,
      init = 0.1,
      chains = 6, parallel_chains = 6)

    fit$save_object(file = paste0("../Fits/NoItemEffects/KFCV/",model,"/fit_",model,"_KFCV_",unique(expdata$exp),"_Folds",fold,".rds"))
  }
}


# Fit No Item Effects, LOP fit ---------------------------
# LP2013_e2
# RSM2009_e1
# KAWY2013_e2
#


for(model in c("2HTM")){
for(u in listexp[c(60,30)]){

  expdata<- alldata[[which(listexp == u)]] %>%
    mutate(rating = as.integer(rating)) %>% # for some data sets it's as character
    mutate(itemInd = "Ind")

  #models: EVSDT, UVSDT, Gumbel, DPSDT, DPRMSDT, DPRM2SDT, M0SDT, MASDT, 2HTM
  #model <- "MASDT"

  print(unique(expdata$exp))
  print(model)

  for(fold in c(1:10)){

    #fold <- 1

    tmp<-selectModel(expdata %>% filter(LOPFolds != fold),model) # all but 1 fold

    # create list with two objects:
    ## tmp$tmpdat = data for stan,
    ## tmp$model = stanmodel file

    tmpdat<- tmp$tmpdat
    #str(tmpdat)

    # make cmdstanr model from file
    stanmodel <- tmp$model %>%
      write_stan_file() %>%
      cmdstan_model() ->
      cmdstanr_model

    tmpdat$exp <- 0 # otherwise crashes

    fit<- stanmodel$sample(
      data = tmpdat,
      adapt_delta = 0.99,
      max_treedepth = 20,
      seed = 667669,
      init = 0.1,
      chains = 6, parallel_chains = 6)

    fit$save_object(file = paste0("../Fits/NoItemEffects/LOP/",unique(expdata$exp),"/fit_",model,"_LOP_",unique(expdata$exp),"_LOPFold",fold,".rds"))
  }
}
}



# Fit No Item Effects, Random Crit, Full/Deviance fit ---------------------------
## in the paper, we only fitted these models to:
## APH2016_e1, D2007_e1b,DW2020_e1, FGR2019_e1, FO2016_e2a,
## KFH2013_e1, KY2011_e1, SB2020_e2, SD2014_e1, ZOL2021_e1

for(u in listexp){

  expdata<- alldata[[which(listexp == u)]] %>%
    mutate(rating = as.integer(rating)) %>% # for some data sets it's as character
    mutate(itemInd = "Indgam")

  #models: EVgamSDT, UVgamSDT, DPgamSDT, M0gamSDT, MAgamSDT
  model <- "EVgamSDT"

  print(unique(expdata$exp))
  print(model)

  tmp<-selectModel(expdata,model)
  # create list with two objects:
  ## tmp$tmpdat = data for stan,
  ## tmp$model = stanmodel file

  tmpdat<- tmp$tmpdat

  # make cmdstanr model from file
  stanmodel <- tmp$model %>%
    write_stan_file() %>%
    cmdstan_model() ->
    cmdstanr_model

  tmpdat$exp <- 0 # otherwise crashes


  # Gumbel and DPRMSDT/DPRM2SDT need init values for some parameters
  # for other models defaults to 0.1
  # six chains
  inits <- list(initfun(tmpdat,model),
                initfun(tmpdat,model),
                initfun(tmpdat,model),
                initfun(tmpdat,model),
                initfun(tmpdat,model),
                initfun(tmpdat,model))



  fit<- stanmodel$sample(
    data = tmpdat,
    adapt_delta = 0.99,
    max_treedepth = 20,
    seed = 667669,
    init = inits,#0.1,
    chains = 6, parallel_chains = 6)

  fit$save_object(file = paste0("../Fits/NoItemEffectsRandomCrit/Full/",unique(expdata$exp),"/fit_",model,"_Full_",unique(expdata$exp),".rds"))

}


# Fit No Item Effects, Random Crit, KFCV fit ---------------------------

for(u in listexp){

  expdata<- alldata[[which(listexp == u)]] %>%
    mutate(rating = as.integer(rating)) %>% # for some data sets it's as character
    mutate(itemInd = "Indgam")

  #models: EVgamSDT, UVgamSDT, DPgamSDT, M0gamSDT, MAgamSDT
  model <- "EVgamSDT"

  print(unique(expdata$exp))
  print(model)

  for(fold in c(1:10)){

    tmp<-selectModel(expdata %>% filter(Folds != fold),model) # all but 1 fold

    # create list with two objects:
    ## tmp$tmpdat = data for stan,
    ## tmp$model = stanmodel file

    tmpdat<- tmp$tmpdat

    # make cmdstanr model from file
    stanmodel <- tmp$model %>%
      write_stan_file() %>%
      cmdstan_model() ->
      cmdstanr_model

    # Gumbel and DPRMSDT/DPRM2SDT need init values for some parameters
    # for other models defaults to 0.1
    # six chains
    inits <- list(initfun(tmpdat,model),
                  initfun(tmpdat,model),
                  initfun(tmpdat,model),
                  initfun(tmpdat,model),
                  initfun(tmpdat,model),
                  initfun(tmpdat,model))


    tmpdat$exp <- 0 # otherwise crashes

    fit<- stanmodel$sample(
      data = tmpdat,
      adapt_delta = 0.99,
      max_treedepth = 20,
      seed = 667669,
      init = inits,#0.1,
      chains = 6, parallel_chains = 6)

    fit$save_object(file = paste0("../Fits/NoItemEffectsRandomCrit/KFCV/",model,"/fit_",model,"_KFCV_",unique(expdata$exp),"_Folds",fold,".rds"))
  }
}


# Fit No Item Effects, Random Crit, LOP fit ---------------------------

for(u in listexp){

  expdata<- alldata[[which(listexp == u)]] %>%
    mutate(rating = as.integer(rating)) %>% # for some data sets it's as character
    mutate(itemInd = "Indgam")

  #models: EVgamSDT, UVgamSDT, DPgamSDT, M0gamSDT, MAgamSDT
  model <- "EVgamSDT"

  print(unique(expdata$exp))
  print(model)

  for(fold in c(1:10)){

    tmp<-selectModel(expdata %>% filter(LOPFolds != fold),model) # all but 1 fold

    # create list with two objects:
    ## tmp$tmpdat = data for stan,
    ## tmp$model = stanmodel file

    tmpdat<- tmp$tmpdat

    # make cmdstanr model from file
    stanmodel <- tmp$model %>%
      write_stan_file() %>%
      cmdstan_model() ->
      cmdstanr_model

    tmpdat$exp <- 0 # otherwise crashes

    # Gumbel and DPRMSDT/DPRM2SDT need init values for some parameters
    # for other models defaults to 0.1
    # six chains
    inits <- list(initfun(tmpdat,model),
                  initfun(tmpdat,model),
                  initfun(tmpdat,model),
                  initfun(tmpdat,model),
                  initfun(tmpdat,model),
                  initfun(tmpdat,model))


    fit<- stanmodel$sample(
      data = tmpdat,
      adapt_delta = 0.99,
      max_treedepth = 20,
      seed = 667669,
      init = inits,#0.1,
      chains = 6, parallel_chains = 6)

    fit$save_object(file = paste0("../Fits/NoItemEffectsRandomCrit/LOP/",model,"/fit_",model,"_LOP_",unique(expdata$exp),"_LOPFold",fold,".rds"))
  }
}




# Fit Item Effects, Random Crit, Full/Deviance fit ---------------------------
## this is only possible for data sets with item information, designated by
## itemInd == "Item" (listexpitem)
## in the paper we only fitted these models to:
## APH2016_e1, D2007_e1b,DW2020_e1, FGR2019_e1, FO2016_e2a,
## KFH2013_e1, KY2011_e1, SB2020_e2, SD2014_e1, ZOL2021_e1

for(u in listexpitem){

  #u <- "FGR2019_e1"
  expdata<- alldata[[which(listexp == u)]] %>%
    mutate(rating = as.integer(rating)) # for some data sets it's as character


  #models: EVgamSDT, UVgamSDT, DPgamSDT, M0gamSDT, MAgamSDT
  model <- "EVgamSDT"

  print(unique(expdata$exp))
  print(model)

  tmp<-selectModel(expdata,model)
  # create list with two objects:
  ## tmp$tmpdat = data for stan,
  ## tmp$model = stanmodel file

  tmpdat<- tmp$tmpdat
  #str(tmpdat)

  # make cmdstanr model from file
  stanmodel <- tmp$model %>%
    write_stan_file() %>%
    cmdstan_model() ->
    cmdstanr_model

  tmpdat$exp <- 0 # otherwise crashes

  fit<- stanmodel$sample(
    data = tmpdat,
    adapt_delta = 0.99,
    max_treedepth = 20,
    seed = 667669,
    init = 0.1,
    chains = 6, parallel_chains = 6)

  # Note the accidentally identical file names. Fits with item effects will include
  # parameters referring to items (sd_2 etc)

  fit$save_object(file = paste0("../Fits/ItemEffects/Full/",unique(expdata$exp),"/fit_",model,"_Full_",unique(expdata$exp),".rds"))

}


# Fit Item Effects, Random Crit, KFCV fit ---------------------------


for(u in listexpitem){

  expdata<- alldata[[which(listexp == u)]] %>%
    mutate(rating = as.integer(rating)) # for some data sets it's as character


  #models: EVgamSDT, UVgamSDT, DPgamSDT, M0gamSDT, MAgamSDT
  model <- "EVgamSDT"

  print(unique(expdata$exp))
  print(model)

  for(fold in c(1:10)){

    tmp<-selectModel(expdata,model) # all but 1 fold

    # create list with two objects:
    ## tmp$tmpdat = data for stan,
    ## tmp$model = stanmodel file

    # leave out ~ 10% of trials
    tmpdat<- makeLO_item(exp = expdata,
                         tmpdat = tmp$tmpdat,
                         fold = fold,
                         model = model,
                         type = "KFCV", #leave out ~10% of trials
                         state = "train")
    #str(tmpdat)


    # make cmdstanr model from file
    stanmodel <- tmp$model %>%
      write_stan_file() %>%
      cmdstan_model() ->
      cmdstanr_model

    tmpdat$exp <- 0 # otherwise crashes

    fit<- stanmodel$sample(
      data = tmpdat,
      adapt_delta = 0.99,
      max_treedepth = 20,
      seed = 667669,
      init = 0.1,
      chains = 6, parallel_chains = 6)

    fit$save_object(file = paste0("../Fits/ItemEffects/KFCV/",model,"/fit_",model,"_KFCV_",unique(expdata$exp),"_Folds",fold,".rds"))
  }
}


# Fit Item Effects, Random Crit, LOP fit ---------------------------

for(u in listexpitem){

  expdata<- alldata[[which(listexp == u)]] %>%
    mutate(rating = as.integer(rating)) # for some data sets it's as character


  #models: EVgamSDT, UVgamSDT, DPgamSDT, M0gamSDT, MAgamSDT
  model <- "EVgamSDT"

  print(unique(expdata$exp))
  print(model)

  for(fold in c(1:10)){

    tmp<-selectModel(expdata %>% filter(LOPFolds != fold),model) # all but 1 fold

    # create list with two objects:
    ## tmp$tmpdat = data for stan,
    ## tmp$model = stanmodel file


    # leave out ~ 10% of trials
    tmpdat<- makeLO_item(exp = expdata,
                         tmpdat = tmp$tmpdat,
                         fold = fold,
                         model = model,
                         type = "LOP", #leave out ~10% of participants
                         state = "train")

    # make cmdstanr model from file
    stanmodel <- tmp$model %>%
      write_stan_file() %>%
      cmdstan_model() ->
      cmdstanr_model

    tmpdat$exp <- 0 # otherwise crashes

    fit<- stanmodel$sample(
      data = tmpdat,
      adapt_delta = 0.99,
      max_treedepth = 20,
      seed = 667669,
      init = 0.1,
      chains = 6, parallel_chains = 6)

    fit$save_object(file = paste0("../Fits/ItemEffects/LOP/",model,"/fit_",model,"_LOP_",unique(expdata$exp),"_LOPFold",fold,".rds"))
  }
}

