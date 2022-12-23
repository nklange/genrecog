# header just in case at the top of all files even where they are sourced.

#library("rstan")
#library("tidybayes")
#library("brms")
library("cmdstanr")
library("Matrix")
library("tidyverse")
#library("bayestestR")

options(contrasts=c('contr.orthonorm', 'contr.poly'))

# path for stanfiles:
stanfilepath <- "Stanfiles/"

# select model by type of implementation ---------------

selectModel <- function(exp,modelname){


  if(unique(exp$itemInd == "Item")){

    information <- ItemModel(exp,modelname)

  } else if (unique(exp$itemInd == "Ind")){

    information <- IndModel(exp,modelname)

  } else if (unique(exp$itemInd == "Indgam")){

    information <- IndgamModel(exp,modelname)
  }

  information
}

# combine standata and stan model file for given data set -------------

IndModel <- function(exp,modelname){

  # makes standata
  expdats <- makeFreqDataInd(exp,modelname)
  tmpdat <- makeTmpDatInd(expdats,modelname)

  # selects stan file
  models <- c("IndStanTest_UVSDT.stan",
              "IndStanTest_UVSDT_None.stan",
              "IndStanTest_EVSDT.stan",
              "IndStanTest_EVSDT_None.stan",
              "IndStanTest_Gumbel.stan",
              "IndStanTest_Gumbel_None.stan",
              "IndStanTest_DPSDT_None.stan",
              "IndStanTest_DPSDT.stan",
              "IndStanTest_DPRMSDT_None.stan",
              "IndStanTest_DPRMSDT.stan",
              "IndStanTest_DPRM2SDT_None.stan",
              "IndStanTest_DPRM2SDT.stan",
              "IndStanTest_M0SDT.stan",
              "IndStanTest_M0SDT_None.stan",
              "IndStanTest_MASDT.stan",
              "IndStanTest_MASDT_None.stan",
              "IndStanTest_2HTM.stan", #multiple sets of old and new items
              "IndStanTest_2HTM_ws.stan", #1 set of new items, multiple sets of old items
              "IndStandTest_2HTM_ws2.stan", #1 set of old items, multiple sets of new items
              "IndStanTest_2HTM_None.stan")
  stanmodel <- models[grepl(modelname,models)]



  if(unique(exp$manipulation) == "None"){

    stanmodel_use <- stanmodel[grepl("None",stanmodel)]

    if(modelname=="2HTM"){

      if(unique(exp$exp) %in% c("LBA2019_e1","LBA2019_e2")){

        stanmodel_use <- "IndStanTest_2HTM_None_neutral.stan"
        # uneven number of confidence ratings, i.e. one middle rating

      } else {
        stanmodel_use <- "IndStanTest_2HTM_None.stan"
      }
    }

  } else {

    if(modelname=="2HTM") {

      if(unique(exp$exp) %in% c("KY2016_e1","GKH1999_e2", "GKH1999_e3","GKH1999_e4",
                                "CSR2015_e1","FBH2013_e1","FO2016_e2a","KAWY2013_e3",
                                "MKG2018_e2","RSP2012_e1","USS2015_e1","WMS2012_e1","ZMD2011_e1") ){

        stanmodel_use <- "IndStanTest_2HTM_ws.stan"

      } else if(unique(exp$exp) %in% c("HUW2015_e1","LM2020_e1")){

        stanmodel_use <- "IndStanTest_2HTM_neutral.stan"

      } else if (unique(exp$exp) %in% c("LP2013_e2")) {

        stanmodel_use <-"IndStanTest_2HTM_ws2.stan"
      }  else {

        stanmodel_use <- "IndStanTest_2HTM.stan"
      }

    } else {
      stanmodel_use <- stanmodel[!grepl("None",stanmodel)]
    }
  }


  list(model = readLines(paste0(stanfilepath,"NoItemEffects/",stanmodel_use)),
       tmpdat = tmpdat)


}
IndgamModel <- function(exp,modelname){

  # makes standata
  expdats <- makeFreqDataInd(exp,modelname)
  tmpdat <- makeTmpDatInd(expdats,modelname)

  # selects stanfile
  # 2HTM has no criteria so 2HTM-Ind == 2HTM-Indgam
  models <- c("IndStanTest_UVgamSDT_None.stan",
              "IndStanTest_UVgamSDT.stan",
              "IndStanTest_EVgamSDT_None.stan",
              "IndStanTest_EVgamSDT.stan",
              "IndStanTest_M0gamSDT_None.stan",
              "IndStanTest_M0gamSDT.stan",
              "IndStanTest_DPgamSDT_None.stan",
              "IndStanTest_DPgamSDT.stan",
              "IndStanTest_MAgamSDT_None.stan",
              "IndStanTest_MAgamSDT.stan",
              "IndStanTest_DPRMgamSDT_None.stan",
              "IndStanTest_DPRMgamSDT.stan",
              "IndStanTest_DPRM2gamSDT_None.stan",
              "IndStanTest_DPRM2gamSDT.stan")

  stanmodel <- models[grepl(modelname,models)]


  if(unique(exp$manipulation) == "None"){

    stanmodel_use <- stanmodel[grepl("None",stanmodel)]

  } else {

    stanmodel_use <- stanmodel[!grepl("None",stanmodel)]

  }


  list(model =  readLines(paste0(stanfilepath,"NoItemEffectsRandomCrit/",stanmodel_use)),
       tmpdat = tmpdat)


}
ItemModel <- function(exp,modelname){

  # in SDT models, 1 set of new items by designation where multiple sets
  # in 2HTM separate parameters for old and new items -> use condition_true
  # instead of the rejigged SDT conditions

  # makes standata
  if(modelname == "2HTM" &
     unique(exp$exp) %in% c("BSG2014_e1","D2007_e1a","D2007_e1b",
                            "KFH2013_e1",
                            "KK2015_e1","RS2009_e1","RS2009_e2",
                            "SBT2018_e1","TMP2014_e1","ZOL2021_e1"

     )){

    if(unique(exp$exp) == "BSG2014_e1"){

      exp <- exp %>% mutate(condition_true = paste0(condition_time,"_",condition_val))
    } else if (unique(exp$exp) == "SBT2018_e1"){

      exp <- exp %>% mutate(condition_true = paste0(condition_time,"_",condition_valence))
    } else if (unique(exp$exp) == "RS2009_e2"){

      exp <- exp %>% mutate(condition_true = paste0(condition_freq,"_",condition_inst))
    } else if (unique(exp$exp) == "KFH2013_e1"){

      exp <- exp %>%
        mutate(condition_true = paste0(condition_acc,"_",condition_freq))
    } else if (unique(exp$exp) == "GKH1999_e1"){

      exp <- exp %>%
        mutate(condition_true = paste0(condition_freq,"_",condition_studyduration))
    }

    exp <- exp %>%
      mutate(isold=isold_true,
             condition = condition_true)
  }


  tmpdat <- makeTmpDatItem(exp,modelname)

  # selects stanfile
  models <- c("ItemStanTest_UVgamSDT.stan",
              "ItemStanTest_UVgamSDT_None.stan",
              "ItemStanTest_EVgamSDT.stan",
              "ItemStanTest_EVgamSDT_None.stan",
              "ItemStanTest_DPgamSDT_None.stan",
              "ItemStanTest_DPgamSDT.stan",
              "ItemStanTest_M0SDT_None.stan",
              "ItemStanTest_M0SDT.stan",
              "ItemStanTest_MASDT.stan",
              "ItemStanTest_MASDT_None.stan",
              "ItemStanTest_2HTM.stan",
              # other 2HTM variations not needed for the select data sets we
              # fitted with item effects
              "ItemStanTest_2HTM_ws.stan",
              "ItemStanTest_2HTM_None.stan")
  stanmodel <- models[grepl(modelname,models)]


  if(unique(exp$manipulation) == "None"){

    stanmodel_use <- stanmodel[grepl("None",stanmodel)]

    if(modelname=="2HTM"){

      if(unique(exp$exp) %in% c("LBA2019_e1","LBA2019_e2")){

        stanmodel_use <- "ItemStanTest_2HTM_None_neutral.stan" # does not exist

      } else {
        stanmodel_use <- "ItemStanTest_2HTM_None.stan"
      }
    }

  } else {

    if(modelname=="2HTM") {

      if(unique(exp$exp) %in% c("KY2016_e1","GKH1999_e2", "GKH1999_e3","GKH1999_e4",
                                "CSR2015_e1","FBH2013_e1","FO2016_e2a","KAWY2013_e3",
                                "MKG2018_e2","RSP2012_e1","USS2015_e1","WMS2012_e1","ZMD2011_e1") ){

        stanmodel_use <- "ItemStanTest_2HTM_ws.stan" # does not exist

      } else if(unique(exp$exp) %in% c("HUW2015_e1","LM2020_e1")){

        stanmodel_use <- "ItemStanTest_2HTM_neutral.stan" #does not exist

      } else {

        stanmodel_use <- "ItemStanTest_2HTM.stan"
      }

    } else {
      stanmodel_use <- stanmodel[!grepl("None",stanmodel)]
    }
  }


  list(model =  readLines(paste0(stanfilepath,"ItemEffects/",stanmodel_use)),
       tmpdat = tmpdat)


}

# Prep stan data ------------------------------------------------------
# like stan file based on brms code, implemented manually because more flexible

# prep data for Ind and Indgam models
# turn by-trial data into multinomial data
makeFreqDataInd <- function(exp,model){

  if(model == "2HTM" &
     unique(exp$exp) %in% c("BSG2014_e1","D2007_e1a","D2007_e1b",
                            "KAWY2013_e2","KFH2013_e1",
                            "KK2015_e1","RS2009_e1","RS2009_e2",
                            "SBT2018_e1","TMP2014_e1","ZOL2021_e1"

     )){

    if(unique(exp$exp) == "BSG2014_e1"){

      exp <- exp %>% mutate(condition_true = paste0(condition_time,"_",condition_val))
    } else if (unique(exp$exp) == "SBT2018_e1"){

      exp <- exp %>% mutate(condition_true = paste0(condition_time,"_",condition_valence))
    } else if (unique(exp$exp) == "RS2009_e2"){

      exp <- exp %>% mutate(condition_true = paste0(condition_freq,"_",condition_inst))
    } else if (unique(exp$exp) == "KFH2013_e1"){

      exp <- exp %>%
        mutate(condition_true = paste0(condition_acc,"_",condition_freq))
    } else if (unique(exp$exp) == "GKH1999_e1"){

      exp <- exp %>%
        mutate(condition_true = paste0(condition_freq,"_",condition_studyduration))
    }

    expdats <- exp %>%
      mutate(isold=isold_true,
             condition = condition_true) %>%
      group_by(exp, manipulation,condition, LOPFolds, id,isold,rating) %>%
      summarize(cums = length(rating)) %>%
      mutate(rating = factor(rating, levels = c(1:max(rating))))%>%
      arrange(rating) %>%
      pivot_wider(id_cols = c("exp","manipulation","condition","LOPFolds","id","isold"),
                  names_from = "rating",
                  values_from = "cums")

    expdats[is.na(expdats)] <- 0


    expdats <- expdats %>% arrange(exp,id,isold,condition)

  } else if (model == "2HTM" & unique(exp$exp) == "GKH1999_e1"){


    expdats <- exp %>%
      mutate(isold=isold_true) %>%
      mutate(condition_studyduration = ifelse(isold == "0","Long",condition_studyduration)) %>%
      group_by(exp, manipulation,condition_freq, condition_studyduration, LOPFolds, id,isold,rating) %>%
      summarize(cums = length(rating)) %>%
      mutate(rating = factor(rating, levels = c(1:max(rating))))%>%
      arrange(rating) %>%
      pivot_wider(id_cols = c("exp","manipulation","condition_freq","condition_studyduration","LOPFolds","id","isold"),
                  names_from = "rating",
                  values_from = "cums")
    expdats[is.na(expdats)] <- 0


    expdats <- expdats %>% arrange(exp,id,isold,condition_freq,condition_studyduration)


  } else if (model == "2HTM" & unique(exp$exp) == "LP2013_e2"){


    expdats <- exp %>%
      mutate(isold=isold_true) %>%
      # mutate(condition = condition_true) %>%
      mutate(condition = ifelse(condition_true == "target","related_lure",condition_true)) %>%

      group_by(exp, manipulation,condition, LOPFolds, id,isold,rating) %>%
      summarize(cums = length(rating)) %>%
      mutate(rating = factor(rating, levels = c(1:max(rating))))%>%
      arrange(rating) %>%
      pivot_wider(id_cols = c("exp","manipulation","condition","LOPFolds","id","isold"),
                  names_from = "rating",
                  values_from = "cums")
    expdats[is.na(expdats)] <- 0


    expdats <- expdats %>% arrange(exp,id,isold,condition)

    } else {


    expdats <- exp %>%
      group_by(exp, manipulation,condition, LOPFolds, id,isold,rating) %>%
      summarize(cums = length(rating)) %>%
      mutate(rating = factor(rating, levels = c(1:max(rating))))%>%
      arrange(rating) %>%
      pivot_wider(id_cols = c("exp","manipulation","condition","LOPFolds","id","isold"),
                  names_from = "rating",
                  values_from = "cums")

    # expdats <- exp %>%
    #   group_by(exp, manipulation,condition_freq, condition_studyduration, LOPFolds, id,isold,rating) %>%
    #   summarize(cums = length(rating)) %>%
    #   mutate(rating = factor(rating, levels = c(1:max(rating))))%>%
    #   arrange(rating) %>%
    #   pivot_wider(id_cols = c("exp","manipulation","condition_freq","condition_studyduration","LOPFolds","id","isold"),
    #               names_from = "rating",
    #               values_from = "cums")

    expdats[is.na(expdats)] <- 0


    expdats <- expdats %>% arrange(exp,id,isold,condition)
  }

  expdats
}

# make standata, largely make design matrices.
makeTmpDatInd <- function(expdats,model){

  tmpdat <- NULL

  tmpdat$exp <- unique(expdats$exp)

  if(unique(expdats$manipulation) != "None"){

    if(model == "2HTM"){

      # 2HTM needs a lot of special cases to deal with multiple sets of new items

      if(unique(expdats$exp) == "GKH1999_e1"){


        newX <-  model.matrix(~0+factor(isold) + factor(isold):condition_freq + factor(isold):condition_freq:condition_studyduration,expdats)

        selectorth <- as_tibble(newX)

      } else {

        newX <-  model.matrix(~0+factor(isold) + factor(isold):condition,expdats)

        selectorth <- as_tibble(newX)
      }

    } else {
      newX <- model.matrix(~isold * condition,expdats)

      selectorth <- as_tibble(newX) %>% select(c("isold",matches("isold:")))
    }

  } else if(unique(expdats$manipulation) == "None"){

    if(model == "2HTM"){

      newX <-model.matrix(~ 0 + isold,expdats %>% mutate(isold = factor(isold)))

      selectorth <- newX


    } else {

      newX <- model.matrix(~isold,expdats)

      selectorth <- as_tibble(newX) %>% select(c("isold",matches("isold:")))
    }

  }


  tmpdat$N <- dim(expdats)[1]

  tmpdat$Y <- as.matrix(expdats[,7:dim(expdats)[2]])


  tmpdat$C <- dim(expdats)[2] - 6
  tmpdat$Chalf <- ifelse(model == "DPRM2SDT",2,floor(tmpdat$C / 2))
  tmpdat$nthres <- tmpdat$C - 1



  if(model == "2HTM"){



    if(unique(expdats$manipulation == "None")){

      tmpdat$K_DO <-1
      tmpdat$K_DN <-1
      tmpdat$K_go <-1
      tmpdat$K_neut <-1

      tmpdat$G_DO <- 1
      tmpdat$G_DN <- 1
      tmpdat$G_go <- 1
      tmpdat$G_neut <-1

      tmpdat$X_DO <- as.matrix(as_tibble(selectorth) %>% select("isold1") )
      tmpdat$X_DN <- as.matrix(as_tibble(selectorth) %>% select("isold0") )
      tmpdat$X_go <- as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("isold0") ))[[1]],ncol=1))
      tmpdat$X_neut <- as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("isold0") ))[[1]],ncol=1))
      tmpdat$Z_DO <- as.matrix(as_tibble(selectorth) %>% select("isold1") )
      tmpdat$Z_DN <- as.matrix(as_tibble(selectorth) %>% select("isold0") )
      tmpdat$Z_go <- as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("isold0") ))[[1]],ncol=1))
      tmpdat$Z_neut <- as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("isold0") ))[[1]],ncol=1))



    } else if (unique(expdats$manipulation == "Between")){


      tmpdat$X_DO <- as.matrix(selectorth %>% select(matches("isold)1")))
      tmpdat$X_DN <-  as.matrix(selectorth %>% select(matches("isold)0")))
      tmpdat$X_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
      tmpdat$X_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


      tmpdat$Z_DO <-  as.matrix(as_tibble(selectorth) %>% select("factor(isold)1") )
      tmpdat$Z_DN <- as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") )
      tmpdat$Z_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
      tmpdat$Z_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


      tmpdat$G_DO <- 1
      tmpdat$G_DN <- 1
      tmpdat$G_go <- 1
      tmpdat$G_neut <- 1
      tmpdat$K_DO <- length(unique(expdats$condition))
      tmpdat$K_DN <- length(unique(expdats$condition))
      tmpdat$K_go <- 1
      tmpdat$K_neut <- 1

    } else if (unique(expdats$manipulation == "Within")){

      if(unique(expdats$exp) %in% c("KY2016_e1","GKH1999_e2",
                                    "GKH1999_e3",
                                    "GKH1999_e4",
                                    "CSR2015_e1",
                                    "FBH2013_e1",
                                    "FO2016_e2a",
                                    "KAWY2013_e3",

                                    "MKG2018_e2",
                                    "RSP2012_e1",
                                    "USS2015_e1",
                                    "WMS2012_e1",
                                    "ZMD2011_e1"

      )){
        #head(expdats)
        tmpdat$X_DO <- as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$X_DN <-  as.matrix(selectorth %>% select("factor(isold)0"))
        tmpdat$X_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$X_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


        tmpdat$Z_DO <-  as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$Z_DN <- as.matrix(selectorth %>% select("factor(isold)0"))
        tmpdat$Z_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$Z_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


        tmpdat$G_DO <- length(unique(expdats$condition))
        tmpdat$G_DN <- 1
        tmpdat$G_go <- 1
        tmpdat$G_neut <- 1
        tmpdat$K_DO <- length(unique(expdats$condition))
        tmpdat$K_DN <- 1
        tmpdat$K_go <- 1
        tmpdat$K_neut <- 1

      } else if(unique(expdats$exp) == "GKH1999_e1"){

        tmpdat$X_DO <- as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$X_DN <-  as.matrix(selectorth %>% select("factor(isold)0","factor(isold)0:condition_freq1"))
        tmpdat$X_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$X_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


        tmpdat$Z_DO <-  as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$Z_DN <- as.matrix(selectorth %>% select("factor(isold)0","factor(isold)0:condition_freq1"))
        tmpdat$Z_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$Z_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


        tmpdat$G_DO <- dim(tmpdat$Z_DO)[2]
        tmpdat$G_DN <- dim(tmpdat$Z_DN)[2]
        tmpdat$G_go <- 1
        tmpdat$G_neut <- 1
        tmpdat$K_DO <- dim(tmpdat$X_DO)[2]
        tmpdat$K_DN <- dim(tmpdat$X_DN)[2]
        tmpdat$K_go <- 1
        tmpdat$K_neut <- 1

        tmpdat$C <- dim(expdats)[2] - 7
        tmpdat$Chalf <- floor(tmpdat$C / 2)
        tmpdat$nthres <- tmpdat$C - 1

        tmpdat$Y <- as.matrix(expdats[,8:dim(expdats)[2]])

      } else if(unique(expdats$exp) %in% c("LP2013_e2")){

        # one set of old items, multiple sets of new items

        tmpdat$X_DO <- as.matrix(selectorth %>% select("factor(isold)1"))
        tmpdat$X_DN <-  as.matrix(selectorth %>% select(matches("isold)0")))
        tmpdat$X_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$X_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


        tmpdat$Z_DO <-  as.matrix(selectorth %>% select("factor(isold)1"))
        tmpdat$Z_DN <- as.matrix(selectorth %>% select(matches("isold)0")))
        tmpdat$Z_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$Z_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


        tmpdat$G_DO <- 1
        tmpdat$G_DN <- length(unique(expdats$condition))
        tmpdat$G_go <- 1
        tmpdat$G_neut <- 1
        tmpdat$K_DO <- 1
        tmpdat$K_DN <- length(unique(expdats$condition))
        tmpdat$K_go <- 1
        tmpdat$K_neut <- 1

      } else {

        # fully crossed within:
        tmpdat$X_DO <- as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$X_DN <-  as.matrix(selectorth %>% select(matches("isold)0")))
        tmpdat$X_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$X_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))

        tmpdat$Z_DO <-  as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$Z_DN <- as.matrix(selectorth %>% select(matches("isold)0")))
        tmpdat$Z_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$Z_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))

        tmpdat$G_DO <- length(unique(expdats$condition))
        tmpdat$G_DN <- length(unique(expdats$condition))
        tmpdat$G_go <- 1
        tmpdat$G_neut <- 1

        tmpdat$K_DO <- length(unique(expdats$condition))
        tmpdat$K_DN <- length(unique(expdats$condition))
        tmpdat$K_go <- 1
        tmpdat$K_neut <- 1

      }

    }

    if (unique(expdats$exp) %in% c("LBA2019_e1","LBA2019_e2","HUW2015_e1","LM2020_e1")){
      tmpdat$M_1 <- tmpdat$G_DO + tmpdat$G_DN + tmpdat$G_go + tmpdat$G_neut
    } else {
      tmpdat$M_1 <- tmpdat$G_DO + tmpdat$G_DN + tmpdat$G_go
    }

  } else {



    # replace contrasts for population-level effects
    tmpdat$N_1 <- length(unique(expdats$id))

    tmpdat$J_1 <- as.numeric(factor(expdats$id))

    tmpdat$K <- length(unique(expdats$condition))


    tmpdat$K_disc <- length(unique(expdats$condition))


    tmpdat$K_nonattmu <- length(unique(expdats$condition))
    tmpdat$X <- as.matrix(selectorth)
    tmpdat$X_disc <- as.matrix(selectorth)
    tmpdat$X_nonattmu <- as.matrix(selectorth) # this is MASDT specific
    tmpdat$X_gamma <- matrix(data = 1, nrow =  tmpdat$N , ncol = tmpdat$nthres) # Indgam
    tmpdat$Z_gamma <- matrix(data = 1, nrow =  tmpdat$N , ncol = tmpdat$nthres) # Indgam

    if(unique(expdats$manipulation == "Within")){

      tmpdat$G <- length(unique(expdats$condition))
      tmpdat$G_disc <- length(unique(expdats$condition))
      tmpdat$G_nonattmu <- length(unique(expdats$condition))
      tmpdat$Z <- as.matrix(selectorth)
      tmpdat$Z_disc <- as.matrix(selectorth)
      tmpdat$Z_nonattmu <- as.matrix(selectorth)


    } else if(unique(expdats$manipulation != "Within")) {

      tmpdat$G <- 1
      tmpdat$G_disc <- 1
      tmpdat$G_nonattmu <- 1
      tmpdat$Z <- as.matrix(selectorth[,1])
      tmpdat$Z_disc <- as.matrix(selectorth[,1])
      tmpdat$Z_nonattmu <- as.matrix(selectorth[,1])

    }

    # for Indgam models: M_1 contains also all threshold parameters

    if(model %in% c("UVSDT","DPSDT","M0SDT","DPRMSDT","DPRM2SDT")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc
    } else if(model %in% c("UVgamSDT","DPRMgamSDT","DPRM2gamSDT","DPgamSDT","M0gamSDT")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc + tmpdat$nthres
    } else if(model %in% c("MAgamSDT")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc + tmpdat$G_nonattmu + tmpdat$nthres
    } else if(model %in% c("MASDT")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc + tmpdat$G_nonattmu
    } else if(model %in% c("EVgamSDT")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$nthres
    } else {
      tmpdat$M_1 <- tmpdat$G # EVSDT, Gumbel

    }
  }



  tmpdat$prior_only <- 0
  tmpdat$prior_alpha_mu <- c(10,3,1,0.1)[1:tmpdat$Chalf] # for DPRM, 2HTM
  tmpdat$prior_alpha_scale <- c(2,1,0.5,0.25)[1:tmpdat$Chalf] # for DPRM, 2HTM
  tmpdat$prior_alpha_mu_a <- c(1,2,3,5)[1:tmpdat$Chalf] # for DPRM, 2HTM
  tmpdat$prior_alpha_scale_a <- c(0.25,.5,1,1.5)[1:tmpdat$Chalf] # for DPRM, 2HTM
  tmpdat$N_1 <- length(unique(expdats$id))
  tmpdat$J_1 <- as.numeric(factor(expdats$id))
  tmpdat
}

# prep data for Item models
makeTmpDatItem <- function(expdats,model){

  tmpdat <- NULL

  tmpdat$exp <- unique(expdats$exp)

  if(unique(expdats$manipulation) != "None"){

    if(model == "2HTM"){

      if(unique(expdats$exp) == "GKH1999_e1"){

        newX <-  model.matrix(~0+factor(isold) + factor(isold):condition_freq + factor(isold):condition_freq:condition_studyduration,expdats)

        selectorth <- as_tibble(newX)

      } else {

        newX <-  model.matrix(~0+factor(isold) + factor(isold):condition,expdats)

        selectorth <- as_tibble(newX)
      }

    } else {
      newX <- model.matrix(~isold * condition,expdats)

      selectorth <- as_tibble(newX) %>% select(c("isold",matches("isold:")))
    }

  } else if(unique(expdats$manipulation) == "None"){

    if(model == "2HTM"){

      newX <-model.matrix(~ 0 + isold,expdats %>% mutate(isold = factor(isold)))

      selectorth <- newX


    } else {

      newX <- model.matrix(~isold,expdats)

      selectorth <- as_tibble(newX) %>% select(c("isold",matches("isold:")))
    }

  }


  tmpdat$N <- dim(expdats)[1]

  tmpdat$Y <-expdats$rating


  tmpdat$C <- length(unique(expdats$rating))
  tmpdat$Chalf <- ifelse(model == "DPRM2SDT",2,floor(tmpdat$C / 2))
  tmpdat$nthres <- tmpdat$C - 1



  if(model == "2HTM"){



    if(unique(expdats$manipulation == "None")){

      tmpdat$K_DO <-1
      tmpdat$K_DN <-1
      tmpdat$K_go <-1
      tmpdat$K_neut <-1

      tmpdat$G_DO <- 1
      tmpdat$G_DN <- 1
      tmpdat$G_go <- 1
      tmpdat$G_neut <-1

      tmpdat$X_DO <- as.matrix(as_tibble(selectorth) %>% select("isold1") )
      tmpdat$X_DN <- as.matrix(as_tibble(selectorth) %>% select("isold0") )
      tmpdat$X_go <- as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("isold0") ))[[1]],ncol=1))
      tmpdat$X_neut <- as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("isold0") ))[[1]],ncol=1))
      tmpdat$Z_DO <- as.matrix(as_tibble(selectorth) %>% select("isold1") )
      tmpdat$Z_DN <- as.matrix(as_tibble(selectorth) %>% select("isold0") )
      tmpdat$Z_go <- as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("isold0") ))[[1]],ncol=1))
      tmpdat$Z_neut <- as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("isold0") ))[[1]],ncol=1))

      tmpdat$Z_2_DO <- rep(1,length=tmpdat$N)
      tmpdat$Z_2_DN <- rep(1,length=tmpdat$N)
      tmpdat$Z_2_go <- rep(1,length=tmpdat$N)
    } else if (unique(expdats$manipulation == "Between")){


      tmpdat$X_DO <- as.matrix(selectorth %>% select(matches("isold)1")))
      tmpdat$X_DN <-  as.matrix(selectorth %>% select(matches("isold)0")))
      tmpdat$X_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
      tmpdat$X_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


      tmpdat$Z_DO <-  as.matrix(as_tibble(selectorth) %>% select("factor(isold)1") )
      tmpdat$Z_DN <- as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") )
      tmpdat$Z_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
      tmpdat$Z_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


      tmpdat$G_DO <- 1
      tmpdat$G_DN <- 1
      tmpdat$G_go <- 1
      tmpdat$G_neut <- 1
      tmpdat$K_DO <- length(unique(expdats$condition))
      tmpdat$K_DN <- length(unique(expdats$condition))
      tmpdat$K_go <- 1
      tmpdat$K_neut <- 1

      tmpdat$Z_2_DO <- rep(1,length=tmpdat$N)
      tmpdat$Z_2_DN <- rep(1,length=tmpdat$N)
      tmpdat$Z_2_go <- rep(1,length=tmpdat$N)

    } else if (unique(expdats$manipulation == "Within")){

      if(unique(expdats$exp) %in% c("KY2016_e1","GKH1999_e2",
                                    "GKH1999_e3",
                                    "GKH1999_e4",
                                    "CSR2015_e1",
                                    "FBH2013_e1",
                                    "FO2016_e2a",
                                    "KAWY2013_e3",

                                    "MKG2018_e2",
                                    "RSP2012_e1",
                                    "USS2015_e1",
                                    "WMS2012_e1",
                                    "ZMD2011_e1"

      )){
        head(expdats)
        tmpdat$X_DO <- as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$X_DN <-  as.matrix(selectorth %>% select("factor(isold)0"))
        tmpdat$X_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$X_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


        tmpdat$Z_DO <-  as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$Z_DN <- as.matrix(selectorth %>% select("factor(isold)0"))
        tmpdat$Z_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$Z_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


        tmpdat$G_DO <- length(unique(expdats$condition))
        tmpdat$G_DN <- 1
        tmpdat$G_go <- 1
        tmpdat$G_neut <- 1
        tmpdat$K_DO <- length(unique(expdats$condition))
        tmpdat$K_DN <- 1
        tmpdat$K_go <- 1
        tmpdat$K_neut <- 1

        tmpdat$Z_2_DO <- rep(1,length=tmpdat$N)
        tmpdat$Z_2_DN <- rep(1,length=tmpdat$N)
        tmpdat$Z_2_go <- rep(1,length=tmpdat$N)

      } else if(unique(expdats$exp) == "GKH1999_e1"){

        tmpdat$X_DO <- as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$X_DN <-  as.matrix(selectorth %>% select("factor(isold)0","factor(isold)0:condition_freq1"))
        tmpdat$X_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$X_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


        tmpdat$Z_DO <-  as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$Z_DN <- as.matrix(selectorth %>% select("factor(isold)0","factor(isold)0:condition_freq1"))
        tmpdat$Z_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$Z_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))


        tmpdat$G_DO <- dim(tmpdat$Z_DO)[2]
        tmpdat$G_DN <- dim(tmpdat$Z_DN)[2]
        tmpdat$G_go <- 1
        tmpdat$G_neut <- 1
        tmpdat$K_DO <- dim(tmpdat$X_DO)[2]
        tmpdat$K_DN <- dim(tmpdat$X_DN)[2]
        tmpdat$K_go <- 1
        tmpdat$K_neut <- 1

        tmpdat$Z_2_DO <- rep(1,length=tmpdat$N)
        tmpdat$Z_2_DN <- rep(1,length=tmpdat$N)
        tmpdat$Z_2_go <- rep(1,length=tmpdat$N)

      } else {

        # fully crossed within:
        tmpdat$X_DO <- as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$X_DN <-  as.matrix(selectorth %>% select(matches("isold)0")))
        tmpdat$X_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$X_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))

        tmpdat$Z_DO <-  as.matrix(selectorth %>% select(matches("isold)1")))
        tmpdat$Z_DN <- as.matrix(selectorth %>% select(matches("isold)0")))
        tmpdat$Z_go <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))
        tmpdat$Z_neut <-  as.matrix(rep(1,dim(as.matrix(as_tibble(selectorth) %>% select("factor(isold)0") ))[[1]],ncol=1))

        tmpdat$G_DO <- length(unique(expdats$condition))
        tmpdat$G_DN <- length(unique(expdats$condition))
        tmpdat$G_go <- 1
        tmpdat$G_neut <- 1

        tmpdat$K_DO <- length(unique(expdats$condition))
        tmpdat$K_DN <- length(unique(expdats$condition))
        tmpdat$K_go <- 1
        tmpdat$K_neut <- 1

        tmpdat$Z_2_DO <- rep(1,length=tmpdat$N)
        tmpdat$Z_2_DN <- rep(1,length=tmpdat$N)
        tmpdat$Z_2_go <- rep(1,length=tmpdat$N)

      }

    }

    if (unique(expdats$exp) %in% c("LBA2019_e1","LBA2019_e2","HUW2015_e1","LM2020_e1")){
      tmpdat$M_1 <- tmpdat$G_DO + tmpdat$G_DN + tmpdat$G_go + tmpdat$G_neut
    } else {
      tmpdat$M_1 <- tmpdat$G_DO + tmpdat$G_DN + tmpdat$G_go
    }

  } else {



    # replace contrasts for population-level effects
    tmpdat$N_1 <- length(unique(expdats$id))

    tmpdat$J_1 <- as.numeric(factor(expdats$id))

    tmpdat$K <- length(unique(expdats$condition))


    tmpdat$K_disc <- length(unique(expdats$condition))


    tmpdat$K_nonattmu <- length(unique(expdats$condition))
    tmpdat$X <- as.matrix(selectorth)
    tmpdat$X_disc <- as.matrix(selectorth)
    tmpdat$X_nonattmu <- as.matrix(selectorth) # this is MASDT specific
    tmpdat$X_gamma <- matrix(data = 1, nrow =  tmpdat$N , ncol = tmpdat$nthres)
    tmpdat$Z_gamma <- matrix(data = 1, nrow =  tmpdat$N , ncol = tmpdat$nthres)

    tmpdat$Z_2_mu <- rep(1,length=tmpdat$N)
    tmpdat$Z_2_disc <- rep(1,length=tmpdat$N)
    tmpdat$Z_2_nonattmu <- rep(1,length=tmpdat$N)

    if(unique(expdats$manipulation == "Within")){

      tmpdat$G <- length(unique(expdats$condition))
      tmpdat$G_disc <- length(unique(expdats$condition))
      tmpdat$G_nonattmu <- length(unique(expdats$condition))
      tmpdat$Z <- as.matrix(selectorth)
      tmpdat$Z_disc <- as.matrix(selectorth)
      tmpdat$Z_nonattmu <- as.matrix(selectorth)

      tmpdat$Z_2_mu <- rep(1,length=tmpdat$N)
      tmpdat$Z_2_disc <- rep(1,length=tmpdat$N)
      tmpdat$Z_2_nonattmu <- rep(1,length=tmpdat$N)


    } else if(unique(expdats$manipulation != "Within")) {

      tmpdat$G <- 1
      tmpdat$G_disc <- 1
      tmpdat$G_nonattmu <- 1
      tmpdat$Z <- as.matrix(selectorth[,1])
      tmpdat$Z_disc <- as.matrix(selectorth[,1])
      tmpdat$Z_nonattmu <- as.matrix(selectorth[,1])

      tmpdat$Z_2_mu <- rep(1,length=tmpdat$N)
      tmpdat$Z_2_disc <- rep(1,length=tmpdat$N)
      tmpdat$Z_2_nonattmu <- rep(1,length=tmpdat$N)

    }
    # participant-effects
    if(model %in% c("UVgamSDT","DPRMgamSDT","DPRM2gamSDT","M0gamSDT","DPgamSDT")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc + tmpdat$nthres
    } else if(model %in% c("EVgamSDT")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$nthres
    }else if(model %in% c("MAgamSDT")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc + tmpdat$G_nonattmu + tmpdat$nthres
    }
  }

  # Number of parameters with item effects
  if(model %in% c("UVgamSDT","EVgamSDT")){
    tmpdat$M_2 <- 1
  } else if(model %in% c("DPgamSDT","M0gamSDT")){
    tmpdat$M_2 <- 2
  } else if(model %in% c("MASDT","2HTM","MAgamSDT")){
    tmpdat$M_2 <- 3
  }

  tmpdat$prior_only <- 0
  tmpdat$prior_alpha_mu <- c(10,3,1,0.1)[1:tmpdat$Chalf]
  tmpdat$prior_alpha_scale <- c(2,1,0.5,0.25)[1:tmpdat$Chalf]
  tmpdat$prior_alpha_mu_a <- c(1,2,3,5)[1:tmpdat$Chalf]
  tmpdat$prior_alpha_scale_a <- c(0.25,.5,1,1.5)[1:tmpdat$Chalf]
  tmpdat$N_1 <- length(unique(expdats$id))
  tmpdat$N_2 <- length(unique(expdats$item))
  tmpdat$J_1 <- as.numeric(factor(expdats$id))
  tmpdat$J_2 <- as.numeric(factor(expdats$item))

  tmpdat
}

# Modify standata for item models for KFCV and LOP ------------------
# for Ind and Indgam models for prediction

prepLOP_Ind <- function(tmpdat,testf,modelpred){

  tmpdat_lop <- tmpdat
  tmpdat_lop$Y <- tmpdat$Y[testf,]
  tmpdat_lop$J_1 <- tmpdat$J_1[testf]
  tmpdat_lop$N <- length(tmpdat_lop$J_1)
  tmpdat_lop$N_1 <- length(unique(tmpdat_lop$J_1))

  if(modelpred=="2HTM") {

    tmpdat_lop$X_DO <- as.matrix(tmpdat$X_DO[testf,])
    tmpdat_lop$X_DN <- as.matrix(tmpdat$X_DN[testf,])
    tmpdat_lop$X_go <- as.matrix(tmpdat$X_go[testf,])
    tmpdat_lop$X_neut <- as.matrix(tmpdat$X_neut[testf,])
    tmpdat_lop$Z_DO <- as.matrix(tmpdat$Z_DO[testf,])
    tmpdat_lop$Z_DN <- as.matrix(tmpdat$Z_DN[testf,])
    tmpdat_lop$Z_go <- as.matrix(tmpdat$Z_go[testf,])
    tmpdat_lop$Z_neut <- as.matrix(tmpdat$Z_neut[testf,])

  } else {


    tmpdat_lop$X <- as.matrix(tmpdat$X[testf,])
    tmpdat_lop$X_disc <- as.matrix(tmpdat$X_disc[testf,])
    tmpdat_lop$Z_disc <- as.matrix(tmpdat$Z_disc[testf,])
    tmpdat_lop$X_gamma <- as.matrix(tmpdat$X_gamma[testf,])
    tmpdat_lop$Z_gamma <- as.matrix(tmpdat$Z_gamma[testf,])
    tmpdat_lop$X_nonattmu <- as.matrix(tmpdat$X_nonattmu[testf,])
    tmpdat_lop$Z_nonattmu <- as.matrix(tmpdat$Z_nonattmu[testf,])
    tmpdat_lop$Z <- as.matrix(tmpdat$Z[testf,])


  }

  tmpdat_lop
}

# for item effects models

makeLO_item <- function(exp,tmpdat,fold,model,type,state){


  if (unique(exp$manipulation) == "Within"){
    newX <- model.matrix(~isold * condition,exp)
    selectorth <- as_tibble(newX) %>% select(c("isold",matches("isold:")))
  }

  ## designate which trials are supposed to be included in the data set
  ## type: LOP, KFCV
  ## state: train, (test)
  if(type == "LOP"){

    if(state == "train"){

      train <- exp$LOPFolds != fold

    } else {

      train <- exp$LOPFolds == fold

    }
  } else if (type=="KFCV") {

    if(state == "train"){

      train <- exp$Folds != fold

    } else {

      train <- exp$Folds == fold

    }
  }

  tmpdat_lop <- tmpdat
  tmpdat_lop$N <- length(train[train==TRUE])
  tmpdat_lop$Y <- tmpdat$Y[train]
  tmpdat_lop$isold <- tmpdat$isold[train,]


  if(model=="2HTM"){

    if(unique(exp$manipulation) != "None"){
      tmpdat_lop$X_DO <- as.matrix(tmpdat$X_DO[train,])
      tmpdat_lop$X_DN <- as.matrix(tmpdat$X_DN[train,])
      tmpdat_lop$X_go <- as.matrix(tmpdat$X_go[train,])
      tmpdat_lop$X_neut <- as.matrix(tmpdat$X_neut[train,])

      tmpdat_lop$Z_DO <- as.matrix(tmpdat$Z_DO[train,])
      tmpdat_lop$Z_DN <- as.matrix(tmpdat$Z_DN[train,])
      tmpdat_lop$Z_go <- as.matrix(tmpdat$Z_go[train,])
      tmpdat_lop$Z_neut <- as.matrix(tmpdat$Z_neut[train,])



    } else{
      tmpdat_lop$X_DO <- as.matrix(tmpdat$X_DO[train])
      tmpdat_lop$X_DN <- as.matrix(tmpdat$X_DN[train])
      tmpdat_lop$X_go <- as.matrix(tmpdat$X_go[train])
      tmpdat_lop$X_neut <- as.matrix(tmpdat$X_neut[train])

      tmpdat_lop$Z_DO <- as.matrix(tmpdat$Z_DO[train])
      tmpdat_lop$Z_DN <- as.matrix(tmpdat$Z_DN[train])
      tmpdat_lop$Z_go <- as.matrix(tmpdat$Z_go[train])
      tmpdat_lop$Z_neut <- as.matrix(tmpdat$Z_neut[train])
    }

    tmpdat_lop$Z_2_DO <- tmpdat$Z_2_DO[train]
    tmpdat_lop$Z_2_DN <- tmpdat$Z_2_DN[train]
    tmpdat_lop$Z_2_go <- tmpdat$Z_2_go[train]

    tmpdat_lop$J_1 <- as.integer(as.factor(tmpdat$J_1[train]))
    tmpdat_lop$N_1 <- length(unique(tmpdat_lop$J_1))
    tmpdat_lop$J_2 <- tmpdat$J_2[train]
    tmpdat_lop$N_2 <- length(unique(tmpdat_lop$J_2))
  } else {



    if (unique(exp$manipulation) != "None"){

      tmpdat_lop$X <-  tmpdat$X[train,]

      if(model %in% c("UVSDT","DPSDT","M0SDT","MASDT","DPRMSDT","DPRM2SDT","DPgamSDT","DPRMgamSDT")){
        tmpdat_lop$X_disc <-  tmpdat$X_disc[train,]
        tmpdat_lop$X_nonattmu <-  as.matrix(tmpdat$X_nonattmu[train,])
      }
    } else {

      tmpdat_lop$X <-  as.matrix(tmpdat$X[train])

      if(model %in% c("UVSDT","DPSDT","M0SDT","MASDT","DPRMSDT","DPRM2SDT","DPgamSDT","DPRMgamSDT")){
        tmpdat_lop$X_disc <-  as.matrix(tmpdat$X_disc[train])
        tmpdat_lop$X_nonattmu <-  as.matrix(tmpdat$X_nonattmu[train])
      }

    }


    tmpdat_lop$J_1 <- as.integer(as.factor(tmpdat$J_1[train]))
    tmpdat_lop$N_1 <- length(unique(tmpdat_lop$J_1))
    tmpdat_lop$J_2 <- tmpdat$J_2[train]

    if (unique(exp$manipulation) == "Within"){

      # replace contrasts for group-level effects
      # replace Zs for mu_o
      # replace Zs for sigma_o
      if(model == "UVSDT"){
        tmpdat_lop[grepl("Z_1_", names(tmpdat_lop))] <- c(selectorth[train,],selectorth[train,])

      } else if(model %in% c("EVSDT","Gumbel")){
        tmpdat_lop[grepl("Z_1_", names(tmpdat_lop))] <- c(selectorth[train,])
      }
      tmpdat_lop$Z_2_mu <-tmpdat$Z_2_mu[train]
      tmpdat_lop$Z_2_disc <-tmpdat$Z_2_disc[train]
      tmpdat_lop$Z_2_nonattmu <-tmpdat$Z_2_nonattmu[train]
    } else  {

      tmpdat_lop$Z_1_1 <-tmpdat$Z_1_1[train]
      if(model %in% c("UVSDT","DPSDT","M0SDT","MASDT","DPRMSDT","DPRM2SDT","DPgamSDT","DPRMgamSDT")){
        tmpdat_lop$Z_1_disc_2 <-tmpdat$Z_1_disc_2[train]
      }
      tmpdat_lop$Z_2_mu <-tmpdat$Z_2_mu[train]
      tmpdat_lop$Z_2_disc <-tmpdat$Z_2_disc[train]
      tmpdat_lop$Z_2_nonattmu <-tmpdat$Z_2_nonattmu[train]
    }
    tmpdat_lop$Z <- as.matrix(tmpdat$Z[train,])
    tmpdat_lop$Z_disc <- as.matrix(tmpdat$Z_disc[train,])
    tmpdat_lop$Z_nonattmu <- as.matrix(tmpdat$Z_nonattmu[train,])

    tmpdat_lop$X_gamma <- as.matrix(tmpdat$X_gamma[train,])
    tmpdat_lop$Z_gamma <- as.matrix(tmpdat$Z_gamma[train,])

  }
  tmpdat_lop
}

# initvariables for Gumbel and DPRMSDT -------------------------------

initfun <- function(tmpdat,model){

  if(model == "Gumbel"){
    return (gumbelinit(tmpdat))
  } else if (model %in% c("DPRM2SDT","DPRMSDT")){
    return (dprminit(tmpdat))
  } else {
    return (0.1)
  }

}

init_makevariable <- function(tmpdat){

  mu_cr = sort(runif(tmpdat$nthres,-3,1))
  sigma_cr = runif(tmpdat$nthres,0,0.1)


  if(tmpdat$K == 1){
    b = rnorm(1,-1,0.1)
  } else {
    b = c(rnorm(1,-1,0.1),rnorm(tmpdat$K - 1,0,0.2))
  }

  list(mu_cr = mu_cr,
       sigma_cr = sigma_cr,
       Intercept = t(apply(mvtnorm::rmvnorm(tmpdat$N_1,mu_cr,diag(sigma_cr)),1,sort)),
       b = b
  )

}

dprminit <- function(tmpdat) {

  # thresholds <- init_makevariable(tmpdat)
  # if(tmpdat$K == 1){
  #   b = rnorm(1,1,0.2)
  #
  # } else {
  #   b = c(rnorm(1,1,0.2),rnorm(tmpdat$K - 1,0,0.01))
  #
  # }

  # if(tmpdat$K_disc == 1){
  #
  #   b_disc = rnorm(1,-0.4,0.1)
  # } else {
  #
  #   b_disc = c(rnorm(1,-0.4,0.1),rnorm(tmpdat$K_disc - 1,0,0.01))
  # }
  #
  list(
    # b =b,
    # b_disc =b_disc,

    # sd_1 = runif(tmpdat$M_1, 0.5, 0.65),
    #
    # z_1 = matrix(rnorm(tmpdat$M_1*tmpdat$N_1, 0, 0.01),
    #              tmpdat$M_1, tmpdat$N_1),
    # L_1 = diag(tmpdat$M_1),
    mu_s =  c(runif(1,10,15),runif(1,3,4)),

    s_m = rdirichlet(tmpdat$N_1,c(10,3))
  )


}
gumbelinit <- function(tmpdat) {

  thresholds <- init_makevariable(tmpdat)

  if(is.null(tmpdat$J_2)){

    list(
      b = thresholds$b,
      mu_cr = thresholds$mu_cr,
      sigma_cr = thresholds$sigma_cr,
      Intercept = thresholds$Intercept,
      sd_1 = runif(tmpdat$M_1, 0.01, 0.05),
      #z_1 = tmpdat$X,
      z_1 = matrix(rnorm(tmpdat$M_1*tmpdat$N_1, 0, 0.01),
                   tmpdat$M_1, tmpdat$N_1),
      L_1 = diag(tmpdat$M_1)
    )
  } else {

    list(
      b = thresholds$b,
      mu_cr = thresholds$mu_cr,
      sigma_cr = thresholds$sigma_cr,
      Intercept = thresholds$Intercept,
      sd_2 = runif(tmpdat$M_2, 0.2, 0.3),
      z_2 = matrix(rnorm(tmpdat$M_2*tmpdat$N_2, 0, 0.01),
                   tmpdat$M_2, tmpdat$N_2),
      sd_1 = runif(tmpdat$M_1, 0.45, 0.5),
      z_1 = tmpdat$X,
      L_1 = diag(tmpdat$M_1)
    )
  }
}

#
# # predicrtive p-values -------------------
#
#
# T1statX2 <- function(n,nhat){
#   (n-nhat)^2/nhat
# }
#
# subsetdf_ind <- function(n,J_1,conditions){
#   data.frame(n) %>% mutate(id = J_1,
#                            cond = conditions) %>%
#     group_by(id,cond) %>%
#     summarize_all(sum) %>%
#     select(-cond) %>% ungroup()
# }
#
# subsetdf_exp <- function(n,J_1,conditions){
#
#   data.frame(n) %>% mutate(id = J_1,
#                            cond = conditions) %>%
#     group_by(id,cond) %>%
#     summarize_all(sum) %>%
#     ungroup() %>%
#     select(-id) %>%
#     group_by(cond) %>%
#     summarize_all(mean) %>%
#     select(-cond) %>% ungroup()
#
# }
#
# T1variouslevel <- function(theta,exp,model){
#
#   if(unique(exp$itemInd == "Item")){
#
#     J_1 <- exp$id # participant number per trial
#     conditions <- paste0(exp$condition,"_",exp$isold) # condition per trial (_isold for within-subject/within-block conditions)
#     Y <- exp$rating # ratings per trial
#
#     # ratings per trial in matrix of trial by rating category
#     nobs <- matrix(0,nrow=dim(theta)[1],ncol=max(Y))
#
#     for(i in seq_along(Y)){
#       nobs[i,Y[i]] <- 1
#     }
#
#     # item level T1:
#
#     nhat <- theta
#     nhat[nhat == 0] <- 1e-8
#
#     npred <- extraDistr::rmnom(dim(theta)[1],1,theta)
#
#     T1itempred <- T1statX2(npred,nhat)
#     T1itemobs <- T1statX2(nobs,nhat)
#
#     item_subj <- bind_cols(itemT1preds = rowSums(T1itempred),
#                            itemT1obs = rowSums(T1itemobs),
#                            id = J_1) %>%
#       group_by(id) %>%
#       summarize_all(sum)
#
#
#   } else {
#
#     prepdata <- selectModel(exp,model)
#     expd <- prepdata$tmpdat
#     expd2 <- makeFreqDataInd(exp,model)
#
#     J_1 <- expd2$id # participant number per trial
#     conditions <- paste0(expd2$condition,"_",expd2$isold) # condition per trial (_isold for within-subject/within-block conditions)
#     Y <- expd$Y # ratings per trial
#
#     theta[theta == 0] <- 1e-8
#
#     nhat <- rowSums(Y) * theta
#     npred <- extraDistr::rmnom(dim(theta)[1],rowSums(Y),theta)
#     nobs <- Y
#   }
#
#   # individual level T1:
#
#   npred_ind <-subsetdf_ind(npred,J_1,conditions)%>% nest(npred = !id)
#   nobs_ind <- subsetdf_ind(nobs,J_1,conditions) %>% nest(nobs = !id)
#   nhat_ind <- subsetdf_ind(nhat,J_1,conditions)%>% nest(nhat = !id)
#
#   ind_subj<-left_join(left_join(npred_ind,nobs_ind),nhat_ind) %>%
#     group_by(id) %>%
#     mutate(indT1preds = map2(npred,nhat, ~sum(T1statX2(.x,.y))),
#            indT1obs = map2(nobs,nhat, ~sum(T1statX2(.x,.y)))) %>%
#     select(-npred,-nobs,-nhat) %>%
#     unnest(cols = c(indT1preds, indT1obs))
#
#
#   #exp-level T1
#
#   npred_exp<-subsetdf_exp(npred,J_1,conditions)
#   nobs_exp <- subsetdf_exp(nobs,J_1,conditions)
#   nhat_exp <- subsetdf_exp(nhat,J_1,conditions)
#
#   expT1preds <- sum(T1statX2(npred_exp,nhat_exp))
#   expT1obs <- sum(T1statX2(nobs_exp,nhat_exp))
#
#
#   if(unique(exp$itemInd) == "Item"){
#     T1 <- left_join(item_subj,ind_subj) %>% bind_cols(expT1preds=expT1preds,
#                                                       expT1obs = expT1obs)
#   } else {
#     T1 <- ind_subj %>% bind_cols(expT1preds=expT1preds,
#                                  expT1obs = expT1obs)
#   }
#
#   return(list(T1 = T1,
#               avpredfreq = data.frame(npred_exp)))
#
# }
# T1variouslevelLOP <- function(theta,exp,model,fold){
#
#   if(unique(exp$itemInd == "Item")){
#
#     J_1 <- exp$id # participant number per trial
#     conditions <- paste0(exp$condition,"_",exp$isold) # condition per trial (_isold for within-subject/within-block conditions)
#     Y <- exp$rating # ratings per trial
#
#     # ratings per trial in matrix of trial by rating category
#     nobs <- matrix(0,nrow=dim(theta)[1],ncol=max(Y))
#
#     for(i in seq_along(Y)){
#       nobs[i,Y[i]] <- 1
#     }
#
#     # item level T1:
#
#     nhat <- theta
#     nhat[nhat == 0] <- 1e-8
#
#     npred <- extraDistr::rmnom(dim(theta)[1],1,theta)
#
#     T1itempred <- T1statX2(npred,nhat)
#     T1itemobs <- T1statX2(nobs,nhat)
#
#     item_subj <- bind_cols(itemT1preds = rowSums(T1itempred),
#                            itemT1obs = rowSums(T1itemobs),
#                            id = J_1) %>%
#       group_by(id) %>%
#       summarize_all(sum)
#
#
#   } else {
#
#     prepdata <- selectModel(exp %>% filter(LOPFolds==fold),model)
#     expd <- prepdata$tmpdat
#     expd2 <- makeFreqDataInd(exp,model) %>% filter(LOPFolds==fold)
#
#     J_1 <- expd2$id # participant number per trial
#     conditions <- paste0(expd2$condition,"_",expd2$isold) # condition per trial (_isold for within-subject/within-block conditions)
#     Y <- expd$Y # ratings per trial
#
#     theta[theta == 0] <- 1e-8
#
#     nhat <- rowSums(Y) * theta
#     npred <- extraDistr::rmnom(dim(theta)[1],rowSums(Y),theta)
#     nobs <- Y
#   }
#
#   # individual level T1:
#
#   npred_ind <-subsetdf_ind(npred,J_1,conditions)%>% nest(npred = !id)
#   nobs_ind <- subsetdf_ind(nobs,J_1,conditions) %>% nest(nobs = !id)
#   nhat_ind <- subsetdf_ind(nhat,J_1,conditions)%>% nest(nhat = !id)
#
#   ind_subj<-left_join(left_join(npred_ind,nobs_ind),nhat_ind) %>%
#     group_by(id) %>%
#     mutate(indT1preds = map2(npred,nhat, ~sum(T1statX2(.x,.y))),
#            indT1obs = map2(nobs,nhat, ~sum(T1statX2(.x,.y)))) %>%
#     select(-npred,-nobs,-nhat) %>%
#     unnest(cols = c(indT1preds, indT1obs))
#
#
#   #exp-level T1
#
#   # npred_exp<-subsetdf_exp(npred,J_1,conditions)
#   # nobs_exp <- subsetdf_exp(nobs,J_1,conditions)
#   # nhat_exp <- subsetdf_exp(nhat,J_1,conditions)
#   #
#   # expT1preds <- sum(T1statX2(npred_exp,nhat_exp))
#   # expT1obs <- sum(T1statX2(nobs_exp,nhat_exp))
#   #
#   #
#   # if(unique(exp$itemInd) == "Item"){
#   #   T1 <- left_join(item_subj,ind_subj) %>% bind_cols(expT1preds=expT1preds,
#   #                                                     expT1obs = expT1obs)
#   # } else {
#   T1 <- ind_subj# %>% bind_cols(expT1preds=expT1preds,
#   #              expT1obs = expT1obs)
#   #  }
#
#   return(T1)
#
# }
#
