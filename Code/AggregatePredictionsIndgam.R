# General information ------------------------

## Aggregate the diagnostics and predictions across data sets and models for
## No Item Effects, Random Crit models
## resultant files in ProcessedPredictions/NoItemEffectsRandomCrit/ folder


library(tidyverse)
library(tidybayes) # for mean_hdi function

rdsfile <- list.files("../ProcessedData/", full.names = TRUE)
d1 <- purrr::map(rdsfile[grepl(".rds",rdsfile)], readRDS)
explist4 <-  d1 %>% purrr::map(. %>% ungroup() %>% distinct(exp)) %>% bind_rows() %>% .$exp

# EXTRACT DEVIANCE/OUT-OF-SAMPLE DEVIANCE --------------

EVSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/EVSDT/", full.names = TRUE)
UVSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/UVSDT/", full.names = TRUE)
DPSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/DPSDT/", full.names = TRUE)
M0SDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/M0SDT/", full.names = TRUE)
MASDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/MASDT/", full.names = TRUE)
# technically these data/predictions already exist in the NoItemEffects section but this was makes
# aggregation a bit easier
H2TMlist <- list.files("../Predictions/NoItemEffects/Full/H2TM/", full.names = TRUE)

# Deviance/Full ------------------
# extract deviances from file for all models / save as separate processed predictions

for(i in seq_along(explist4)){

  exp<-unique(d1[[i]]$exp)

  for(model in c("EVgamSDT","UVgamSDT","M0gamSDT","2HTM","DPgamSDT")){


    if(model == "EVgamSDT"){
      extract_expn <- EVSDTlist[grepl(exp,EVSDTlist)]
    } else if(model == "UVgamSDT"){
      extract_expn <- UVSDTlist[grepl(exp,UVSDTlist)]
    } else if(model == "DPgamSDT"){
      extract_expn <- DPSDTlist[grepl(exp,DPSDTlist)]
    } else if(model == "2HTM"){
      extract_expn <- H2TMlist[grepl(exp,H2TMlist)]
    } else if(model == "M0gamSDT"){
      extract_expn <- M0SDTlist[grepl(exp,M0SDTlist)]
    } else if(model == "MAgamSDT"){
      extract_expn <- MASDTlist[grepl(exp,MASDTlist)]
    }


    divs <- purrr::map(extract_expn[grepl(model,extract_expn)], readRDS)# %>% bind_rows()


    print(exp)
    if(length(divs) > 0){

      Deviance <- lapply(divs,function(x) extractDevs(x))


      saveRDS(Deviance[[1]] %>% mutate(experiment = exp) %>%  mutate(model = ifelse(model=="DPRM2SDT","DPRMSDT",model))
              ,file=paste0("../ProcessedPredictions/NoItemEffectsRandomCrit/Full/deviances_",model,"_",exp,"_Full.rds"))
    }
  }

}

# Deviance/Full by ID

EVSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/IndividualFull/", full.names = TRUE)

for(i in seq_along(d1)){

  exp<-unique(d1[[i]]$exp)


  for(model in c("EVgamSDT","UVgamSDT","M0gamSDT","2HTM","DPgamSDT")){


    extract_expn <- EVSDTlist[grepl(exp,EVSDTlist)]
    extract_mod <- extract_expn[grepl(model,extract_expn)]


    divs <- purrr::map(extract_mod[grepl(model,extract_mod)], readRDS)# %>% bind_rows()


    print(exp)
    if(length(divs) > 0){

      Deviance <- lapply(divs,function(x) extractDevsInd(x))


      saveRDS(Deviance[[1]]# %>% mutate(experiment = exp)
              ,file=paste0("../ProcessedPredictions/NoItemEffectsRandomCrit/IndividualFull/deviances_",model,"_",exp,"_Fullindividual.rds"))
    }
  }

}

# KFCV -----------------------------------------

EVSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/KFCV/EVSDT/", full.names = TRUE)
UVSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/KFCV/UVSDT/", full.names = TRUE)
DPSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/KFCV/DPSDT/", full.names = TRUE)
M0SDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/KFCV/M0SDT/", full.names = TRUE)
H2TMlist <- list.files("../Predictions/NoItemEffects/KFCV/H2TM/", full.names = TRUE)
MASDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/KFCV/MASDT/", full.names = TRUE)

for(i in seq_along(d1)){

  exp <-unique(d1[[i]]$exp)

  for(model in c("EVgamSDT","UVgamSDT","M0gamSDT","2HTM","DPgamSDT")){

    if(model == "EVgamSDT"){
      extract_expn <- EVSDTlist[grepl(exp,EVSDTlist)]
    } else if(model == "UVgamSDT"){
      extract_expn <- UVSDTlist[grepl(exp,UVSDTlist)]
    } else if(model == "DPgamSDT"){
      extract_expn <- DPSDTlist[grepl(exp,DPSDTlist)]
    } else if(model == "2HTM"){
      extract_expn <- H2TMlist[grepl(exp,H2TMlist)]
    } else if(model == "M0gamSDT"){
      extract_expn <- M0SDTlist[grepl(exp,M0SDTlist)]
    } else if(model == "MAgamSDT"){
      extract_expn <- MASDTlist[grepl(exp,MASDTlist)]
    }




    divs <- purrr::map(extract_expn[grepl(modeln,extract_expn)], readRDS)# %>% bind_rows()

    print(exp)
    if(length(divs) > 0){



      DevianceInd <-  extractDevs_KFCV(divs[[1]]) %>%
        group_by(experiment,type,id,model,Fold) %>%
        mutate(nume = row_number()) %>%
        group_by(experiment,type,id,model,nume) %>%
        summarize(sumsDev = sum(Dev,na.rm=T)) %>%
        rename("Dev" = sumsDev) %>%
        ungroup() %>%
        select(-nume)

      Deviance <-  extractDevs_KFCV(divs[[1]]) %>%
        group_by(experiment,type,id,model,Fold) %>%
        mutate(nume = row_number()) %>%
        group_by(experiment,type,model,nume) %>%
        summarize(sumsDev = sum(Dev,na.rm=T)) %>%
        rename("Dev" = sumsDev) %>%
        ungroup() %>%
        select(-nume)



      saveRDS(Deviance %>% mutate(model = ifelse(model=="DPRM2SDT","DPRMSDT",model))#%>% mutate(experiment = exp)
              ,file=paste0("../ProcessedPredictions/NoItemEffectsRandomCrit/KFCV/deviances_",model,"_",exp,"_KFCV.rds"))
      saveRDS(DevianceInd %>% mutate(model = ifelse(model=="DPRM2SDT","DPRMSDT",model))#%>% mutate(experiment = exp)
              ,file=paste0("../ProcessedPredictions/NoItemEffectsRandomCrit/IndividualKFCV/deviances_",model,"_",exp,"_individualKFCV.rds"))
    }
  }

}


# LOP ---------------------------------------------------

# across all ppts LOP

EVSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/LOP/EVSDT/", full.names = TRUE)
MASDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/LOP/MASDT/", full.names = TRUE)
UVSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/LOP/UVSDT/", full.names = TRUE)
DPSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/LOP/DPSDT/", full.names = TRUE)
M0SDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/LOP/M0SDT/", full.names = TRUE)
H2TMlist<- list.files("../Predictions/NoItemEffects/LOP/H2TM/", full.names = TRUE)

for(i in seq_along(d1)){

  exp<-unique(d1[[i]]$exp)
  for(model in c("EVgamSDT","UVgamSDT","M0gamSDT","2HTM","DPgamSDT")){

    if(model == "EVgamSDT"){
      extract_expn <- EVSDTlist[grepl(exp,EVSDTlist)]
    } else if(model == "UVgamSDT"){
      extract_expn <- UVSDTlist[grepl(exp,UVSDTlist)]
    } else if(model == "DPgamSDT"){
      extract_expn <- DPSDTlist[grepl(exp,DPSDTlist)]
    } else if(model == "2HTM"){
      extract_expn <- H2TMlist[grepl(exp,H2TMlist)]
    } else if(model == "M0gamSDT"){
      extract_expn <- M0SDTlist[grepl(exp,M0SDTlist)]
    } else if(model == "MAgamSDT"){
      extract_expn <- MASDTlist[grepl(exp,MASDTlist)]
    }




      divs <- purrr::map(extract_expn[grepl(modeln,extract_expn)], readRDS)# %>% bind_rows()





      Deviance <-  extractDevs_LOP(divs[[1]]) %>%
        bind_rows() %>%
        group_by(experiment,type,model,Fold) %>%
        mutate(nume = row_number()) %>%
        group_by(experiment,type,model,nume) %>%
        summarize(sumsDev = sum(Dev,na.rm=T)) %>%
        rename("Dev" = sumsDev) %>%
        ungroup() %>%
        select(-nume)


      saveRDS(Deviance
              ,file=paste0("../ProcessedPredictions/NoItemEffectsRandomCrit/LOP/deviances_",model,"_",exp,"_LOP.rds"))
    }

}

# Individual LOP

for(i in c(4,10,12,14,15,33,41,62,67,78)){



  exp<-unique(d1[[i]]$exp)

    for(model in c("EVgamSDT","UVgamSDT","DPgamSDT","2HTM","M0gamSDT","MAgamSDT")){

    if(model == "EVgamSDT"){
      extract_expn <- EVSDTlist[grepl(exp,EVSDTlist)]
    } else if(model == "UVgamSDT"){
      extract_expn <- UVSDTlist[grepl(exp,UVSDTlist)]
    } else if(model == "DPgamSDT"){
      extract_expn <- DPSDTlist[grepl(exp,DPSDTlist)]
    } else if(model == "2HTM"){
      extract_expn <- H2TMlist[grepl(exp,H2TMlist)]
    } else if(model == "M0gamSDT"){
      extract_expn <- M0SDTlist[grepl(exp,M0SDTlist)]
    } else if(model == "MAgamSDT"){
      extract_expn <- MASDTlist[grepl(exp,MASDTlist)]
    }
    if(length(extract_expn) > 0){




      divs <- purrr::map(extract_expn[grepl(modeln,extract_expn)], readRDS)# %>% bind_rows()



      Deviance <-  extractDevs_LOPind(divs[[1]],d1[[i]])


      saveRDS(Deviance %>% mutate(model = ifelse(model=="DPRM2SDT","DPRMSDT",model))
              ,file=paste0("../ProcessedPredictions/NoItemEffectsRandomCrit/IndividualLOP/deviances_",model,"_",exp,"_LOPind.rds"))
    }
  }

}



# Summary statistics for deviance -----------------------------------

# For plots and shiny apps, the deviance distributions are summarized
# and aggregated in a single file for each fitting type (Full, KFCV, LOP)
# to minimize file sizes

# Full/Deviance

Devlist <-list.files("../ProcessedPredictions/NoItemEffectsRandomCrit/Full", full.names = TRUE)

dats <- purrr::map(Devlist[grepl("deviances",Devlist)], readRDS) %>% bind_rows()

# extract mean with high density interval at .95,.8,.5

hdilevels <- c(.95,.8,.5)
meanhdi <- NULL
for(i in hdilevels){

  onemean <- dats %>% group_by(experiment,model) %>%
    mean_hdi(Dev,.width=i) %>%

    pivot_longer(cols=c("Dev",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
    mutate(quantdesc = factor(quantdesc,levels=c("Dev",".lower",".upper"),
                              labels=c("Mean",paste((1-i)/2),paste(i + ((1-i)/2))))) %>%
    select(-.width,-.point,-.interval) %>% mutate(quantdesc = as.character(quantdesc))

  meanhdi <- meanhdi %>% bind_rows(onemean)

}

forplot <- meanhdi %>% distinct() %>%
  mutate(quantdesc = factor(quantdesc,levels=c("0.025","0.1","0.25","Mean","0.75","0.9","0.975"),
                            labels=c("aout","bmid","cin","Mean","cin","bmid","aout"))) %>%
  mutate(type="Full")

saveRDS(forplot,file=paste0("../Predictions/NoItemEffectsRandomCrit/processedDev.rds"))

# Full/Deviance - by individual ID

Devlist <-list.files("../ProcessedPredictions/NoItemEffectsRandomCrit/IndividualFull", full.names = TRUE)


dats <- purrr::map(Devlist[grepl("deviances",Devlist)], readRDS) %>% bind_rows()

hdilevels <- c(.95,.8,.5)
meanhdi <- NULL
for(i in hdilevels){

  onemean <- dats %>% group_by(exp,id,model) %>%
    mean_hdi(Dev,.width=i) %>%

    pivot_longer(cols=c("Dev",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
    mutate(quantdesc = factor(quantdesc,levels=c("Dev",".lower",".upper"),
                              labels=c("Mean",paste((1-i)/2),paste(i + ((1-i)/2))))) %>%
    select(-.width,-.point,-.interval) %>% mutate(quantdesc = as.character(quantdesc))

  meanhdi <- meanhdi %>% bind_rows(onemean)

}


forplot <- meanhdi %>% distinct() %>%
  mutate(quantdesc = factor(quantdesc,levels=c("0.025","0.1","0.25","Mean","0.75","0.9","0.975"),
                            labels=c("aout","bmid","cin","Mean","cin","bmid","aout"))) %>%
  mutate(type="Full")

saveRDS(forplot,file=paste0("../Predictions/NoItemEffectsRandomCrit/processedIndividualFull.rds"))

# KFCV

Devlist <-list.files("../ProcessedPredictions/NoItemEffectsRandomCrit/KFCV", full.names = TRUE)

dats <- purrr::map(Devlist[grepl("deviances",Devlist)], readRDS) %>% bind_rows()

# extract mean with high density interval at .95,.8,.5

hdilevels <- c(.95,.8,.5)
meanhdi <- NULL
for(i in hdilevels){

  onemean <- dats %>% group_by(experiment,model) %>%
    mean_hdi(Dev,.width=i) %>%

    pivot_longer(cols=c("Dev",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
    mutate(quantdesc = factor(quantdesc,levels=c("Dev",".lower",".upper"),
                              labels=c("Mean",paste((1-i)/2),paste(i + ((1-i)/2))))) %>%
    select(-.width,-.point,-.interval) %>% mutate(quantdesc = as.character(quantdesc))

  meanhdi <- meanhdi %>% bind_rows(onemean)

}

forplot <- meanhdi %>% distinct() %>%
  mutate(quantdesc = factor(quantdesc,levels=c("0.025","0.1","0.25","Mean","0.75","0.9","0.975"),
                            labels=c("aout","bmid","cin","Mean","cin","bmid","aout"))) %>%
  mutate(type="KFCV")

saveRDS(forplot,file=paste0("../Predictions/NoItemEffects/processedKFCV.rds"))

# LOP

Devlist <-list.files("../ProcessedPredictions/NoItemEffectsRandomCrit/LOP", full.names = TRUE)

dats <- purrr::map(Devlist[grepl("deviances",Devlist)], readRDS) %>% bind_rows()

# extract mean with high density interval at .95,.8,.5

hdilevels <- c(.95,.8,.5)
meanhdi <- NULL
for(i in hdilevels){

  onemean <- dats %>% group_by(experiment,model) %>%
    mean_hdi(Dev,.width=i) %>%

    pivot_longer(cols=c("Dev",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
    mutate(quantdesc = factor(quantdesc,levels=c("Dev",".lower",".upper"),
                              labels=c("Mean",paste((1-i)/2),paste(i + ((1-i)/2))))) %>%
    select(-.width,-.point,-.interval) %>% mutate(quantdesc = as.character(quantdesc))

  meanhdi <- meanhdi %>% bind_rows(onemean)

}

forplot <- meanhdi %>% distinct() %>%
  mutate(quantdesc = factor(quantdesc,levels=c("0.025","0.1","0.25","Mean","0.75","0.9","0.975"),
                            labels=c("aout","bmid","cin","Mean","cin","bmid","aout"))) %>%
  mutate(type="Full")

saveRDS(forplot,file=paste0("../Predictions/NoItemEffectsRandomCrit/processedLOP.rds"))

# LOP - by individual ID

Devlist <-list.files("../ProcessedPredictions/NoItemEffectsRandomCrit/IndividualFull", full.names = TRUE)


dats <- purrr::map(Devlist[grepl("deviances",Devlist)], readRDS) %>% bind_rows()

hdilevels <- c(.95,.8,.5)
meanhdi <- NULL
for(i in hdilevels){

  onemean <- dats %>% group_by(exp,id,model) %>%
    mean_hdi(Dev,.width=i) %>%

    pivot_longer(cols=c("Dev",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
    mutate(quantdesc = factor(quantdesc,levels=c("Dev",".lower",".upper"),
                              labels=c("Mean",paste((1-i)/2),paste(i + ((1-i)/2))))) %>%
    select(-.width,-.point,-.interval) %>% mutate(quantdesc = as.character(quantdesc))

  meanhdi <- meanhdi %>% bind_rows(onemean)

}


forplot <- meanhdi %>% distinct() %>%
  mutate(quantdesc = factor(quantdesc,levels=c("0.025","0.1","0.25","Mean","0.75","0.9","0.975"),
                            labels=c("aout","bmid","cin","Mean","cin","bmid","aout"))) %>%
  mutate(type="Full")

saveRDS(forplot,file=paste0("../Predictions/NoItemEffectsRandomCrit/processedIndividualLOP.rds"))



#EXTRACT DIVERGENT TRANSITIONS/pps ------------------------
# Deviance/Full --------------------------------------------
EVSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/EVSDT/", full.names = TRUE)
UVSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/UVSDT/", full.names = TRUE)
DPSDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/DPSDT/", full.names = TRUE)
M0SDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/M0SDT/", full.names = TRUE)
MASDTlist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/MASDT/", full.names = TRUE)
H2TMlist <- list.files("../Predictions/NoItemEffects/Full/H2TM/", full.names = TRUE)

for(i in c(30,59,37,64,77,6,58,34,33,18,70,60)){#seq_along(d1)[1:78]){

  exp<-unique(d1[[i]]$exp)


  for(model in c("EVgamSDT","UVgamSDT","DPgamSDT","2HTM","M0gamSDT","MAgamSDT")){

    if(model == "EVgamSDT"){
      extract_expn <- EVSDTlist[grepl(exp,EVSDTlist)]
    } else if(model == "UVgamSDT"){
      extract_expn <- UVSDTlist[grepl(exp,UVSDTlist)]
    } else if(model == "DPgamSDT"){
      extract_expn <- DPSDTlist[grepl(exp,DPSDTlist)]
    } else if(model == "2HTM"){
      extract_expn <- H2TMlist[grepl(exp,H2TMlist)]
    } else if(model == "M0gamSDT"){
      extract_expn <- M0SDTlist[grepl(exp,M0SDTlist)]
    } else if(model == "MAgamSDT"){
      extract_expn <- MASDTlist[grepl(exp,MASDTlist)]
    }

    print(exp)
    print(model)

    divs <- purrr::map(extract_expn[grepl(exp,extract_expn)], readRDS)# %>% bind_rows()

    #if(!is.null(divs)){
    if(length(divs) > 0){

      T1s <- lapply(divs,function(x) extractT1s_ind(x)) %>% bind_rows()
      Divergent <- lapply(divs,function(x) extractDiv(x)) %>% bind_rows()


      divergtable <- T1s %>%
        mutate(pitemT1 = 1,
               pindT1 = prop_indT1,
               partraw = n_indT1greater,
               partall =n_exp,
               pexpT1 = prop_expT1) %>%

        select(model,pitemT1,pindT1,pexpT1,partraw,partall) %>%
        #select(model,Div) %>%
        pivot_longer(names_to = "type",values_to = "value",
                     cols = c(pitemT1,pindT1,pexpT1,partraw,partall)) %>%

        mutate(experiment = exp)%>%
        mutate(model = model)

      divergtable2 <- Divergent %>%

        mutate(Div = Div * 100) %>%
        # mutate(Div = paste0(round(Div,2),"%")) %>%
        select(model,Div) %>%
        #select(model,Div) %>%
        pivot_longer(names_to = "type",values_to = "value",
                     cols = c(Div)) %>%

        mutate(experiment = exp) %>%
        mutate(model = model)

      # divergent transitions
      saveRDS(divergtable,file=paste0("../ProcessedPredictions/NoItemEffectsRandomCrit/Full/divergtable_",model,"_",exp,"_Full.rds"))
      # ppps
       saveRDS(divergtable2,file=paste0("../ProcessedPredictions/NoItemEffectsRandomCrit/Full/diverg2table_",model,"_",exp,"_Full.rds"))
    }
  }

}

# Summarize percentage divergent transitions

filelist <-list.files("D:/HS01/ProcessedPredictions/NoItemEffectsRandomCrit/Full/", full.names = TRUE)

dats <- purrr::map(filelist[grepl("diverg2table",filelist)], readRDS) %>% bind_rows()

saveRDS(dats %>% filter(type == "Div"),file=paste0("../Predictions/NoItemEffectsRandomCrit/processedDev_Div.rds"))


# KFCV -----------------------

Divs1 <- NULL
for (i in seq_along(d1)){

  exp<-unique(d1[[i]]$exp)
  for(model in c("EVgamSDT","UVgamSDT","DPgamSDT","2HTM","M0gamSDT","MAgamSDT")){

    if(model == "2HTM"){
      modeln1 = "H2TM"
    } else {
      modeln1 = model
    }

    KFCVlist <-list.files(paste0("../Predictions/NoItemEffectsRandomCrit/KFCV/",modeln1,"/"), full.names = TRUE)
    LOPS <- KFCVlist[grepl(exp,KFCVlist)]


    kfcvs<- purrr::map(LOPS, readRDS) #%>% bind_rows()

    if(length(kfcvs)>0){
      DivKFCV <-  lapply(kfcvs,function(x) extractDiv_KFCV(x))  %>%
        bind_rows() %>%
        group_by(experiment,model) %>%
        summarize(value = mean(Div,na.rm=T)) %>%
        ungroup() %>%
        mutate(type = "DivKFCV") %>%
        mutate(value = value * 100) %>%
        mutate(value = paste0(round(value,2),"%"))


      Divs1 <- Divs1 %>% bind_rows(DivKFCV) %>% mutate(model = ifelse(model=="DPRM2SDT","DPRMSDT",model))
    }
  }

}

saveRDS(Divs1,file=paste0("../Predictions/NoItemEffectsRandomCrit/processedKFCV_Div.rds"))

# LOP -----------------------


Divs1 <- NULL
for (i in seq_along(d1)){

  exp<-unique(d1[[i]]$exp)
  for(model in c("EVgamSDT","UVgamSDT","DPgamSDT","2HTM","M0gamSDT","MAgamSDT")){


    if(model == "2HTM"){
      modeln1 = "H2TM"
    } else {
      modeln1 = model
    }


    LOPlist <-list.files(paste0("../Predictions/NoItemEffectsRandomCrit/LOP/",modeln1,"/"), full.names = TRUE)
    LOPS <- LOPlist[grepl(exp,LOPlist)]

    lops<- purrr::map(LOPS[grepl(model,LOPS)], readRDS) #%>% bind_rows()

    if(length(lops)>0){
      DivLOP <-  lapply(lops,function(x) extractDiv_LOP(x[[1]]))  %>%
        bind_rows() %>%
        group_by(experiment,model) %>%
        summarize(value = mean(Div,na.rm=T)) %>%
        ungroup() %>%
        mutate(type = "DivLOP") %>%
        mutate(value = value * 100) %>%
        mutate(value = paste0(round(value,2),"%"))


      Divs1 <- Divs1 %>% bind_rows(DivLOP)%>% mutate(model = ifelse(model=="DPRM2SDT","DPRMSDT",model))
    }
  }

}

saveRDS(Divs1,file=paste0("../Predictions/NoItemEffectsRandomCrit/processedLOP_Div.rds"))



# Summarize posterior p-values -----------------
# only calculated for Deviance/Full fits

filelist <-list.files("../ProcessedPredictions/NoItemEffectsRandomCrit/Full/", full.names = TRUE)

dats <- purrr::map(filelist[grepl("divergtable",filelist)], readRDS) %>% bind_rows()

pps2 <- dats %>%
  filter(type != "Div") %>%
  group_by(experiment,model) %>%
  distinct() %>%

  pivot_wider(id_cols = c("experiment","model"),names_from="type",values_from="value") %>%
  mutate(pitemT1 = "N/A") %>%
  mutate(pindT1 = paste0(ifelse( round(pindT1,2) > 0,sub("^0+", "",  round(pindT1,2)), 0) , " (",partraw,"/",partall,")")) %>%
  mutate(pexpT1 = ifelse( round(pexpT1,2) > 0,sub("^0+", "",  round(pexpT1,2)), "0")) %>%
  select(-partraw,-partall)

saveRDS(pps2,"../Predictions/NoItemEffectsRandomCrit/processedDev_ppp.rds")

# Summarize LOOIC ------------------------------
# for Full/Deviance fits

# Deviance/Full by ID

filelist <- list.files("../Predictions/NoItemEffectsRandomCrit/Full/LOOIC/", full.names = TRUE)

LOOIC <- NULL
pareto <- NULL

for(i in seq_along(d1)){

  exp<-unique(d1[[i]]$exp)


  for(model in c("EVgamSDT","UVgamSDT","DPgamSDT","2HTM","M0gamSDT","MAgamSDT")){


    extract_expn <- filelist[grepl(exp,filelist)]



    divs <- purrr::map(extract_expn[grepl(model,extract_expn)], readRDS)


    print(exp)
    if(length(divs) > 0){

      looic <- extractLOOIC(divs %>% .[[1]])


      LOOIC <- LOOIC %>% bind_rows(looic$LOOIC)
      pareto <- pareto %>% bind_rows(looic$pareto_k)

          }
  }

}

saveRDS(LOOIC
        ,file=paste0("../Predictions/NoItemEffectsRandomCrit/processedLOOIC.rds"))
saveRDS(pareto
        ,file=paste0("../Predictions/NoItemEffectsRandomCrit/processedLOOIC_paretok.rds"))



# Helper functions -----------------------------
# PROCESS DEVIANCES
extractDevs <- function(div){

  Devs <- tibble(Dev = div$deviances$Dev) %>%
    mutate(model = div$information$model,
           type = div$information$type
    )
}


extractDevs_LOP <- function(div){

  divtot <- list()
  for(i in c(1:length(div))){

    if(dim(div[[i]])[[2]] > 8){



      divtot[[i]] <- tibble(Dev = rowSums(div[[i]] %>% select(starts_with("X")))) %>%
        mutate(model = div[[i]]$model,
               experiment = div[[i]]$experiment,
               type = div[[i]]$type,
               Fold = div[[i]]$Fold)

    } else {

      divtot[[i]] <- tibble(Dev = rowSums(div[[i]] %>% select(starts_with("matrix")))) %>%
        mutate(model = div[[i]]$model,
               experiment = div[[i]]$experiment,
               type = div[[i]]$type,
               Fold = div[[i]]$Fold)

    }
  }

  divtot %>% bind_rows()

}
extractDevs_KFCV <- function(div){

  Devs <- tibble(Dev = div$deviances$Dev) %>%
    mutate(model = div$information$model,
           experiment = div$information$experiment,
           type = div$information$type,
           Fold = div$deviances$Folds)


  return(Devs)

}
extractT1s_ind<- function(div){

  div<-divs[[1]]

  id1 <-div$T1s[[1]]$id[1]
  #unique(div$T1s %>% bind_rows() %>% .$id)[1]

  gb <- div$T1s %>%
    bind_rows()

  T1s <- gb %>%
    mutate(
      indT1 = indT1preds - indT1obs) %>%
    group_by(id) %>%
    summarize(
      p_indT1 = mean(indT1 > 0)) %>%
    ungroup() %>%
    summarize(
      meanp_ind = mean(p_indT1),

      sdp_ind = sd(p_indT1),

      n_indT1greater = length(p_indT1[p_indT1 > .05]),
      n_exp = length(p_indT1)) %>%
    mutate(model = div$information$model) %>%
    mutate(
      prop_indT1 = n_indT1greater/n_exp)


  T1sexp <- gb %>%
    mutate(expT1 = expT1preds - expT1obs) %>%
    filter(id == id1) %>%
    select(id,expT1) %>%
    summarize(prop_expT1 = length(expT1[expT1 > 0])/length(expT1))

  return(T1s %>% bind_cols(T1sexp))

}


# Process divergent transitions

extractDiv <- function(div){

  tibble(Div = unique(div$deviances$DivChainsPercent)) %>%
    mutate(model = div$information$model,
           experiment = div$information$experiment
    )
}

extractDiv_KFCV <- function(div){

  tibble(Div = div$deviances$DivChainsPercent) %>%
    mutate(model = div$information$model,
           experiment = div$information$experiment,
           Fold = div$deviances$Folds
    ) %>%
    distinct()
}

extractDiv_LOP <- function(div){

  lapply(div, function(x)  tibble(Div = div$DivChainsPercent) %>%
           mutate(model = div$model,
                  experiment = div$experiment,
                  Fold = div$Fold) %>%
           distinct()) %>%
    bind_rows()


}

# Process LOOIC

extractLOOIC <- function(LOOIC){

# diagnostics - information sampling blubs

partoks <- LOOIC$LOOIC$diagnostics$pareto_k

pareto_k <- tibble(pareto_labels = c("good","ok","bad","verybad"),
       values = c(length(paretoks[paretoks <= 0.5]),
                  length(paretoks[paretoks > 0.5 & paretoks <= .7]),
                  length(paretoks[paretoks > .7 & paretoks <= 1]),
                  length(paretoks[paretoks > 1])),
       percentage = round(values/sum(values) * 100,2))

pareto_out <- bind_cols(pareto_k,LOOIC$information)

# could also extract p_loo to get an idea of
# misspecification  vs flexibility when pareto-k high
# https://discourse.mc-stan.org/t/a-quick-note-what-i-infer-from-p-loo-and-pareto-k-values/3446/7
# would then take some parameter calculations specifically for all data sets


# LOOIC
# using the same labels as for the summary of deviance
# assuming normal distribution here, using the SE to calculate
# the same 'bounds' we do for deviance/out-of sample deviance
# i.e. inner 50%, 80%, 95%
# for an approximation of the same scale (use high-density interval for deviance, not equal-density though)
# (we use inner 50% to determine, whether models are equally good at predicting data,
# alternative would be to use SE scale (+/- 1 SE ~ inner 68%)


looic <- tibble(quantdesc = c("Mean","SE","cin","cin","bmid","bmid","aout","aout"),
                quants = c(as.numeric(LOOIC$LOOIC$estimates["looic",]["Estimate"]),
                           as.numeric(LOOIC$LOOIC$estimates["looic",]["SE"]),
                           as.numeric(LOOIC$LOOIC$estimates["looic",]["Estimate"]) -
                             abs(as.numeric(LOOIC$LOOIC$estimates["looic",]["SE"]) * qnorm(.25)),
                           as.numeric(LOOIC$LOOIC$estimates["looic",]["Estimate"]) +
                             abs(as.numeric(LOOIC$LOOIC$estimates["looic",]["SE"]) * qnorm(.25)),
                           as.numeric(LOOIC$LOOIC$estimates["looic",]["Estimate"]) -
                             abs(as.numeric(LOOIC$LOOIC$estimates["looic",]["SE"]) * qnorm(.1)),
                           as.numeric(LOOIC$LOOIC$estimates["looic",]["Estimate"]) +
                             abs(as.numeric(LOOIC$LOOIC$estimates["looic",]["SE"]) * qnorm(.1)),
                           as.numeric(LOOIC$LOOIC$estimates["looic",]["Estimate"]) -
                             abs(as.numeric(LOOIC$LOOIC$estimates["looic",]["SE"]) * qnorm(.025)),
                           as.numeric(LOOIC$LOOIC$estimates["looic",]["Estimate"]) +
                             abs(as.numeric(LOOIC$LOOIC$estimates["looic",]["SE"]) * qnorm(.025))),
                paretowarning = ifelse(sum(length(paretoks[paretoks > .7 & paretoks <= 1]),
                                         length(paretoks[paretoks > .7 & paretoks <= 1])) >= 1, "warning","ok" ))
LOOICout <- bind_cols(looic,LOOIC$information)


list(pareto_k = pareto_out,
     LOOIC = LOOICout)
}


makeFreqDataInd <- function(exp,model){

  if(model == "2HTM" &
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

extractDevs_LOPind <- function(div,exp){



  divtot <- list()
  for(i in c(1:length(div))){
    uniqueids <- unique(makeFreqDataInd(exp %>% filter(LOPFolds==i),"EVSDT") %>% .$id)
    if(dim(div[[i]])[[2]] > 8){

      divtot[[i]] <- div[[i]] %>%
        pivot_longer(cols = starts_with("X"),values_to = "value",names_to = "id") %>%
        mutate(id = str_remove(id,"X")) %>%
        mutate(id = factor(id,labels = uniqueids)) %>%
        mutate(it = rep(c(1:6000),each=length(uniqueids)))




    } else {



      divtot[[i]] <-  div[[i]] %>%
        pivot_longer(cols = starts_with("matrix"),values_to = "value",names_to = "id") %>%
        mutate(id =1) %>%
        mutate(id = factor(id,labels = uniqueids)) %>%
        mutate(it = rep(c(1:6000),each=length(uniqueids)))

    }
  }

  divtot %>% bind_rows()

}
