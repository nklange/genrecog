#library(shinydashboard)
#library(plotly)

library(DT)
library(tidyverse)
library(cowplot)
library(tidybayes)
library(shinythemes)
library(shinyWidgets)
library(shinyjs)
library("scales")
library(zip)

rds_file <- list.files("ProcessedData/", full.names = TRUE)
rds_file3 <- list.files("DownloadData/", full.names=TRUE)

d1 <- purrr::map(rds_file[grepl(".rds",rds_file)], readRDS)

d3 <- purrr::map(rds_file3[grepl(".rds",rds_file3)], readRDS)

explist <-  d1 %>% purrr::map(. %>% ungroup() %>% distinct(exp)) %>% bind_rows() %>% .$exp

explistfull <- d3 %>% purrr::map(. %>% ungroup() %>% distinct(exp)) %>% bind_rows() %>% .$exp

expinfo <- read.csv("ROCexpinformation.csv")


exps1 <-c("AP2007_e1",   "AP2007_e2",   "AP2007_e3",   "APH2016_e1",  "BKSR2013_e1", "BSG2014_e1",
          "BTL2013_e1",  "CSR2015_e1",  "D2007_e1a",   "D2007_e1b",   "DR2012_e1b",  "DW2020_e1",
          "FBH2013_e1",  "FGR2019_e1",  "FO2016_e2a",  "FO2016_e2b",  "FO2016_e3",   "GKH1999_e1",
          "GKH1999_e2",  "GKH1999_e3",  "GKH1999_e4",  "HDM2006_e1",  "HDM2006_e2",  "HUW2015_e1",
          "JCD2012_e1a", "JCD2012_e1b", "JCM2019_e1",  "JW2019_e1" ,  "JWH2009_e1",  "KAWY2013_e2",
          "KAWY2013_e3", "KAWY2013_e4", "KFH2013_e1",  "KK2015_e1",   "KL2012_e1",   "KL2012_e2",
          "KL2012_e3",   "KL2012_e4",   "KUO2017_e2",  "KY2010_e1",   "KY2011_e1",   "KY2016_e1",
          "LBA2019_e1",  "LBA2019_e2",  "LM2020_e1",   "LP2006_e2",   "LP2013_e1",   "LP2013_e2",
          "MKG2018_e2",
          "NFR2013_e1",  "NFR2013_e2",  "NFR2013_e3",  "OBD2017_e1",  "OZH2010_e1",  "PCM2006_e1",
          "PRM2010_e1",  "QGM2021_e1",  "RS2009_e1",   "RS2009_e2",   "RSM2009_e1",  "RSP2012_e1",
          "SB2020_e1",   "SB2020_e2",   "SB2020_e3",   "SBT2018_e1",  "SCR2019_e1",  "SD2004_e2",
          "SD2014_e1",   "SD2014_e2",   "SHJ2005_e1",  "TMP2014_e1",  "TR2017_e1",   "USS2015_e1",
          "WHD2018_e1",  "WKH2020_e1",  "WKS2020_e1",  "WMS2012_e1",  "ZMD2011_e1",  "ZOL2021_e1" )





facet_strip_bigger <- function(gp, size){
  if(missing(gp)){
    print("this function needs a facet_wrap ggplotly object")
  }
  if(missing(size)){
    print("this function needs 'size' argument to be specified as integer. 80 will be introduced as default")
    size <- 80
  }

  n_facets <- c(1:length(gp[["x"]][["layout"]][["shapes"]]))

  for(i in n_facets){
    if(n_facets[i] %% 2 == 0){
      gp[["x"]][["layout"]][["shapes"]][[i]][["y0"]] <- + as.numeric(size)
      gp[["x"]][["layout"]][["shapes"]][[i]][["y1"]] <- 0
    }
  }

  return(gp)
}

makeROC <- function(exp){

  if(is.null(exp$condition_true)){


    if(unique(exp$exp) %in% c(c("CSR2015_e1","FO2016_e2a",
                                "KY2016_e1","MKG2018_e2",
                                "MKG2018_e2","JW2019_e1","USS2015_e1",
                                "ZMD2011_e1","KAWY2013_e3","KAWY2013_e2",
                                "KAWY2013_e4","FBH2013_e1",
                                "WMS2012_e1","RSP2012_e1"),
                              c("BSG2014_e1","KFH2013_e1",
                                "LP2013_e2","MKG2019_e1",
                                "SBT2018_e1","RS2009_e2",
                                "GKH1999_e1","GKH1999_e2",
                                "GKH1999_e3","GKH1999_e4"))){


      expOld <- exp %>% filter(isold == 1)
      expNew <- exp %>% filter(isold == 0)

      expInd <- list()
      for(i in seq_along(unique(expOld$condition))){

        expInd[[i]] <- expOld %>% filter(condition == unique(expOld$condition)[[i]]) %>%
          bind_rows(expNew) %>% mutate(condition = unique(expOld$condition)[[i]])
      }

      ex <- bind_rows(expInd) %>%
        ungroup() %>%
        select(exp,id,condition,isold,rating)


      # } else if (unique(exp$exp) %in% c("BSG2014_e1","KFH2013_e1",
      #              "LP2013_e2","MKG2019_e1",
      #              "SBT2018_e1","RS2009_e2")){
      #
      #
      #   expOld <- exp %>% filter(isold == 1)
      #   expNew <- exp %>% filter(isold == 0)
      #
      #   expInd <- list()
      #   for(i in seq_along(unique(expOld$condition))){
      #
      #     expInd[[i]] <- expOld %>% filter(condition == unique(expOld$condition)[[i]]) %>%
      #       bind_rows(expNew) %>% mutate(condition = unique(expOld$condition)[[i]])
      #   }
      #
      #   ex <- bind_rows(expInd) %>%
      #     ungroup() %>%
      #     select(exp,id,condition,isold,rating)

    } else {



      ex <- exp %>%
        mutate(rating = as.numeric(as.character(rating))) %>%
        ungroup() %>%
        select(exp,id,condition,isold,rating)

    }

  } else {



    expOld <- exp %>% filter(isold == 1)
    expNew <- exp %>% filter(isold == 0)

    expInd <- list()
    for(i in seq_along(unique(expOld$condition))){

      expInd[[i]] <- expOld %>% filter(condition == unique(expOld$condition)[[i]]) %>%
        bind_rows(expNew) %>% mutate(condition = unique(expOld$condition)[[i]])
    }

    ex <- bind_rows(expInd) %>%
      ungroup() %>%
      select(exp,id,condition,isold,rating)

  }



  makerocdata <- ex %>%
    group_by(exp,id,condition,isold,rating) %>%
    summarize(cumresp = length(rating)) %>%
    mutate(rating = factor(rating, levels = c(max(ex$rating):1))) %>%
    arrange(exp,id,condition,isold,rating) %>%
    mutate(isold = factor(isold,levels=c(0,1),labels=c("INew","IOld"))) %>%
    group_by(exp,id,condition,isold) %>%
    complete(rating,fill = list(cumresp = 0)) %>%
    summarize(cumprop = cumsum(cumresp/sum(cumresp)),
              rating = rating) %>%
    ungroup() %>%
    pivot_wider(names_from="isold",values_from="cumprop",
                id_cols = c("exp","id","condition","rating")) %>%
    filter(rating != 1)


  ggplot(makerocdata,aes(x = INew,y=IOld,group=id,text=id))+
    facet_wrap(.~condition,ncol=4)+

    #geom_point(alpha=0.2)+
    geom_line(alpha=0.2)+
    geom_abline(intercept=0,slope=1,color="blue",linetype="dashed")+
    theme_bw()+
    coord_fixed(xlim=c(0,1),ylim=c(0,1))+
    #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
    scale_shape_manual(values=c(16,17),name="Data")+
    scale_x_continuous(name="p('old'|new)")+
    scale_y_continuous(name = "p('old'|old)")+
    theme_bw()+
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14),
          strip.background = element_rect(fill="black"),
          strip.text = element_text(color="white",face = "bold",size=16))
}

makeHead <- function(exp){

   as.data.frame(d3[which(explistfull == exp)])
}
