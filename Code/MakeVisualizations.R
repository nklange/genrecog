# General information ----------------------------------------------------------

# Create visualizations for publications
# more specific/experiment-specific visualizations/diagnostics are in shiny apps.

library(tidyverse)
library(waffle)


# No item effects (79 data sets) -----------------------------------------------------

KFCVlistp <- readRDS("../Predictions/NoItemEffects/processedKFCV.rds")
LOPlistp <- readRDS("../Predictions/NoItemEffects/processedLOP.rds")
Devlistp <- readRDS("../Predictions/NoItemEffects/processedDev.rds")
LOOIClistp <- readRDS("../Predictions/NoItemEffects/processedLOOIC.rds")


allDev <- bind_rows(Devlistp,LOOICs,KFCVlistp,LOPlistp) %>%
  mutate(type = factor(type, levels = c("Full","LOOIC","KFCV","LOP")))



ranges <- allDev %>%
  filter(quantdesc == "cin") %>%
  #filter(experiment %in% fulls) %>%
  #filter(type!="LOP") %>%
  group_by(type,experiment) %>%
  #filter(experiment=="AP2007_e1") %>%
  mutate(topbottom = rep(c("start","end"),n()/2)) %>%
  group_by(type,experiment,topbottom) %>%
  arrange(type,experiment,topbottom,quants) %>%
  mutate(ids = row_number()) %>%
  group_by(type,experiment) %>%
  mutate(refend = quants[1]) %>%
  mutate(simplecomp = ifelse((quants - quants[1]) < 0,"in","out")) %>%
  # mutate(simplecomp2 = ifelse(model %in% c("EVSDT","Gumbel") & simplecomp=="in","simplejoint","unjoint")) %>%
  #  mutate(simplecomp3 = any(simplecomp2 == "simplejoint")) %>%
  group_by(type,experiment,topbottom) %>%
  filter((topbottom == "end" & ids == 1) |
           (topbottom == "start" & ids == 2)) %>%
  ungroup() %>%
  group_by(type,experiment) %>%
  mutate(model = as.character(model)) %>%
  #filter(!experiment %in% c("KL2012_e3","RSP2012_e1","TMP2014_e1","SBT2018_e1","RS2009_e2")) %>%
  mutate(modelwin = ifelse(diff(quants) > 0,"single","undecided")) %>%
  filter(ids==1) %>%
  mutate(win = paste0(model,"_",modelwin))

unique(ranges$win)

full <- bind_rows(
                  ranges  %>%
                    group_by(type,win) %>% count()) %>%
  mutate(model = word(win,1,sep = "\\_"),
         decide = word(win,2,sep = "\\_"))%>%
  mutate(type = factor(type,levels=c("Full","LOOIC","KFCV","LOP"),
                       labels = c("Deviance","LOO-IC","out-of-sample K-Fold","out-of-sample Participant")),

         win = factor(win, levels = c("2HTM_single",
                                      "2HTM_undecided",
                                      "MASDT_single" ,
                                      "MASDT_undecided",
                                      "DPRMSDT_single" ,
                                      "DPRMSDT_undecided",
                                      "M0SDT_single",
                                      "M0SDT_undecided",
                                      "DPSDT_single",
                                      "DPSDT_undecided",
                                      "UVSDT_single",
                                      "UVSDT_undecided" ,
                                      "Gumbel_single" ,
                                      "Gumbel_undecided",
                                      "EVSDT_single",
                                      "EVSDT_undecided"  )

         )) %>%
  arrange(type,win)


restrictedwaffle <- ggplot(full, aes(fill = win, values = n,color=decide)) +
  geom_waffle(n_rows = 8,height=0.9,width=0.9) +
  facet_wrap(.~type)+

  theme_bw()+
  scale_fill_manual(name = NULL,


                    # this is done manually depending on models/winning actually included


                    values = c("#E69F00",alpha("#E69F00",.4), "#56B4E9", alpha("#56B4E9",.4),
                               "#D55E00", alpha("#D55E00",.4),"#F0E442",alpha("#F0E442",.4),
                               alpha("#CC79A7",.4),
                               "#009E73", alpha("#009E73",.4),
                               "#0072B2",alpha("#0072B2",.4),
                               "#999999", alpha("#999999",.4))) +
  coord_equal() +
  scale_color_manual(name = NULL,
                     values=c("black",alpha("white",.9)))+


  theme(axis.title=element_blank(),

        strip.text = element_text(color="black",size=12),
        strip.background = element_rect(fill="transparent"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "None")

legendplot <- tibble(model = rep(c("2HTM",
                                   "MASDT",
                                   "DPRMSDT",
                                   "M0SDT"  ,
                                   "DPSDT",
                                   "UVSDT" ,
                                   "Gumbel",
                                   "EVSDT"),each=2),
                     decide = rep(c("single","undecided"),8)) %>%
  mutate(win = paste0(model,"_",decide)) %>%
  mutate(win = factor(win, levels = c("2HTM_single",
                                      "2HTM_undecided",
                                      "MASDT_single" ,
                                      "MASDT_undecided",
                                      "DPRMSDT_single" ,
                                      "DPRMSDT_undecided",
                                      "M0SDT_single",
                                      "M0SDT_undecided"  ,
                                      "DPSDT_single",
                                      "DPSDT_undecided",
                                      "UVSDT_single",
                                      "UVSDT_undecided" ,
                                      "Gumbel_single" ,
                                      "Gumbel_undecided",
                                      "EVSDT_single",
                                      "EVSDT_undecided"  ))) %>%

  mutate(model = factor(model,levels=c("2HTM","MASDT","DPRMSDT","M0SDT","DPSDT","UVSDT","Gumbel","EVSDT"),
                        labels=c("MPT","MSDT","DPSDT^RM","MSDT^0","DPSDT","UVSDT","Gumbel","EVSDT"))) %>%
  mutate(decide = factor(decide,levels=c("single","undecided"),
                         labels=c("unique","joint")))

# make legend via separate plot as it doesn't work otherwise

legendp <- ggplot(legendplot, aes(fill = win, x=decide, color=decide)) +
  geom_point(y=0,shape=22,size=6) +
  facet_grid(model~.,labeller=label_parsed)+
  scale_fill_manual(name = NULL,



                    values = c("#E69F00",alpha("#E69F00",.4), "#56B4E9", alpha("#56B4E9",.4),
                               "#D55E00", alpha("#D55E00",.4),"#F0E442",alpha("#F0E442",.4),
                               "#CC79A7",alpha("#CC79A7",.4),
                               "#009E73", alpha("#009E73",.4),
                               "#0072B2",alpha("#0072B2",.4),
                               "#999999", alpha("#999999",.4))) +
  #coord_equal() +
  scale_color_manual(name = NULL,
                     values=c("black","white")) +
  scale_x_discrete(position="top",name="winning")+
  scale_y_continuous(limits = c(-0.5,0.5))+
  coord_equal()+
  theme_minimal()+
  theme(legend.position = "None",
        strip.text.y = element_text(angle=0,hjust=0,size=12),
        panel.grid = element_blank(),
        #axis.text.x = element_text(angle=90),
        axis.text.y  = element_blank())



plot_grid(restrictedwaffle,legendp,nrow=1,rel_widths=c(.8,.2))

# Combine all model types (10 data sets) -------------------------------------------------

exps <- c("APH2016_e1","D2007_e1b","DW2020_e1","FGR2019_e1","FO2016_e2a",
          "KFH2013_e1","KY2011_e1","SB2020_e2","SD2014_e1","ZOL2021_e1")
models <- c("EVSDT","UVSDT","M0SDT","MASDT","DPSDT","2HTM")
# no item effects

Ind_allDev <- bind_rows(Devlistp,LOOIClistp,KFCVlistp,LOPlistp) %>%
  mutate(type = factor(type, levels = c("Full","LOOIC","KFCV","LOP"))) %>%
  filter(experiment %in% exps) %>%
  filter(model %in% models) %>%
  mutate(model = factor(model,levels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT"),
                        labels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT"))) %>%
  mutate(datatype = "a")


# no item effects, random crit

KFCVlist_indgam <- readRDS("../Predictions/NoItemEffectsRandomCrit/processedKFCV.rds")
LOPlist_indgam <- readRDS("../Predictions/NoItemEffectsRandomCrit/processedLOP.rds")
Devlist_indgam <- readRDS("../Predictions/NoItemEffectsRandomCrit/processedDev.rds")
LOOICs_indgam <- readRDS("../Predictions/NoItemEffectsRandomCrit/processedLOOIC.rds")

Indgam_allDev <- bind_rows(Devlist_indgam,LOOICs_indgam,KFCVlist_indgam,LOPlist_indgam) %>%
  mutate(type = factor(type, levels = c("Full","LOOIC","KFCV","LOP"))) %>%
  filter(experiment %in% exps) %>%
  mutate(model = factor(model,levels=c("2HTM","MAgamSDT","M0gamSDT","DPgamSDT","UVgamSDT","EVgamSDT"),
                        labels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT")))%>%
  mutate(datatype = "b")

# item effects, random crit

KFCVlist_item <- readRDS("../Predictions/ItemEffectsRandomCrit/processedKFCV.rds")
LOPlist_item <- readRDS("../Predictions/ItemEffectsRandomCrit/processedLOP.rds")
Devlist_item <- readRDS("../Predictions/ItemEffectsRandomCrit/processedDev.rds")
LOOICs_item <- readRDS("../Predictions/ItemEffectsRandomCrit/processedLOOIC.rds")

Item_allDev <- bind_rows(Devlist_item,LOOICs_item,KFCVlist_item,LOPlist_item) %>%
  mutate(type = factor(type, levels = c("Full","LOOIC","KFCV","LOP"))) %>%
  filter(experiment %in% exps)  %>%
  mutate(model = factor(model,levels=c("2HTM","MAgamSDT","M0gamSDT","DPgamSDT","UVgamSDT","EVgamSDT"),
                        labels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT")))%>%
  mutate(datatype = "c")


allDev <- bind_rows(Ind_allDev,
                    Indgam_allDev,
                    Item_allDev)

allDev <- bind_rows(allDevItem,allDevInd,allDevCrit) %>%
  mutate(model = factor(model,levels=c("EVSDT","UVSDT","DPSDT","M0SDT","MASDT","2HTM"))) %>%
  mutate(modeltype = factor(modeltype,levels=c("a","b","c"),
                            labels=c("Part","Part + Crit", "Part + Item + Crit")))

ind1plot <- function(Devs){
  ggplot(Devs,aes(x=quants,y=model))+
    #geom_point()+
    geom_line(data=Devs %>% filter(!(quantdesc %in% c("dmed","Mean","HDI"))),
              aes(x=quants,y=model,group=interaction(quantdesc,model),color=quantdesc),
              size=5)+

    scale_color_manual(values=c("#440154",
                                "#22A884",
                                "#FDE725"),
                       labels=c("95%","80%","50%"))+
    # geom_point(data=Devs %>% filter(quantdesc == "dmed"),
    #            aes(x=quants,y=model),
    #            color="black",
    #            size=2,
    #            position=position_nudge(y = +0.1, x = 0))+
    # geom_line(data=Devs %>% filter(!(quantdesc == c("Mean","dmed","aout","bmid","cin"))),
    #           aes(x = quants,y=model),
    #           size=1,
    #           position = position_nudge(y=-0.1,x=0))+
    geom_point(data=Devs %>% filter(quantdesc == "Mean"),
               aes(x = quants,y=model,
                   fill="Mean with 95% HDI"),
               shape=16,
               size=1)+
    #geom_hline(data = Divpercent,aes(yintercept=model,alpha=fillval),color="white",size=3)+
    scale_alpha_manual(values=c(.8,0)) +
    scale_x_continuous(name = "out-of-sample Deviance (-2 * LL)")+
    labs(color="Mean with HDI")+
    labs(fill = "")+
    # facet_wrap(.~experiment,scales ="free_x")+
    facet_wrap(modeltype~experiment,ncol=1)+
    theme_bw()+
    theme(axis.title.y=element_blank(),
          axis.title.x=element_text(size=6),
          strip.background = element_rect(fill="white"),
          strip.text = element_text(color="black",size=6),
          axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=10),
          legend.position="None",
          panel.grid = element_blank(),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8))
}


indplot <- function(Devs){
  ggplot(Devs,aes(x=quants,y=model))+
    #geom_point()+
    geom_line(data=Devs %>% filter(!(quantdesc %in% c("dmed","Mean","HDI"))),
              aes(x=quants,y=model,group=interaction(quantdesc,model),color=quantdesc),
              size=5)+

    scale_color_manual(values=c("#440154",
                                "#22A884",
                                "#FDE725"),
                       labels=c("95%","80%","50%"))+
    # geom_point(data=Devs %>% filter(quantdesc == "dmed"),
    #            aes(x=quants,y=model),
    #            color="black",
    #            size=2,
    #            position=position_nudge(y = +0.1, x = 0))+
    # geom_line(data=Devs %>% filter(!(quantdesc == c("Mean","dmed","aout","bmid","cin"))),
    #           aes(x = quants,y=model),
    #           size=1,
    #           position = position_nudge(y=-0.1,x=0))+
    geom_point(data=Devs %>% filter(quantdesc == "Mean"),
               aes(x = quants,y=model,
                   fill="Mean with 95% HDI"),
               shape=16,
               size=1)+
    #geom_hline(data = Divpercent,aes(yintercept=model,alpha=fillval),color="white",size=3)+
    scale_alpha_manual(values=c(.8,0)) +
    scale_x_continuous(name = "out-of-sample Deviance (-2 * LL)")+
    labs(color="Mean with HDI")+
    labs(fill = "")+
    # facet_wrap(.~experiment,scales ="free_x")+
    facet_wrap(modeltype~experiment,ncol=1)+
    theme_bw()+
    theme(axis.title.y=element_blank(),
          axis.title.x=element_text(size=6),
          strip.background = element_rect(fill="white"),
          strip.text = element_text(color="black",size=6),
          axis.text.x = element_text(size=6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="None",
          panel.grid = element_blank(),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8))
}


plind<-list()
exps <- c("APH2016_e1","D2007_e1b","DW2020_e1","FGR2019_e1","FO2016_e2a",
          "KFH2013_e1","KY2011_e1","SB2020_e2","SD2014_e1","ZOL2021_e1")


for(i in seq_along(exps)){

  if(i == 1){
    plind[[i]] <- ind1plot(Devs <- allDev %>% filter(type=="LOP") %>% filter(experiment==exps[i]))
  } else {
    plind[[i]] <- indplot(Devs <- allDev %>% filter(type=="LOP") %>% filter(experiment==exps[i]))
  }

}


library(cowplot)
plot_grid(plind[[1]],plind[[2]],
          plind[[3]],plind[[4]],
          plind[[5]],plind[[6]],
          plind[[7]],plind[[8]],
          plind[[9]],plind[[10]],nrow=1,rel_widths=c(.13,rep(.09,9)))


## Waffle plots -------
ranges <- allDev %>%
  filter(quantdesc == "cin") %>%
  #filter(experiment %in% fulls) %>%
  #filter(type!="LOP") %>%
  group_by(datatype,type,experiment) %>%
  #filter(experiment=="AP2007_e1") %>%
  mutate(topbottom = rep(c("start","end"),n()/2)) %>%
  group_by(datatype,type,experiment,topbottom) %>%
  arrange(datatype,type,experiment,topbottom,quants) %>%
  mutate(ids = row_number()) %>%
  group_by(datatype,type,experiment) %>%
  mutate(refend = quants[1]) %>%
  mutate(simplecomp = ifelse((quants - quants[1]) < 0,"in","out")) %>%
  # mutate(simplecomp2 = ifelse(model %in% c("EVSDT","Gumbel") & simplecomp=="in","simplejoint","unjoint")) %>%
  #  mutate(simplecomp3 = any(simplecomp2 == "simplejoint")) %>%
  group_by(datatype,type,experiment,topbottom) %>%
  filter((topbottom == "end" & ids == 1) |
           (topbottom == "start" & ids == 2)) %>%
  ungroup() %>%
  group_by(datatype,type,experiment) %>%
  mutate(model = as.character(model)) %>%
  #filter(!experiment %in% c("KL2012_e3","RSP2012_e1","TMP2014_e1","SBT2018_e1","RS2009_e2")) %>%
  mutate(modelwin = ifelse(diff(quants) > 0,"single","undecided")) %>%
  filter(ids==1) %>%
  mutate(win = paste0(model,"_",modelwin))

unique(ranges$win)

full <- bind_rows(
                  ranges  %>%
                    group_by(datatype,type,win) %>% count()) %>%
  mutate(model = word(win,1,sep = "\\_"),
         decide = word(win,2,sep = "\\_"))%>%
  mutate(datatype = factor(datatype,levels=c("a","b","c"),
                           labels=c("Ppt","Ppt/Crit","Ppt/Crit + Item"))) %>%

  mutate(type = factor(type,levels=c("Full","LOOIC","KFCV","LOP"),
                       labels = c("Deviance","LOO-IC","out-of-sample K-Fold","out-of-sample Participant")),

         win = factor(win, levels = c("2HTM_single",
                                      "2HTM_undecided",
                                      "MASDT_single" ,
                                      "MASDT_undecided",

                                      "M0SDT_single",
                                      "M0SDT_undecided",
                                      "DPSDT_single",
                                      "DPSDT_undecided",
                                      "UVSDT_single",
                                      "UVSDT_undecided" ,

                                      "EVSDT_single",
                                      "EVSDT_undecided"  )

         )) %>%
  arrange(type,win) %>%
  mutate(win = ifelse(type=="KFCV","NA",win))


restrictedwaffle <- ggplot(full, aes(fill = win, values = n,color=decide)) +
  geom_waffle(n_rows = 2,height=0.9,width=0.9) +
  facet_grid(datatype~type,switch="y")+

  theme_bw()+
  scale_fill_manual(name = NULL,



                    values = c("#E69F00",alpha("#E69F00",.4), #2HTM
                               "#56B4E9", alpha("#56B4E9",.4),#MASDT
                               #"#D55E00", alpha("#D55E00",.4),#DPRMSDT
                               "#F0E442",alpha("#F0E442",.4), #M0SDT
                               "#CC79A7", #DPSDT
                               "#009E73", alpha("#009E73",.4), #UVSDT
                               #"#0072B2",alpha("#0072B2",.4), #Gumbel
                               "#999999",alpha("#999999",.4))) + #EVSDT
  coord_equal() +
  scale_color_manual(name = NULL,
                     values=c("black",alpha("white",.9)))+


  theme(axis.title=element_blank(),

        strip.text = element_text(color="black",size=12),
        strip.background = element_rect(fill="transparent"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "None")



legendplot <- tibble(model = rep(c("2HTM",
                                   "MASDT",
                                   # "DPRMSDT",
                                   "M0SDT"  ,
                                   "DPSDT",
                                   "UVSDT" ,
                                   #"Gumbel",
                                   "EVSDT"),each=2),
                     decide = rep(c("single","undecided"),6)) %>%
  mutate(win = paste0(model,"_",decide)) %>%
  mutate(win = factor(win, levels = c("2HTM_single",
                                      "2HTM_undecided",
                                      "MASDT_single" ,
                                      "MASDT_undecided",
                                      # "DPRMSDT_single" ,
                                      #  "DPRMSDT_undecided",
                                      "M0SDT_single",
                                      "M0SDT_undecided"  ,
                                      "DPSDT_single",
                                      "DPSDT_undecided",
                                      "UVSDT_single",
                                      "UVSDT_undecided" ,
                                      # "Gumbel_single" ,
                                      #  "Gumbel_undecided",
                                      "EVSDT_single",
                                      "EVSDT_undecided"  ))) %>%

  mutate(model = factor(model,levels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT"),
                        labels=c("MPT","MSDT","MSDT^0","DPSDT","UVSDT","EVSDT"))) %>%
  mutate(decide = factor(decide,levels=c("single","undecided"),
                         labels=c("unique","joint")))

legendp <- ggplot(legendplot, aes(fill = win, x=decide, color=decide)) +
  geom_point(y=0,shape=22,size=6) +
  facet_grid(model~.,labeller=label_parsed)+
  scale_fill_manual(name = NULL,



                    values = c("#E69F00",alpha("#E69F00",.4),
                               "#56B4E9", alpha("#56B4E9",.4),
                               #"#D55E00", alpha("#D55E00",.4),
                               "#F0E442",alpha("#F0E442",.4),
                               "#CC79A7",alpha("#CC79A7",.4),
                               "#009E73", alpha("#009E73",.4),
                               #"#0072B2",alpha("#0072B2",.4),
                               "#999999", alpha("#999999",.4))) +
  #coord_equal() +
  scale_color_manual(name = NULL,
                     values=c("black","white")) +
  scale_x_discrete(position="top",name="winning")+
  scale_y_continuous(limits = c(-0.5,0.5))+
  coord_equal()+
  theme_minimal()+
  theme(legend.position = "None",
        strip.text.y = element_text(angle=0,hjust=0,size=12),
        panel.grid = element_blank(),
        #axis.text.x = element_text(angle=90),
        axis.text.y  = element_blank())

unique(full$win)

plot_grid(restrictedwaffle,legendp,nrow=1,rel_widths=c(.85,.15))


# By Individual (10 data sets) ---------------------------------------------------------
# probably for supplementary material / sanity checking
# LOP


DevCrit <- readRDS("../Predictions/NoItemEffectsRandomCrit/processedIndividualLOP.rds") %>%
  mutate(model = factor(model,levels=c("2HTM","MAgamSDT","M0gamSDT","DPgamSDT","UVgamSDT","EVgamSDT"),
                        labels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT"))) %>%
  mutate(type="b")


DevInd <- readRDS("../Predictions/NoItemEffects/processedIndividualLOP.rds") %>%
  mutate(model = factor(model,levels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT"),
                        labels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT"))) %>%
  mutate(type = "a")

DevItem <- readRDS("../Predictions/ItemEffectsRandomCrit/processedIndividualLOP.rds") %>%
  mutate(model = factor(model,levels=c("2HTM","MAgamSDT","M0gamSDT","DPgamSDT","UVgamSDT","EVgamSDT"),
                        labels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT"))) %>%
  mutate(type = "c")

Devs <- bind_rows(DevCrit,DevInd,DevItem) %>%
  ungroup() %>%
  #filter(model %in% c("2HTM","UVSDT")) %>%
  mutate(type = factor(type,labels=c("Ppt","Ppt/Crit","Ppt/Crit + Item")))


allDev <- Devs

ranges <- allDev %>%
  filter(quantdesc == "cin") %>%
  #filter(experiment %in% fulls) %>%
  #filter(type!="LOP") %>%
  group_by(type,experiment,id) %>%
  #filter(experiment=="AP2007_e1") %>%
  mutate(topbottom = rep(c("start","end"),n()/2)) %>%
  group_by(type,experiment,id,topbottom) %>%
  arrange(type,experiment,id,topbottom,quants) %>%
  mutate(ids = row_number()) %>%
  group_by(type,experiment,id) %>%
  mutate(refend = quants[1]) %>%
  mutate(simplecomp = ifelse((quants - quants[1]) < 0,"in","out")) %>%
  #mutate(simplecomp2 = ifelse(model %in% c("EVSDT","Gumbel") & simplecomp=="in","simplejoint","unjoint")) %>%
  #mutate(simplecomp3 = any(simplecomp2 == "simplejoint")) %>%
  group_by(type,experiment,id,topbottom) %>%
  filter((topbottom == "end" & ids == 1) |
           (topbottom == "start" & ids == 2)) %>%
  ungroup() %>%
  group_by(type,experiment,id) %>%
  mutate(model = as.character(model)) %>%
  #filter(!experiment %in% c("KL2012_e3","RSP2012_e1","TMP2014_e1","SBT2018_e1","RS2009_e2")) %>%
  mutate(modelwin = ifelse(diff(quants) > 0,"single","undecided")) %>%
  filter(ids==1) %>%
  mutate(win = paste0(model,"_",modelwin))


full <- bind_rows(
  ranges  %>%
    group_by(type,experiment,win) %>% count()) %>%
  mutate(model = word(win,1,sep = "\\_"),
         decide = word(win,2,sep = "\\_")) %>%
  mutate(#type = factor(type,levels=c("Participant","Participant + Crit","Participant + Crit + Item")),
    win = factor(win, levels=c("2HTM_single",
                               "2HTM_undecided",
                               "MASDT_single",
                               "MASDT_undecided",
                               "M0SDT_single",
                               "M0SDT_undecided",
                               "DPSDT_single",
                               "DPSDT_undecided",
                               "UVSDT_single",
                               "UVSDT_undecided",
                               "EVSDT_single",
                               "EVSDT_undecided"))) %>%
  arrange(type,experiment,win)


restrictedwaffleExpLOP <- ggplot(full, aes(fill = win, values = n,color=win)) +
  geom_waffle( n_rows = 8,height=0.8,width=0.8) +
  facet_grid(type~experiment,switch="y")+
  scale_fill_manual(name = NULL,
                    #                   # labels=rev(c("EVSDT","UVSDT","DPSDT",
                    #                   #              expression(paste(M[0],"SDT")),"MSDT","2HTM")),
                    values = c("#E69F00",alpha("#E69F00",.6),
                               alpha("#56B4E9",.6),
                               #"#D55E00", alpha("#D55E00",.4),
                               alpha("#F0E442",.6),
                               alpha("#CC79A7",.6),
                               alpha("#009E73",.6),
                               #"#0072B2",alpha("#0072B2",.4),
                               "#999999",alpha("#999999",.6)))+
  # scale_fill_manual(name = NULL,
  #                   #                   # labels=rev(c("EVSDT","UVSDT","DPSDT",
  #                   #                   #              expression(paste(M[0],"SDT")),"MSDT","2HTM")),
  #                   values = c("#E69F00","#56B4E9","#009E73","#0072B2","#F0E442",
  #                              "#CC79A7"))+
  coord_equal() +
  scale_color_manual(name = NULL,
                     values=c("black",alpha("white",.9),alpha("white",.9),alpha("white",.9),
                              alpha("white",.9),alpha("white",.9),"black",alpha("white",.9)))+
  # scale_color_manual(name = NULL,
  #                    values=c("black",alpha("white",.9)))+


  theme_bw()+
  theme(axis.title=element_blank(),
        aspect.ratio=1,
        strip.background.x = element_rect(fill="transparent",color="transparent"),
        strip.background.y = element_rect(fill="white",color="black"),
        strip.text = element_text(color="black",size=10),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position= "None")


plot_grid(restrictedwaffleExpLOP,legendp,nrow=1,rel_widths=c(.9,.1))

# Full


DevCrit <- readRDS("D:/HS01/Predictions/NoItemEffectsRandomCrit/processedIndividualFull.rds") %>%
  mutate(model = factor(model,levels=c("2HTM","MAgamSDT","M0gamSDT","DPgamSDT","UVgamSDT","EVgamSDT"),
                        labels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT"))) %>%
  mutate(type="b")


DevInd <- readRDS("D:/HS01/Predictions/NoItemEffects/processedIndividualFull.rds") %>%
  mutate(model = factor(model,levels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT"),
                        labels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT"))) %>%
  mutate(type = "a")

DevItem <- readRDS("D:/HS01/Predictions/ItemEffectsRandomCrit/processedIndividualFull.rds") %>%
  mutate(model = factor(model,levels=c("2HTM","MAgamSDT","M0gamSDT","DPgamSDT","UVgamSDT","EVgamSDT"),
                        labels=c("2HTM","MASDT","M0SDT","DPSDT","UVSDT","EVSDT"))) %>%
  mutate(type = "c")

Devs <- bind_rows(DevCrit,DevInd,DevItem) %>%
  ungroup() %>%
  #filter(model %in% c("2HTM","UVSDT")) %>%
  mutate(type = factor(type,labels=c("Ppt","Ppt/Crit","Ppt/Crit + Item"))) %>%
  rename("experiment" = exp)


allDev <- Devs

ranges <- allDev %>%
  filter(quantdesc == "cin") %>%
  #filter(experiment %in% fulls) %>%
  #filter(type!="LOP") %>%
  group_by(type,experiment,id) %>%
  #filter(experiment=="AP2007_e1") %>%
  mutate(topbottom = rep(c("start","end"),n()/2)) %>%
  group_by(type,experiment,id,topbottom) %>%
  arrange(type,experiment,id,topbottom,quants) %>%
  mutate(ids = row_number()) %>%
  group_by(type,experiment,id) %>%
  mutate(refend = quants[1]) %>%
  mutate(simplecomp = ifelse((quants - quants[1]) < 0,"in","out")) %>%
  #mutate(simplecomp2 = ifelse(model %in% c("EVSDT","Gumbel") & simplecomp=="in","simplejoint","unjoint")) %>%
  #mutate(simplecomp3 = any(simplecomp2 == "simplejoint")) %>%
  group_by(type,experiment,id,topbottom) %>%
  filter((topbottom == "end" & ids == 1) |
           (topbottom == "start" & ids == 2)) %>%
  ungroup() %>%
  group_by(type,experiment,id) %>%
  mutate(model = as.character(model)) %>%
  #filter(!experiment %in% c("KL2012_e3","RSP2012_e1","TMP2014_e1","SBT2018_e1","RS2009_e2")) %>%
  mutate(modelwin = ifelse(diff(quants) > 0,"single","undecided")) %>%
  filter(ids==1) %>%
  mutate(win = paste0(model,"_",modelwin))


full <- bind_rows(
  ranges  %>%
    group_by(type,experiment,win) %>% count()) %>%
  mutate(model = word(win,1,sep = "\\_"),
         decide = word(win,2,sep = "\\_")) %>%
  mutate(#type = factor(type,levels=c("Participant","Participant + Crit","Participant + Crit + Item")),
    win = factor(win, levels=c("2HTM_single",
                               "2HTM_undecided",
                               "MASDT_single",
                               "MASDT_undecided",
                               "M0SDT_single",
                               "M0SDT_undecided",
                               "DPSDT_single",
                               "DPSDT_undecided",
                               "UVSDT_single",
                               "UVSDT_undecided",
                               "EVSDT_single",
                               "EVSDT_undecided"))) %>%
  arrange(type,experiment,win)
unique(full$win)

library(waffle)

restrictedwaffleExpFull <- ggplot(full, aes(fill = win, values = n,color=win)) +
  geom_waffle( n_rows = 8,height=0.8,width=0.8) +
  facet_grid(type~experiment,switch="y")+
  scale_fill_manual(name = NULL,
                    #                   # labels=rev(c("EVSDT","UVSDT","DPSDT",
                    #                   #              expression(paste(M[0],"SDT")),"MSDT","2HTM")),
                    values = c("#E69F00",alpha("#E69F00",.6),
                               "#56B4E9", alpha("#56B4E9",.6),
                               #"#D55E00", alpha("#D55E00",.4),
                               "#F0E442",alpha("#F0E442",.6),
                               "#CC79A7",alpha("#CC79A7",.6),
                               "#009E73", alpha("#009E73",.6),
                               #"#0072B2",alpha("#0072B2",.4),
                               alpha("#999999",.6)))+
  # scale_fill_manual(name = NULL,
  #                   #                   # labels=rev(c("EVSDT","UVSDT","DPSDT",
  #                   #                   #              expression(paste(M[0],"SDT")),"MSDT","2HTM")),
  #                   values = c("#E69F00","#56B4E9","#009E73","#0072B2","#F0E442",
  #                              "#CC79A7"))+
  coord_equal() +
  scale_color_manual(name = NULL,
                     values=c("black",alpha("white",.9),
                              "black",alpha("white",.9),
                              "black",alpha("white",.9),
                              "black", alpha("white",.9),
                              "black", alpha("white",.9),
                              alpha("white",.9)))+
  # scale_color_manual(name = NULL,
  #                    values=c("black",alpha("white",.9)))+


  theme_bw()+
  theme(axis.title=element_blank(),
        aspect.ratio=1,
        strip.background.x = element_rect(fill="transparent",color="transparent"),
        strip.background.y = element_rect(fill="white",color="black"),
        strip.text = element_text(color="black",size=10),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position= "None")


plot_grid(restrictedwaffleExpFull,legendp,nrow=1,rel_widths=c(.9,.1))

