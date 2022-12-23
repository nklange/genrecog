library("rstan")
library("tidybayes")
library("brms")
library("cmdstanr")
library("tidyverse")
library("bayestestR")
options(contrasts=c('contr.orthonorm', 'contr.poly'))
rds_file_pred <- list.files("ProcessedData/", full.names = TRUE)
rds_file_pred <- rds_file_pred[grepl(".rds",rds_file_pred)]
explist_pred <- str_remove(str_remove(rds_file_pred,"ProcessedData/"),".rds")
exps <- explist_pred[!(explist_pred %in% c("SSW2015_e3"))]


d1pred <- purrr::map(rds_file_pred, readRDS)


#source("LOPThetaInd.R")
#source("MakeAdjustStanModel.R")
#source("MakeAdjustIndStanModel.R")

pathfit <- "D:/HS01/Fits/ItemEffectsRandomCrit/Full/"


  summarizeparameters <- function(expstring,model,coll){




  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  if(model!="2HTM"){
  coef <- ncoef - nthres
  } else {
    coef <- ncoef
  }
  if(unique(tidydat$manipulation) == "Between"){
    condlabels<-"None"
  } else {
    condlabels <-unique(tidydat$condition)

  }

  if(model %in% c("EVgamSDT","Gumbel")){

    modelpar1 <- c("mu","crit")

    modelpar <- c(rep(modelpar1[[1]],coef),rep(modelpar1[[2]],nthres))

    allSD <- summarizesds(coll)
    participantsd <-  allSD$participantsd %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar = factor(par,labels = modelpar))
    grouplevel_item_sd <- allSD$itemsd %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar = factor(par,labels=modelpar1[!modelpar1=="crit"]))


    grouplevel_noncrit_sd <- participantsd %>% filter(specpar !="crit") %>%
      mutate(condition = parse_number(par)) %>%
      mutate(condition = factor(condition,labels = condlabels))

   grouplevel_crit_sd <- participantsd %>% filter(specpar =="crit") %>%
     mutate(condition = parse_number(par)) %>%
     mutate(condition = factor(condition,labels = c(paste0("cr",c(1:nthres)))))

   population_weights <- bind_rows(summarizeweights(coll[[2]]$b, "b") %>%
                                     mutate(par = as.character(par)) %>%
                                     mutate(condition =parse_number(par)) %>%
                                     mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                     mutate(specpar = modelpar1[[1]]))

   crit_weights <- summarizeweights(coll[[2]]$b_gamma,"b_gamma") %>%
     mutate(par = as.character(par)) %>%
     mutate(specpar =modelpar1[[2]])

   itemeffects <- summarizeitemeffects(coll,tidydat)


   distpar_all <- summarizedistpar(coll,tidydat,model)

   correlationsSDalpha <- summarizecorrelations(coll,coll[[2]]$L_1,"L_1")


   list(grouplevel_noncrit_sd = grouplevel_noncrit_sd,
        grouplevel_crit_sd = grouplevel_crit_sd,
        grouplevel_item_sd = grouplevel_item_sd,
        population_weights = population_weights,
        crit_weights = crit_weights,
        individual_items = itemeffects,
        individual_distpar = distpar_all$individual_dist,
        population_distpar = distpar_all$population_dist,
        individual_crit = distpar_all$individual_crit,
        population_crit = distpar_all$population_crit,
        summary_correlations_ppt = correlationsSDalpha$summary_corr,
        raw_correlations_ppt = correlationsSDalpha$raw_corr,
        Exp = expstring)

  }else if (model %in% c("UVgamSDT")){

    modelpar1 <- c("mu","disc","crit")

    modelpar <- c(rep(modelpar1[[1]],coef/2),rep(modelpar1[[2]],coef/2),rep(modelpar1[[3]],nthres))

    allSD <- summarizesds(coll)
    participantsd <-  allSD$participantsd %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar = factor(par,
                              levels = paste0("sd",c(1:length(modelpar))),
                              labels = modelpar))


      grouplevel_item_sd <- allSD$itemsd %>%
        mutate(par = as.character(par)) %>%
        mutate(specpar = factor(par,labels=modelpar1[modelpar1=="mu"]))



      seqb <-  coef/length(modelpar1[modelpar1!="crit"])
      bound <- seq(1,coef,seqb)
      #bound <- ifelse(bound[2]==2 | length(bound) == 1,1,bound[2])

      if(max(bound)==2){
        grouplevel_noncrit_sd <- participantsd %>% filter(specpar !="crit") %>%
          mutate(par = as.character(par)) %>%
          mutate(condition_tmp = parse_number(par)) %>%
          mutate(specpar = case_when(condition_tmp == 1 ~ modelpar[[1]],
                                     condition_tmp == 2 ~ modelpar[[2]],
                                     TRUE ~ "NA")) %>%
          mutate(condition = condlabels) %>%
          select(-condition_tmp)
      } else {

        grouplevel_noncrit_sd <- participantsd %>% filter(specpar !="crit")  %>%
          mutate(par = as.character(par)) %>%
          mutate(condition_tmp = parse_number(par)) %>%
          mutate(specpar = case_when(condition_tmp < bound[[2]] ~ modelpar[[1]],
                                     condition_tmp >= bound[[2]] ~ modelpar[[2]],
                                     TRUE ~ "NA")) %>%
          mutate(condition = case_when(condition_tmp < bound[[2]] ~ condition_tmp,
                                       condition_tmp >= bound[[2]] ~ condition_tmp - length(condlabels) ,
                                       TRUE ~ 0)) %>%
          mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%

          select(-condition_tmp)
      }


    grouplevel_crit_sd <- participantsd %>% filter(specpar =="crit") %>%
      mutate(condition = parse_number(par)) %>%
      mutate(condition = factor(condition,labels = c(paste0("cr",c(1:nthres)))))

    population_weights <- bind_rows(summarizeweights(coll[[2]]$b, "b") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar1[[1]]),
                                    summarizeweights(coll[[2]]$b_disc, "b_disc") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar1[[2]]))

    crit_weights <- summarizeweights(coll[[2]]$b_gamma,"b_gamma") %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar =modelpar1[[2]])

    # itemeffects_mu <-summarizeitemeffects(coll,coll[[2]]$r_2_mu,tidydat) %>%
    #   mutate(specpar = modelpar1[[1]])



    distpar_all <- summarizedistpar(coll,tidydat,model)



    correlationsSDalpha <- summarizecorrelations(coll,coll[[2]]$L_1,"L_1")


    list(grouplevel_noncrit_sd = grouplevel_noncrit_sd,
         grouplevel_crit_sd = grouplevel_crit_sd,
         grouplevel_item_sd = grouplevel_item_sd,
         population_weights = population_weights,
         crit_weights = crit_weights,

         individual_distpar = distpar_all$individual_dist,
         population_distpar = distpar_all$population_dist,
         individual_crit = distpar_all$individual_crit,
         population_crit = distpar_all$population_crit,
         summary_correlations_ppt = correlationsSDalpha$summary_corr,
         raw_correlations_ppt = correlationsSDalpha$raw_corr,
         Exp = expstring)

  } else if (model %in% c("DPgamSDT","M0gamSDT","DPRMSDT","DPRM2SDT")){

    modelpar1 <- c("mu","disc","crit")

      modelpar <- c(rep(modelpar1[[1]],coef/2),rep(modelpar1[[2]],coef/2),rep(modelpar1[[3]],nthres))

    allSD <- summarizesds(coll)
    participantsd <-  allSD$participantsd %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar = factor(par,
                              levels = paste0("sd",c(1:length(modelpar))),
                              labels = modelpar))


    grouplevel_item_sd <- allSD$itemsd %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar = factor(par,labels=modelpar1[modelpar1!="crit"]))



    seqb <-  coef/length(modelpar1[modelpar1!="crit"])
    bound <- seq(1,coef,seqb)
    #bound <- ifelse(bound[2]==2 | length(bound) == 1,1,bound[2])

    if(max(bound)==2){
      grouplevel_noncrit_sd <- participantsd %>% filter(specpar !="crit") %>%
        mutate(par = as.character(par)) %>%
        mutate(condition_tmp = parse_number(par)) %>%
        mutate(specpar = case_when(condition_tmp == 1 ~ modelpar[[1]],
                                   condition_tmp == 2 ~ modelpar[[2]],
                                   TRUE ~ "NA")) %>%
        mutate(condition = condlabels) %>%
        select(-condition_tmp)
    } else {

      grouplevel_noncrit_sd <- participantsd %>% filter(specpar !="crit")  %>%
        mutate(par = as.character(par)) %>%
        mutate(condition_tmp = parse_number(par)) %>%
        mutate(specpar = case_when(condition_tmp < bound[[2]] ~ modelpar[[1]],
                                   condition_tmp >= bound[[2]] ~ modelpar[[2]],
                                   TRUE ~ "NA")) %>%
        mutate(condition = case_when(condition_tmp < bound[[2]] ~ condition_tmp,
                                     condition_tmp >= bound[[2]] ~ condition_tmp - length(condlabels) ,
                                     TRUE ~ 0)) %>%
        mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%

        select(-condition_tmp)
    }


    grouplevel_crit_sd <- participantsd %>% filter(specpar =="crit") %>%
      mutate(condition = parse_number(par)) %>%
      mutate(condition = factor(condition,labels = c(paste0("cr",c(1:nthres)))))

    population_weights <- bind_rows(summarizeweights(coll[[2]]$b, "b") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar1[[1]]),
                                    summarizeweights(coll[[2]]$b_disc, "b_disc") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar1[[2]]))

    crit_weights <- summarizeweights(coll[[2]]$b_gamma,"b_gamma") %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar =modelpar1[[3]])

    # itemeffects_mu <-summarizeitemeffects(coll,coll[[2]]$r_2_mu,tidydat) %>%
    #                            mutate(specpar = modelpar1[[1]])
    #


    distpar_all <- summarizedistpar(coll,tidydat,model)


    correlationsSDalpha <- summarizecorrelations(coll,coll[[2]]$L_1,"L_1")
    correlationsSDbeta <- summarizecorrelations(coll,coll[[2]]$L_2,"L_2")

    list(grouplevel_noncrit_sd = grouplevel_noncrit_sd,
         grouplevel_crit_sd = grouplevel_crit_sd,
         grouplevel_item_sd = grouplevel_item_sd,
         population_weights = population_weights,
         crit_weights = crit_weights,
         #individual_items_mu = itemeffects_mu,
         #individual_items_sigma = itemeffects_sigma,
         individual_distpar = distpar_all$individual_dist,
         population_distpar = distpar_all$population_dist,
         individual_crit = distpar_all$individual_crit,
         population_crit = distpar_all$population_crit,
         summary_correlations_ppt = correlationsSDalpha$summary_corr,
         raw_correlations_ppt = correlationsSDalpha$raw_corr,
         summary_correlations_item = correlationsSDbeta$summary_corr,
         raw_correlations_item = correlationsSDbeta$raw_corr,
         Exp = expstring)

  } else if (model %in% c("MASDT")){
    modelpar1 <- c("mu","disc","nonattf","crit")


    modelpar <- c(rep(modelpar1[[1]],coef/3),rep(modelpar1[[2]],coef/3),rep(modelpar1[[3]],coef/3),rep(modelpar1[[4]],nthres))

    allSD <- summarizesds(coll)
    participantsd <-  allSD$participantsd %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar = factor(par,
                              levels = paste0("sd",c(1:length(modelpar))),
                              labels = modelpar))


    grouplevel_item_sd <- allSD$itemsd %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar = factor(par,labels=modelpar1[modelpar1!="crit"]))



    seqb <-  coef/length(modelpar1[modelpar1!="crit"])
    bound <- seq(1,coef,seqb)
    #bound <- ifelse(bound[2]==2 | length(bound) == 1,1,bound[2])

    if(max(bound)==3){
      grouplevel_noncrit_sd <- participantsd %>% filter(specpar !="crit") %>%
        mutate(par = as.character(par)) %>%
        mutate(condition_tmp = parse_number(par)) %>%
        mutate(specpar = case_when(condition_tmp == 1 ~ modelpar[[1]],
                                   condition_tmp == 2 ~ modelpar[[2]],
                                   condition_tmp == 3 ~ modelpar[[3]],
                                   TRUE ~ "NA")) %>%
        mutate(condition = condlabels) %>%
        select(-condition_tmp)
    } else {

      grouplevel_noncrit_sd <- participantsd %>% filter(specpar !="crit")  %>%
        mutate(par = as.character(par)) %>%
        mutate(condition_tmp = parse_number(par)) %>%
        mutate(specpar = case_when(condition_tmp < bound[[2]] ~ modelpar[[1]],
                                   condition_tmp < bound[[3]] ~ modelpar[[2]],
                                   condition_tmp >= bound[[3]] ~ modelpar[[3]],
                                   TRUE ~ "NA")) %>%
        mutate(condition = case_when(condition_tmp < bound[[2]] ~ condition_tmp,
                                     condition_tmp < bound[[3]] ~ condition_tmp - length(condlabels),
                                     condition_tmp >= bound[[3]] ~  condition_tmp - (length(condlabels) *2) ,
                                     TRUE ~ 0)) %>%
        mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
        select(-condition_tmp)
    }


    grouplevel_crit_sd <- participantsd %>% filter(specpar =="crit") %>%
      mutate(condition = parse_number(par)) %>%
      mutate(condition = factor(condition,labels = c(paste0("cr",c(1:nthres)))))

    population_weights <- bind_rows(summarizeweights(coll[[2]]$b, "b") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar1[[1]]),
                                    summarizeweights(coll[[2]]$b_disc, "b_disc") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar1[[2]]),
                                    summarizeweights(coll[[2]]$b_nonattmu, "b_nonattmu") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar1[[3]]))

    crit_weights <- summarizeweights(coll[[2]]$b_gamma,"b_gamma") %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar =modelpar1[[4]])

    # itemeffects_mu <-summarizeitemeffects(coll,coll[[2]]$r_2_mu,tidydat) %>%
    #                            mutate(specpar = modelpar1[[1]])
    #


    distpar_all <- summarizedistpar(coll,tidydat,model)


    correlationsSDalpha <- summarizecorrelations(coll,coll[[2]]$L_1,"L_1")
    correlationsSDbeta <- summarizecorrelations(coll,coll[[2]]$L_2,"L_2")

    list(grouplevel_noncrit_sd = grouplevel_noncrit_sd,
         grouplevel_crit_sd = grouplevel_crit_sd,
         grouplevel_item_sd = grouplevel_item_sd,
         population_weights = population_weights,
         crit_weights = crit_weights,
         #individual_items_mu = itemeffects_mu,
         #individual_items_sigma = itemeffects_sigma,
         individual_distpar = distpar_all$individual_dist,
         population_distpar = distpar_all$population_dist,
         individual_crit = distpar_all$individual_crit,
         population_crit = distpar_all$population_crit,
         summary_correlations_ppt = correlationsSDalpha$summary_corr,
         raw_correlations_ppt = correlationsSDalpha$raw_corr,
         summary_correlations_item = correlationsSDbeta$summary_corr,
         raw_correlations_item = correlationsSDbeta$raw_corr,
         Exp = expstring)

  } else if (model %in% c("2HTM")){
    modelpar <- c("DO","DN","go")

    modelpar2 <- c(rep("DO",coll[[1]]$G_DO),rep("DN",coll[[1]]$G_DN),"go")

    seqb <-  (ncoef-1)/2
    bound <- seq(1,ncoef,seqb)
    #bound <- ifelse(bound[2]==2,1,bound[2])

    allSD <- summarizesds(coll)
    participantsd <-  allSD$participantsd %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar = factor(par,
                              levels = paste0("sd",c(1:length(modelpar2))),
                              labels = modelpar2))


    grouplevel_item_sd <- allSD$itemsd %>%
      mutate(par = as.character(par)) %>%
      mutate(specpar = factor(par,labels=modelpar))


      if(max(bound) == 3){
        grouplevel_noncrit_sd <-participantsd %>%
          mutate(par = as.character(par)) %>%
          mutate(condition_tmp = parse_number(par)) %>%
          mutate(specpar = case_when(condition_tmp == 1 ~ modelpar[[1]],
                                     condition_tmp == 2 ~ modelpar[[2]],
                                     condition_tmp == 3 ~ modelpar[[3]],
                                     TRUE ~ "NA")) %>%
          mutate(condition = condlabels) %>%
          select(-condition_tmp)
      } else {

        if(expstring=="GKH1999_e1"){
          tidydat <- tidydat %>% mutate(condition = paste0(condition_freq,"_",condition_studyduration))
          condlabels <- unique(tidydat$condition)
        }

        grouplevel_noncrit_sd <- participantsd %>%
          mutate(par = as.character(par)) %>%
          mutate(condition_tmp = parse_number(par)) %>%
          mutate(specpar = case_when(condition_tmp <= coll[[1]]$G_DO~ modelpar[[1]],
                                     condition_tmp <= bound[[3]] - 1 ~ modelpar[[2]],
                                     condition_tmp == bound[[3]] ~ modelpar[[3]],
                                     TRUE ~ "NA")) %>%
          mutate(condition = case_when(condition_tmp <= coll[[1]]$G_DO ~ condition_tmp,
                                       condition_tmp <= (bound[[3]] - 1) ~ condition_tmp - length(condlabels),
                                       condition_tmp == bound[[3]] ~ length(unique(tidydat$condition))+1,
                                       TRUE ~ 0)) %>%
          mutate(condition = factor(condition,labels = c(unique(tidydat$condition),"None"))) %>%
          select(-condition_tmp)
      }




    rmprobs <- bind_rows(summarizermprobs(coll,coll[[2]]$s_m,"s_m"),
                         summarizermprobs(coll,coll[[2]]$s_m,"a_m_o"),
                         summarizermprobs(coll,coll[[2]]$s_m,"a_m_n"))

    rmweights <- bind_rows(summarizeweights(coll[[2]]$mu_s,"mu_s"),
                           summarizeweights(coll[[2]]$mu_a_o,"mu_a_o"),
                           summarizeweights(coll[[2]]$mu_a_n,"mu_a_n"))


      population_weights <- bind_rows(summarizeweights(coll[[2]]$b_DO,"b_DO") %>%
                                        mutate(par = as.character(par)) %>%
                                        mutate(condition =parse_number(par)) %>%
                                        mutate(condition = factor(condition,labels = unique(tidydat%>% filter(isold==1) %>% .$condition))) %>%
                                        mutate(specpar = modelpar[[1]]),
                                      summarizeweights(coll[[2]]$b_DN,"b_DN") %>%
                                        mutate(par = as.character(par)) %>%
                                        mutate(condition =  parse_number(par)) %>%
                                        mutate(condition = factor(condition,labels = unique(tidydat %>% filter(isold==0) %>% .$condition))) %>%
                                        mutate(specpar = modelpar[[2]]),
                                      summarizeweights(coll[[2]]$b_go,"b_go") %>%

                                        mutate(condition =  "None") %>%

                                        mutate(specpar = modelpar[[3]])
      )



    distpar_all <- summarizedistpar(coll,tidydat,model)
    correlationsSDalpha <- summarizecorrelations(coll,coll[[2]]$L_1,"L_1")
    correlationsSDbeta <- summarizecorrelations(coll,coll[[2]]$L_2,"L_2")

    list(grouplevel_noncrit_sd = grouplevel_noncrit_sd,

         grouplevel_item_sd = grouplevel_item_sd,
         population_weights = population_weights,
         rm_weights = rmweights,
         #individual_items_mu = itemeffects_mu,
         #individual_items_sigma = itemeffects_sigma,
         individual_distpar = distpar_all$individual_dist,
         population_distpar = distpar_all$population_dist,
         rm_probs = rmprobs,

         summary_correlations_ppt = correlationsSDalpha$summary_corr,
         raw_correlations_ppt = correlationsSDalpha$raw_corr,
         summary_correlations_item = correlationsSDbeta$summary_corr,
         raw_correlations_item = correlationsSDbeta$raw_corr,
         Exp = expstring)
  }



}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
softmax <- function(par){
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk)))
  }
  val <- exp(par - Lk)
  return(val)
}

prepTheta <- function(tmpdat,fit, model){

  if (model %in% c("UVgamSDT")){

    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_gamma =  apply(fit$draws("b_gamma")[,,],3,rbind),

                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),

                       Intercept = apply(fit$draws("crits")[,,],3,rbind),

                       L_1 = apply(fit$draws("L_1")[,,],3,rbind),
                       #mu_cr = apply(fit$draws("mu_cr")[,,],3,rbind),
                       #sigma_cr = apply(fit$draws("sigma_cr")[,,],3,rbind),

                       sd_2 = apply(fit$draws("sd_2")[,,],3,rbind),
                       r_2_mu = apply(fit$draws("r_2_mu")[,,],3,rbind))

    # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1


    extractdat <- list(J_1 = tmpdat$J_1,
                       nthres = tmpdat$nthres,
                       Y = tmpdat$Y,
                       X = tmpdat$X,
                       X_disc = tmpdat$X_disc,
                       X_gamma = tmpdat$X_gamma,

                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       disc_Z_1 = tmpdat$Z_disc,
                       gamma_Z_1 = tmpdat$Z_gamma,
                       J_2 = tmpdat$J_2,
                       Z_2_mu = tmpdat$Z_2_mu)

  } else if (model %in% c("DPRM2SDT","DPRMSDT")){
    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("Intercept")[,,],3,rbind),
                       s_m = apply(fit$draws("s_m")[,,],3,rbind),
                       mu_s = apply(fit$draws("mu_s")[,,],3,rbind),
                       L_1 = apply(fit$draws("L_1")[,,],3,rbind),
                       mu_cr = apply(fit$draws("mu_cr")[,,],3,rbind),
                       sigma_cr = apply(fit$draws("sigma_cr")[,,],3,rbind))

    # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1


    extractdat <- list(J_1 = tmpdat$J_1,
                       nthres = tmpdat$nthres,
                       Y = tmpdat$Y,
                       X = tmpdat$X,
                       X_disc = tmpdat$X_disc,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       sigma_Z_1 = tmpdat$Z_disc,
                       M2=tmpdat$M_2,
                       J_2 = tmpdat$J_2)

  } else  if (model %in% c("DPgamSDT")){

    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_gamma = apply(fit$draws("b_gamma")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("crits")[,,],3,rbind),
                       L_1 = apply(fit$draws("L_1")[,,],3,rbind),
                       #mu_cr = apply(fit$draws("mu_cr")[,,],3,rbind),
                       #sigma_cr = apply(fit$draws("sigma_cr")[,,],3,rbind),
                       L_2 = apply(fit$draws("L_2")[,,],3,rbind),
                       sd_2 = apply(fit$draws("sd_2")[,,],3,rbind),
                       r_2_mu = apply(fit$draws("r_2_mu")[,,],3,rbind),
                       r_2_disc = apply(fit$draws("r_2_disc")[,,],3,rbind) )

    # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1


    extractdat <-  list(J_1 = tmpdat$J_1,
                        nthres = tmpdat$nthres,
                        Y = tmpdat$Y,
                        X = tmpdat$X,
                        X_disc = tmpdat$X_disc,
                        X_gamma = tmpdat$X_gamma,
                        M2=tmpdat$M_2,
                        ncoef=ncoef,
                        mu_Z_1 = tmpdat$Z,
                        disc_Z_1 = tmpdat$Z_disc,
                        gamma_Z_1 = tmpdat$Z_gamma,
                        J_2 = tmpdat$J_2,
                        Z_2_mu = tmpdat$Z_2_mu)

  } else  if (model %in% c("M0gamSDT")){

    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_gamma =  apply(fit$draws("b_gamma")[,,],3,rbind),

                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("crits")[,,],3,rbind),

                       L_1 = apply(fit$draws("L_1")[,,],3,rbind),
                       L_2 = apply(fit$draws("L_2")[,,],3,rbind),
                       sd_2 = apply(fit$draws("sd_2")[,,],3,rbind),
                       r_2_mu = apply(fit$draws("r_2_mu")[,,],3,rbind),
                       r_2_disc = apply(fit$draws("r_2_disc")[,,],3,rbind) )

    # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1


    extractdat <- list(J_1 = tmpdat$J_1,
                       nthres = tmpdat$nthres,
                       Y = tmpdat$Y,
                       X = tmpdat$X,
                       X_disc = tmpdat$X_disc,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       disc_Z_1 = tmpdat$Z_disc,
                       J_2 = tmpdat$J_2,
                       M2=tmpdat$M_2,
                       Z_2_mu = tmpdat$Z_2_mu,
                       Z_2_disc = tmpdat$Z_2_disc)

  } else if (model %in% c("MAgamSDT")){
    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_nonattmu = apply(fit$draws("b_nonattmu")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       b_gamma =  apply(fit$draws("b_gamma")[,,],3,rbind),
                       Intercept = apply(fit$draws("crits")[,,],3,rbind),
                       L_2 = apply(fit$draws("L_2")[,,],3,rbind),
                       sd_2 = apply(fit$draws("sd_2")[,,],3,rbind),
                       r_2_mu = apply(fit$draws("r_2_mu")[,,],3,rbind),
                       r_2_disc = apply(fit$draws("r_2_disc")[,,],3,rbind),
                       r_2_nonattmu = apply(fit$draws("r_2_nonattmu")[,,],3,rbind),
                       L_1 = apply(fit$draws("L_1")[,,],3,rbind))

    # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1


    extractdat <- list(J_1 = tmpdat$J_1,
                       nthres = tmpdat$nthres,
                       Y = tmpdat$Y,
                       X = tmpdat$X,
                       X_disc = tmpdat$X_disc,
                       X_nonattmu = tmpdat$X_nonattmu,
                       G = tmpdat$G,
                       G_disc = tmpdat$G_disc,
                       G_nonattmu = tmpdat$G_nonattmu,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       disc_Z_1 = tmpdat$Z_disc,
                       nonattemu_Z_1 = tmpdat$Z_nonattmu,
                       J_2 = tmpdat$J_2,
                       M2=tmpdat$M_2,
                       Z_2_mu = tmpdat$Z_2_mu,
                       Z_2_disc = tmpdat$Z_2_disc,

                       Z_2_nonattmu = tmpdat$Z_2_nonattmu)
  } else if (model %in% c("EVgamSDT","Gumbel")){

    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_gamma =  apply(fit$draws("b_gamma")[,,],3,rbind),

                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),

                       Intercept = apply(fit$draws("crits")[,,],3,rbind),

                       L_1 = apply(fit$draws("L_1")[,,],3,rbind),
                       #mu_cr = apply(fit$draws("mu_cr")[,,],3,rbind),
                       #sigma_cr = apply(fit$draws("sigma_cr")[,,],3,rbind),

                       sd_2 = apply(fit$draws("sd_2")[,,],3,rbind),
                       r_2_mu = apply(fit$draws("r_2_mu")[,,],3,rbind))

    # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1


    extractdat <- list(J_1 = tmpdat$J_1,
                       nthres = tmpdat$nthres,
                       Y = tmpdat$Y,
                       X = tmpdat$X,
                       X_gamma = tmpdat$X_gamma,
                       M2 = tmpdat$M_2,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       gamma_Z_1 = tmpdat$Z_gamma,
                       J_2 = tmpdat$J_2,
                       Z_2_mu = tmpdat$Z_2_mu)
  } else if (model %in% c("2HTM")){
    extractfit <- list(b_DO =  apply(fit$draws("b_DO")[,,],3,rbind),
                       b_DN =  apply(fit$draws("b_DN")[,,],3,rbind),
                       b_go = apply(fit$draws("b_go")[,,],3,rbind),

                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       mu_s = apply(fit$draws("mu_s")[,,],3,rbind),
                       mu_a_o = apply(fit$draws("mu_a_o")[,,],3,rbind),
                       mu_a_n = apply(fit$draws("mu_a_n")[,,],3,rbind),
                       s_m = apply(fit$draws("s_m")[,,],3,rbind),
                       a_m_o = apply(fit$draws("a_m_o")[,,],3,rbind),
                       a_m_n = apply(fit$draws("a_m_n")[,,],3,rbind),

                       L_1 = apply(fit$draws("L_1")[,,],3,rbind),
                       L_2 = apply(fit$draws("L_2")[,,],3,rbind),
                       sd_2 = apply(fit$draws("sd_2")[,,],3,rbind),
                       r_2_DO = apply(fit$draws("r_2_DO")[,,],3,rbind),

                       r_2_DN = apply(fit$draws("r_2_DN")[,,],3,rbind),

                       r_2_go = apply(fit$draws("r_2_go")[,,],3,rbind))

    # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1


    extractdat <- list(J_1 = tmpdat$J_1,
                       nthres = tmpdat$nthres,
                       Y = tmpdat$Y,
                       X_DO = tmpdat$X_DO,
                       X_DN = tmpdat$X_DN,
                       X_go = tmpdat$X_go,
                       G_DO = tmpdat$G_DO,
                       G_DN = tmpdat$G_DN,
                       G_go = tmpdat$G_go,
                       K_DO = tmpdat$K_DO,
                       K_DN = tmpdat$K_DN,
                       K_go = tmpdat$K_go,
                       ncoef=ncoef,
                       lenrm = tmpdat$Chalf,
                       M2=tmpdat$M_2,
                       J_2 = tmpdat$J_2,
                       Z_DO = tmpdat$Z_DO,
                       Z_DN = tmpdat$Z_DN,
                       Z_go = tmpdat$Z_go,
                       Z_2_DO = tmpdat$Z_2_DO,
                       Z_2_DN = tmpdat$Z_2_DN,
                       Z_2_go = tmpdat$Z_2_go)

  }

  list(extractdat,extractfit)

}
makeFreqDataItem <- function(exp,modelname){


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
    exp
  }
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
  tmpdat$Chalf <- floor(tmpdat$C / 2)
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
    if(model %in% c("UVSDT","DPSDT","M0SDT","DPRMSDT","DPRM2SDT")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc
    } else if(model %in% c("UVgamSDT","DPRMgamSDT","DPRM2gamSDT","DPgamSDT","M0gamSDT")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc + tmpdat$nthres
    }  else if(model %in% c("EVgamSDT","Gumbelgam")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$nthres
    } else if(model == "MAgamSDT"){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc + tmpdat$G_disc+ tmpdat$nthres
    } else {
      tmpdat$M_1 <- tmpdat$G

    }
  }

  if(model %in% c("UVSDT","Gumbel","EVSDT","UVgamSDT","EVgamSDT")){
    tmpdat$M_2 <- 1
  } else if(model %in% c("DPSDT","M0SDT","DPRMSDT","DPRM2SDT","DPgamSDT","M0gamSDT","DPRMgamSDT")){
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
summarizesds <- function(coll){

  # participant-effect, group-level distributions

  sds <- as_tibble(coll[[2]]$sd_1)

  sds <- bind_rows(as_tibble(coll[[2]]$sd_1) %>%

                     magrittr::set_colnames(paste0("sd",1:dim(sds)[2])) %>%
                     pivot_longer(cols=paste0("sd",1:dim(sds)[2]),
                                  names_to="par",values_to="value") %>%
                     mutate(par = factor(par,levels=paste0("sd",1:dim(sds)[2]))) %>%
                     group_by(par) %>%
                     mean_hdi(value,.width = c(.50,.8,.95))  %>%
                     pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                     mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                               labels=c("Mean","HDI","HDI"))) %>%
                     mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                     filter(quantdesc2!="Mean0.5",
                            quantdesc2!="Mean0.8"),
                   as_tibble(coll[[2]]$sd_1) %>%

                     magrittr::set_colnames(paste0("sd",1:dim(sds)[2])) %>%
                     pivot_longer(cols=paste0("sd",1:dim(sds)[2]),
                                  names_to="par",values_to="value") %>%
                     mutate(par = factor(par,levels=paste0("sd",1:dim(sds)[2]))) %>%
                     group_by(par) %>%
                     median_qi(value,.width = c(.50,.8,.95))  %>%
                     pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                     mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                               labels=c("Median","QI","QI"))) %>%
                     mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                     filter(quantdesc2!="Median0.5",
                            quantdesc2!="Median0.8"),

                   as_tibble(coll[[2]]$sd_1) %>%

                     magrittr::set_colnames(paste0("sd",1:dim(sds)[2])) %>%
                     pivot_longer(cols=paste0("sd",1:dim(sds)[2]),
                                  names_to="par",values_to="value") %>%
                     mutate(par = factor(par,levels=paste0("sd",1:dim(sds)[2]))) %>%
                     group_by(par) %>%
                     summarize(quants = sd(value)) %>%
                     mutate(.width=NA,
                            .point=NA,
                            .interval=NA,
                            quantdesc = "SD",
                            quantdesc2 = "SD"))

  sds2 <-  as_tibble(coll[[2]]$sd_2)

  sds2 <- bind_rows(as_tibble(coll[[2]]$sd_2) %>%

                     magrittr::set_colnames(paste0("sd",1:dim(sds2)[2])) %>%
                     pivot_longer(cols=paste0("sd",1:dim(sds2)[2]),
                                  names_to="par",values_to="value") %>%
                     mutate(par = factor(par,levels=paste0("sd",1:dim(sds2)[2]))) %>%
                     group_by(par) %>%
                     mean_hdi(value,.width = c(.50,.8,.95))  %>%
                     pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                     mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                               labels=c("Mean","HDI","HDI"))) %>%
                     mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                     filter(quantdesc2!="Mean0.5",
                            quantdesc2!="Mean0.8"),
                   as_tibble(coll[[2]]$sd_2) %>%

                     magrittr::set_colnames(paste0("sd",1:dim(sds2)[2])) %>%
                     pivot_longer(cols=paste0("sd",1:dim(sds2)[2]),
                                  names_to="par",values_to="value") %>%
                     mutate(par = factor(par,levels=paste0("sd",1:dim(sds2)[2]))) %>%
                     group_by(par) %>%
                     median_qi(value,.width = c(.50,.8,.95))  %>%
                     pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                     mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                               labels=c("Median","QI","QI"))) %>%
                     mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                     filter(quantdesc2!="Median0.5",
                            quantdesc2!="Median0.8"),

                   as_tibble(coll[[2]]$sd_2) %>%

                     magrittr::set_colnames(paste0("sd",1:dim(sds2)[2])) %>%
                     pivot_longer(cols=paste0("sd",1:dim(sds2)[2]),
                                  names_to="par",values_to="value") %>%
                     mutate(par = factor(par,levels=paste0("sd",1:dim(sds2)[2]))) %>%
                     group_by(par) %>%
                     summarize(quants = sd(value)) %>%
                     mutate(.width=NA,
                            .point=NA,
                            .interval=NA,
                            quantdesc = "SD",
                            quantdesc2 = "SD"))

  return(list(participantsd = sds,
              itemsd = sds2))
}
summarizecrits <- function(coll){

  # criteria, group-level distributions
nthres<-coll[[1]]$nthres
 tes<-bind_rows(bind_cols(as_tibble(coll[[2]]$mu_cr),
                     as_tibble(coll[[2]]$sigma_cr)
  ) %>%

    # mutate(it = i) %>%
    magrittr::set_colnames(c(paste0("mu_cr",c(1:nthres)),
                             paste0("sigma_cr",c(1:nthres)))) %>%
    pivot_longer(cols=c(paste0("mu_cr",c(1:nthres)),
                        paste0("sigma_cr",c(1:nthres))),
                 names_to="par",values_to="value") %>%
    group_by(par) %>%
    mean_hdi(value,.width = c(.50,0.8,.95))  %>%
    pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
    mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                              labels=c("Mean","HDI","HDI"))) %>%
    mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
    filter(quantdesc2!="Mean0.5") %>%
    filter(quantdesc2!="Mean0.8"),
  bind_cols(as_tibble(coll[[2]]$mu_cr),
            as_tibble(coll[[2]]$sigma_cr)
  ) %>%

    # mutate(it = i) %>%
    magrittr::set_colnames(c(paste0("mu_cr",c(1:nthres)),
                             paste0("sigma_cr",c(1:nthres)))) %>%
    pivot_longer(cols=c(paste0("mu_cr",c(1:nthres)),
                        paste0("sigma_cr",c(1:nthres))),
                 names_to="par",values_to="value") %>%
    group_by(par) %>%
    median_qi(value,.width = c(.50,0.8,.95))  %>%
    pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
    mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                              labels=c("Median","QI","QI"))) %>%
    mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
    filter(quantdesc2!="Median0.5") %>%
    filter(quantdesc2!="Median0.8"),
  bind_cols(as_tibble(coll[[2]]$mu_cr),
            as_tibble(coll[[2]]$sigma_cr)
  ) %>%

    # mutate(it = i) %>%
    magrittr::set_colnames(c(paste0("mu_cr",c(1:nthres)),
                             paste0("sigma_cr",c(1:nthres)))) %>%
    pivot_longer(cols=c(paste0("mu_cr",c(1:nthres)),
                        paste0("sigma_cr",c(1:nthres))),
                 names_to="par",values_to="value") %>%
    group_by(par) %>%
    summarize(quants = sd(value)) %>%
    mutate(.width=NA,
           .point=NA,
           .interval=NA,
           quantdesc = "SD",
           quantdesc2 = "SD"))


  pargr <- tes %>% distinct() %>%
    mutate(pargroup = case_when(grepl("mu_cr",par) ~ "mu_cr",
                                grepl("sigma_cr",par) ~ "sigma_cr")) %>%
    mutate(element = substring(par,nchar(par),nchar(par)))


  return(pargr)
}

summarizeweights <-function(speccoll,specstring){

  bs <- as_tibble(speccoll)

  bs <- bind_rows(as_tibble(speccoll) %>%

                    magrittr::set_colnames(paste0(specstring,1:dim(bs)[2])) %>%
                    pivot_longer(cols=paste0(specstring,1:dim(bs)[2]),
                                 names_to="par",values_to="value") %>%
                    group_by(par) %>%
                    mean_hdi(value,.width = c(.50,.8,.95))  %>%
                    pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                    mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                              labels=c("Mean","HDI","HDI"))) %>%
                    mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                    filter(quantdesc2!="Mean0.5",
                           quantdesc2!="Mean0.8"),
                  as_tibble(speccoll) %>%

                    magrittr::set_colnames(paste0(specstring,1:dim(bs)[2])) %>%
                    pivot_longer(cols=paste0(specstring,1:dim(bs)[2]),
                                 names_to="par",values_to="value") %>%
                    group_by(par) %>%
                    median_qi(value,.width = c(.50,.8,.95))  %>%
                    pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                    mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                              labels=c("Median","QI","QI"))) %>%
                    mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                    filter(quantdesc2!="Median0.5",
                           quantdesc2!="Median0.8"),

                  as_tibble(speccoll) %>%

                    magrittr::set_colnames(paste0(specstring,1:dim(bs)[2])) %>%
                    pivot_longer(cols=paste0(specstring,1:dim(bs)[2]),
                                 names_to="par",values_to="value") %>%
                    group_by(par) %>%
                    summarize(quants = sd(value)) %>%
                    mutate(.width=NA,
                           .point=NA,
                           .interval=NA,
                           quantdesc = "SD",
                           quantdesc2 = "SD"))

  return(bs)


}
summarizedistpar <- function(coll,tidydat,model){

  if (model %in% c("EVgamSDT","Gumbel")){

    res <- EVSDT_distpar(coll,tidydat)

  } else if (model == "UVgamSDT"){

    res <- UVSDT_distpar(coll,tidydat)

  } else if (model == "DPgamSDT"){

    res <- DPSDT_distpar(coll,tidydat)
  } else if (model=="M0gamSDT"){
    res <- M0SDT_distpar(coll,tidydat)
  } else if (model=="MAgamSDT"){
    res <- MASDT_distpar(coll,tidydat)
  } else if (model == "2HTM"){
   res <- H2TM_distpar(coll,tidydat)
  }

  res
}

summarizermprobs <- function(coll,speccoll,specstring){



  allSM <- NULL
  nthres<-coll[[1]]$nthres
  for(j in c(1:6000)){
    s_m_use <- matrix(speccoll[j,],
                      ncol = coll[[1]]$lenrm,byrow=F)

    sms <- as_tibble(s_m_use) %>% set_names(paste0(specstring,c(1:coll[[1]]$lenrm))) %>%
      mutate(J_1 = c(1:length(unique(coll[[1]]$J_1))),
             it = j)

    allSM <- allSM %>% bind_rows(sms)

  }



  r1shdi <-bind_rows(as_tibble(allSM) %>%
                       pivot_longer(cols=paste0(specstring,c(1:coll[[1]]$lenrm)),names_to="par",values_to = "value") %>%
                       group_by(J_1,par) %>%
                       mean_hdi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5",
                              quantdesc2!="Mean0.8"),

                     as_tibble(allSM) %>%
                       pivot_longer(cols=paste0(specstring,c(1:coll[[1]]$lenrm)),names_to="par",values_to = "value") %>%
                       group_by(J_1,par) %>%
                       median_qi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5",
                              quantdesc2!="Median0.8"),
                     as_tibble(allSM) %>%
                       pivot_longer(cols=paste0(specstring,c(1:coll[[1]]$lenrm)),names_to="par",values_to = "value") %>%
                       group_by(J_1,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par))

  r1shdi


}
summarizeitemeffects <- function(coll,speccoll,tidydat){

  r2_mu <-  apply(speccoll,1,function(x) x[coll[[1]]$J_2])

  lengthsam<-6000
  r2s <- NULL
  for (i in c(1:lengthsam)){

   samp <-   tibble(J_2 = coll[[1]]$J_2,
                    item = tidydat$item,
                    value=r2_mu[,i]) %>%
     distinct() %>%
     mutate(it = i)

   r2s <- r2s %>% bind_rows(samp)

  }


  r2shdi <-bind_rows(as_tibble(r2s) %>%

                       group_by(J_2,item) %>%
                       mean_hdi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5",
                              quantdesc2!="Mean0.8"),

                     as_tibble(r2s) %>%

                       group_by(J_2,item) %>%
                       median_qi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5",
                              quantdesc2!="Median0.8"),
                     as_tibble(r2s) %>%

                       group_by(J_2,item) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(item = factor(item))



}

EVSDT_distpar <- function(coll,tidydat){

  mu_o <- as.matrix(coll[[1]]$X)  %*% t(coll[[2]]$b)


  Y <- coll[[1]]$Y
  J_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(J_1))
  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  mutmpr1 <- asplit(tmpr_1[,c(1:((ncoef-nthres) * nJ_1))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)))

  gamtmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres) * nJ_1) + 1):dim(tmpr_1)[[2]])],1)
  gam_r_1 <- lapply(gamtmpr1,function(x) matrix(x,ncol=nthres))

  gamma <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma[[m]] <- tmp
  }


  lengthsam <- 6000

  r1s <- NULL
  grp <- NULL
  cr1s <- NULL
  crgrp <- NULL
  for (i in c(1:lengthsam)){


    mu = mu_o[,i]

    mur1 <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    gamma_use <- gamma[[i]] + gam_r_1[[i]]

    crit_tmp <- matrix(ncol=nthres,nrow=nJ_1)

    for(j in c(1:nJ_1)){
    crit_tmp[j,] <- 2*qnorm(head(cumsum(softmax(c(gamma_use[j,],0))),nthres))


    }

    crit_samp <- as_tibble(crit_tmp,.name_repair = "minimal") %>%
      setNames(., paste0("cr",c(1:nthres))) %>%
    mutate(J_1 = c(1:nJ_1)) %>%
      mutate(it = i)

    grcrit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma[[i]][1,],0))),nthres))

    crit_samp2 <- as_tibble(t(grcrit_tmp),.name_repair = "minimal") %>%
      setNames(., paste0("cr",c(1:nthres))) %>%
      mutate(it = i)



    samp <- tibble(J_1 = J_1,
                   mu_r1=mu + mur1,

                   isold = coll[[1]]$X[,1],
                   condition = tidydat$condition) %>%
      distinct() %>%
      #filter(isold==1) %>%
      mutate(it = i)

    samp2 <- tibble(J_1 = J_1,
                    mu_r1 = mu,

                    isold = coll[[1]]$X[,1],
                    condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      select(-J_1) %>%
      distinct() %>%
      mutate(it = i)

    r1s <- r1s %>% bind_rows(samp)
    grp <- grp %>% bind_rows(samp2)

    cr1s <- cr1s %>% bind_rows(crit_samp)
    crgrp <- crgrp %>% bind_rows(crit_samp2)
  }


  if(unique(tidydat$exp) == "GKH1999_e1"){

    r1s <- r1s %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
    grp <- grp %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
    }

  r1shdi <-bind_rows(as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%

                       group_by(J_1,condition,par) %>%
                       mean_hdi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5",
                              quantdesc2!="Mean0.8"),

                     as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%

                       group_by(J_1,condition,par) %>%
                       median_qi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5",
                              quantdesc2!="Median0.8"),
                     as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(J_1,condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par,levels=c("mu_r1"),
                        labels=c("mu")))


  critr1s <-bind_rows(as_tibble(cr1s) %>%
                       pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                       group_by(J_1,par) %>%

                       mean_hdi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5",
                              quantdesc2!="Mean0.8"),

                     as_tibble(cr1s) %>%
                       pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                       group_by(J_1,par) %>%

                       median_qi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5",
                              quantdesc2!="Median0.8"),
                     as_tibble(cr1s) %>%
                       pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                       group_by(J_1,par) %>%

                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par))


  grphdi <-bind_rows(as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       mean_hdi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5") %>%
                       filter(quantdesc2!="Mean0.8"),
                     as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       median_qi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5") %>%
                       filter(quantdesc2!="Median0.8"),
                     as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par,levels=c("mu_r1"),
                        labels=c("mu")))


  critgrp <-bind_rows(as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                       mean_hdi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5") %>%
                       filter(quantdesc2!="Mean0.8"),
                      as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                       median_qi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5") %>%
                       filter(quantdesc2!="Median0.8"),
                      as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD")) %>%
    mutate(par = factor(par))


  list(individual_dist = r1shdi,
       population_dist = grphdi,
       individual_crit = critr1s,
       population_crit = critgrp)
}
UVSDT_distpar <- function(coll,tidydat){

  # full <- bind_cols(as_tibble(as.matrix(coll[[1]]$X)) %>% set_names(c("m1","m2","m3")),
  #                   as_tibble(as.matrix(coll[[1]]$X_disc)) %>% set_names(c("d1","d2","d3")),
  #                   as_tibble(as.matrix(coll[[1]]$J_1)) %>% set_names(c("J_1")),
  #                   as_tibble(as.matrix(coll[[1]]$mu_Z_1)) %>% set_names(c("mz1")),
  #                   as_tibble(as.matrix(coll[[1]]$disc_Z_1)) %>% set_names(c("dz1")),
  #                   tidydat %>% ungroup() %>% select(condition)) %>%
  #   distinct()



   mu_o <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
   sigma_o <- as.matrix(coll[[1]]$X_disc) %*% t(coll[[2]]$b_disc)
  #mu_o <- as.matrix(full %>% select(m1,m2,m3)) %*% t(coll[[2]]$b)
  #sigma_o <- as.matrix(full %>% select(d1,d2,d3)) %*% t(coll[[2]]$b_disc)

  Y <- coll[[1]]$Y
  J_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(J_1))
  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  mutmpr1 <- asplit(tmpr_1[,c(1:((ncoef-nthres)/2 * nJ_1))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  sigmatmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/2 * nJ_1) + 1):((ncoef-nthres)/2 * nJ_1 * 2))],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  gamtmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/2 * nJ_1 * 2) + 1):dim(tmpr_1)[[2]])],1)
  gam_r_1 <- lapply(gamtmpr1,function(x) matrix(x,ncol=nthres))

  gamma <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma[[m]] <- tmp
  }

  lengthsam <- 6000

  r1s <- NULL
  grp <- NULL
  cr1s <- NULL
  crgrp <- NULL

  for (i in c(1:lengthsam)){


    mu = mu_o[,i]
    sigma = sigma_o[,i]
    # mur1 <- rowSums(apply(mu_r_1[[i]],2,function(x) x[full$J_1])  * as.matrix(full %>% select(mz1)))
    # sigmar1 <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[full$J_1])  * as.matrix(full %>% select(dz1)))
    #
    mur1 <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1 <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$disc_Z_1)



    samp <- tibble(J_1 = J_1,
                   mu_r1=mu + mur1,
                   sigma_r1 =1/exp(sigma + sigmar1),
                   isold =coll[[1]]$X[,1],
                   condition = tidydat$condition) %>%
      distinct() %>%
      #filter(isold==1) %>%
      mutate(it = i)

    samp2 <- tibble(J_1 = J_1,
                    mu_r1 = mu,
                    sigma_r1 =1/exp(sigma),
                    isold = coll[[1]]$X[,1],
                    condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      select(-J_1) %>%
      distinct() %>%
      mutate(it = i)

    r1s <- r1s %>% bind_rows(samp)
    grp <- grp %>% bind_rows(samp2)

    gamma_use <- gamma[[i]] + gam_r_1[[i]]

    crit_tmp <- matrix(ncol=nthres,nrow=nJ_1)

    for(j in c(1:nJ_1)){
      crit_tmp[j,] <- 2*qnorm(head(cumsum(softmax(c(gamma_use[j,],0))),nthres))


    }

    crit_samp <- as_tibble(crit_tmp,.name_repair = "minimal") %>%
      setNames(., paste0("cr",c(1:nthres))) %>%
      mutate(J_1 = c(1:nJ_1)) %>%
      mutate(it = i)

    grcrit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma[[i]][1,],0))),nthres))

    crit_samp2 <- as_tibble(t(grcrit_tmp),.name_repair = "minimal") %>%
      setNames(., paste0("cr",c(1:nthres))) %>%
      mutate(it = i)

    cr1s <- cr1s %>% bind_rows(crit_samp)
    crgrp <- crgrp %>% bind_rows(crit_samp2)


  }


  if(unique(tidydat$exp) == "GKH1999_e1"){

    r1s <- r1s %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
    grp <- grp %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
  }

  r1shdi <-bind_rows(as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%

                       group_by(J_1,condition,par) %>%
                       mean_hdi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5",
                              quantdesc2!="Mean0.8"),

                     as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%

                       group_by(J_1,condition,par) %>%
                       median_qi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5",
                              quantdesc2!="Median0.8"),
                     as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(J_1,condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par,levels=c("mu_r1","sigma_r1"),
                        labels=c("mu","sigma")))

  grphdi <-bind_rows(as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       mean_hdi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5") %>%
                       filter(quantdesc2!="Mean0.8"),
                     as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       median_qi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5") %>%
                       filter(quantdesc2!="Median0.8"),
                     as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par,levels=c("mu_r1","sigma_r1"),
                        labels=c("mu","sigma")))

  critr1s <-bind_rows(as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        mean_hdi(value,.width = c(.50,.80,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Mean","HDI","HDI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Mean0.5",
                               quantdesc2!="Mean0.8"),

                      as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        median_qi(value,.width = c(.50,.80,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Median","QI","QI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Median0.5",
                               quantdesc2!="Median0.8"),
                      as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        summarize(quants = sd(value)) %>%
                        mutate(.width=NA,
                               .point=NA,
                               .interval=NA,
                               quantdesc = "SD",
                               quantdesc2 = "SD"))%>%
    mutate(par = factor(par))

  critgrp <-bind_rows(as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        mean_hdi(value,.width = c(.50,.8,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Mean","HDI","HDI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Mean0.5") %>%
                        filter(quantdesc2!="Mean0.8"),
                      as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        median_qi(value,.width = c(.50,.8,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Median","QI","QI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Median0.5") %>%
                        filter(quantdesc2!="Median0.8"),
                      as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        summarize(quants = sd(value)) %>%
                        mutate(.width=NA,
                               .point=NA,
                               .interval=NA,
                               quantdesc = "SD",
                               quantdesc2 = "SD")) %>%
    mutate(par = factor(par))


  list(individual_crit = critr1s,
       population_crit = critgrp,
       individual_dist = r1shdi,
       population_dist = grphdi)
}
DPSDT_distpar <- function(coll,tidydat){

  mu_o <- as.matrix(coll[[1]]$X)  %*% t(coll[[2]]$b)
  sigma_o <- as.matrix(coll[[1]]$X_disc)  %*% t(coll[[2]]$b_disc)

  Y <- coll[[1]]$Y
  J_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(J_1))
  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  mutmpr1 <- asplit(tmpr_1[,c(1:((ncoef-nthres)/2 * nJ_1))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  sigmatmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/2 * nJ_1) + 1):((ncoef-nthres)/2 * nJ_1 * 2))],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))


  gamtmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/2 * nJ_1 * 2) + 1):dim(tmpr_1)[[2]])],1)
  gam_r_1 <- lapply(gamtmpr1,function(x) matrix(x,ncol=nthres))

  gamma <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma[[m]] <- tmp
  }


  lengthsam <- 6000

  r1s <- NULL
  grp <- NULL
  cr1s <- NULL
  crgrp <- NULL
  for (i in c(1:lengthsam)){


    mu = mu_o[,i]
    sigma = sigma_o[,i]
    mur1 <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1 <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$disc_Z_1)
    gamma_use <- gamma[[i]] + gam_r_1[[i]]

    crit_tmp <- matrix(ncol=nthres,nrow=nJ_1)

    for(j in c(1:nJ_1)){
      crit_tmp[j,] <- 2*qnorm(head(cumsum(softmax(c(gamma_use[j,],0))),nthres))


    }

    crit_samp <- as_tibble(crit_tmp,.name_repair = "minimal") %>%
      setNames(., paste0("cr",c(1:nthres))) %>%
      mutate(J_1 = c(1:nJ_1)) %>%
      mutate(it = i)

    grcrit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma[[i]][1,],0))),nthres))

    crit_samp2 <- as_tibble(t(grcrit_tmp),.name_repair = "minimal") %>%
      setNames(., paste0("cr",c(1:nthres))) %>%
      mutate(it = i)



    samp <- tibble(J_1 = J_1,
                   mu_r1=mu + mur1,
                   sigma_r1 =coll[[1]]$X[,1] * pnorm(sigma + sigmar1),
                   isold = coll[[1]]$X[,1],
                   condition = tidydat$condition) %>%
      distinct() %>%
      #filter(isold==1) %>%
      mutate(it = i)

    samp2 <- tibble(J_1 = J_1,
                    mu_r1 = mu,
                    sigma_r1 =coll[[1]]$X[,1] * pnorm(sigma),
                    isold = coll[[1]]$X[,1],
                    condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      select(-J_1) %>%
      distinct() %>%
      mutate(it = i)

    r1s <- r1s %>% bind_rows(samp)
    grp <- grp %>% bind_rows(samp2)

    cr1s <- cr1s %>% bind_rows(crit_samp)
    crgrp <- crgrp %>% bind_rows(crit_samp2)
  }


  if(unique(tidydat$exp) == "GKH1999_e1"){

    r1s <- r1s %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
    grp <- grp %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
  }

  r1shdi <-bind_rows(as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%

                       group_by(J_1,condition,par) %>%
                       mean_hdi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5",
                              quantdesc2!="Mean0.8"),

                     as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%

                       group_by(J_1,condition,par) %>%
                       median_qi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5",
                              quantdesc2!="Median0.8"),
                     as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(J_1,condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par,levels=c("mu_r1","sigma_r1"),
                        labels=c("mu","sigma")))


  critr1s <-bind_rows(as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        mean_hdi(value,.width = c(.50,.80,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Mean","HDI","HDI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Mean0.5",
                               quantdesc2!="Mean0.8"),

                      as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        median_qi(value,.width = c(.50,.80,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Median","QI","QI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Median0.5",
                               quantdesc2!="Median0.8"),
                      as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        summarize(quants = sd(value)) %>%
                        mutate(.width=NA,
                               .point=NA,
                               .interval=NA,
                               quantdesc = "SD",
                               quantdesc2 = "SD"))%>%
    mutate(par = factor(par))


  grphdi <-bind_rows(as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       mean_hdi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5") %>%
                       filter(quantdesc2!="Mean0.8"),
                     as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       median_qi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5") %>%
                       filter(quantdesc2!="Median0.8"),
                     as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par,levels=c("mu_r1","sigma_r1"),
                        labels=c("mu","sigma")))


  critgrp <-bind_rows(as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        mean_hdi(value,.width = c(.50,.8,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Mean","HDI","HDI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Mean0.5") %>%
                        filter(quantdesc2!="Mean0.8"),
                      as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        median_qi(value,.width = c(.50,.8,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Median","QI","QI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Median0.5") %>%
                        filter(quantdesc2!="Median0.8"),
                      as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        summarize(quants = sd(value)) %>%
                        mutate(.width=NA,
                               .point=NA,
                               .interval=NA,
                               quantdesc = "SD",
                               quantdesc2 = "SD")) %>%
    mutate(par = factor(par))


  list(individual_dist = r1shdi,
       population_dist = grphdi,
       individual_crit = critr1s,
       population_crit = critgrp)
}
M0SDT_distpar <- function(coll,tidydat){


  mu_o <- as.matrix(coll[[1]]$X)  %*% t(coll[[2]]$b)
  sigma_o <- as.matrix(coll[[1]]$X_disc)  %*% t(coll[[2]]$b_disc)

  Y <- coll[[1]]$Y
  J_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(J_1))
  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  mutmpr1 <- asplit(tmpr_1[,c(1:((ncoef-nthres)/2 * nJ_1))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  sigmatmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/2 * nJ_1) + 1):((ncoef-nthres)/2 * nJ_1 * 2))],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))


  gamtmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/2 * nJ_1 * 2) + 1):dim(tmpr_1)[[2]])],1)
  gam_r_1 <- lapply(gamtmpr1,function(x) matrix(x,ncol=nthres))

  gamma <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma[[m]] <- tmp
  }


  lengthsam <- 6000

  r1s <- NULL
  grp <- NULL
  cr1s <- NULL
  crgrp <- NULL
  for (i in c(1:lengthsam)){


    mu = mu_o[,i]
    sigma = sigma_o[,i]
    mur1 <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1 <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$disc_Z_1)
    gamma_use <- gamma[[i]] + gam_r_1[[i]]

    crit_tmp <- matrix(ncol=nthres,nrow=nJ_1)

    for(j in c(1:nJ_1)){
      crit_tmp[j,] <- 2*qnorm(head(cumsum(softmax(c(gamma_use[j,],0))),nthres))


    }

    crit_samp <- as_tibble(crit_tmp,.name_repair = "minimal") %>%
      setNames(., paste0("cr",c(1:nthres))) %>%
      mutate(J_1 = c(1:nJ_1)) %>%
      mutate(it = i)

    grcrit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma[[i]][1,],0))),nthres))

    crit_samp2 <- as_tibble(t(grcrit_tmp),.name_repair = "minimal") %>%
      setNames(., paste0("cr",c(1:nthres))) %>%
      mutate(it = i)



    samp <- tibble(J_1 = J_1,
                   mu_r1=mu + mur1,
                   sigma_r1 =coll[[1]]$X[,1] * pnorm(sigma + sigmar1),
                   isold = coll[[1]]$X[,1],
                   condition = tidydat$condition) %>%
      distinct() %>%
      #filter(isold==1) %>%
      mutate(it = i)

    samp2 <- tibble(J_1 = J_1,
                    mu_r1 = mu,
                    sigma_r1 =coll[[1]]$X[,1] * pnorm(sigma),
                    isold = coll[[1]]$X[,1],
                    condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      select(-J_1) %>%
      distinct() %>%
      mutate(it = i)

    r1s <- r1s %>% bind_rows(samp)
    grp <- grp %>% bind_rows(samp2)

    cr1s <- cr1s %>% bind_rows(crit_samp)
    crgrp <- crgrp %>% bind_rows(crit_samp2)
  }


  if(unique(tidydat$exp) == "GKH1999_e1"){

    r1s <- r1s %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
    grp <- grp %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
  }

  r1shdi <-bind_rows(as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%

                       group_by(J_1,condition,par) %>%
                       mean_hdi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5",
                              quantdesc2!="Mean0.8"),

                     as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%

                       group_by(J_1,condition,par) %>%
                       median_qi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5",
                              quantdesc2!="Median0.8"),
                     as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(J_1,condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par,levels=c("mu_r1","sigma_r1"),
                        labels=c("mu","sigma")))


  critr1s <-bind_rows(as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        mean_hdi(value,.width = c(.50,.80,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Mean","HDI","HDI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Mean0.5",
                               quantdesc2!="Mean0.8"),

                      as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        median_qi(value,.width = c(.50,.80,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Median","QI","QI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Median0.5",
                               quantdesc2!="Median0.8"),
                      as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        summarize(quants = sd(value)) %>%
                        mutate(.width=NA,
                               .point=NA,
                               .interval=NA,
                               quantdesc = "SD",
                               quantdesc2 = "SD"))%>%
    mutate(par = factor(par))


  grphdi <-bind_rows(as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       mean_hdi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5") %>%
                       filter(quantdesc2!="Mean0.8"),
                     as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       median_qi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5") %>%
                       filter(quantdesc2!="Median0.8"),
                     as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par,levels=c("mu_r1","sigma_r1"),
                        labels=c("mu","sigma")))


  critgrp <-bind_rows(as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        mean_hdi(value,.width = c(.50,.8,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Mean","HDI","HDI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Mean0.5") %>%
                        filter(quantdesc2!="Mean0.8"),
                      as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        median_qi(value,.width = c(.50,.8,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Median","QI","QI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Median0.5") %>%
                        filter(quantdesc2!="Median0.8"),
                      as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        summarize(quants = sd(value)) %>%
                        mutate(.width=NA,
                               .point=NA,
                               .interval=NA,
                               quantdesc = "SD",
                               quantdesc2 = "SD")) %>%
    mutate(par = factor(par))


  list(individual_dist = r1shdi,
       population_dist = grphdi,
       individual_crit = critr1s,
       population_crit = critgrp)
}
MASDT_distpar <- function(coll,tidydat){


  mu_o <- as.matrix(coll[[1]]$X)  %*% t(coll[[2]]$b)
  sigma_o <- as.matrix(coll[[1]]$X_disc) %*% t(coll[[2]]$b_disc)
  nonattmu_o <- as.matrix(coll[[1]]$X_nonattmu) %*% t(coll[[2]]$b_nonattmu)

  Y <- coll[[1]]$Y
  J_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(J_1))
  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  mutmpr1 <- asplit(tmpr_1[,c(1:((ncoef-nthres)/3 * nJ_1))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/3))

  sigmatmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/3 * nJ_1) + 1):((ncoef-nthres)/3 * nJ_1 * 2))],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/3))

  nonattmutmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/3 * nJ_1 * 2) + 1):((ncoef-nthres)/3 * nJ_1 * 3))],1)
  nonattmu_r_1 <- lapply(nonattmutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/3))

  gamtmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/3 * nJ_1 * 3) + 1):dim(tmpr_1)[[2]])],1)
  gam_r_1 <- lapply(gamtmpr1,function(x) matrix(x,ncol=nthres))

  gamma <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma[[m]] <- tmp
  }



  lengthsam <- 6000
  mur1 <- matrix(ncol = lengthsam,nrow=length(Y))
  sigmar1 <- matrix(ncol = lengthsam,nrow=length(Y))
  nonattmur1 <- matrix(ncol = lengthsam,nrow=length(Y))

  r1s <- NULL
  grp <- NULL
  cr1s <- NULL
  crgrp <- NULL
  for (i in c(1:lengthsam)){

    mu = mu_o[,i]
    sigo = sigma_o[,i]
    nonatt = nonattmu_o[,i]
    mur1 <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1 <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$disc_Z_1)
    nonattmur1 <- rowSums(apply(nonattmu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$nonattemu_Z_1)
    gamma_use <- gamma[[i]] + gam_r_1[[i]]

    crit_tmp <- matrix(ncol=nthres,nrow=nJ_1)

    for(j in c(1:nJ_1)){
      crit_tmp[j,] <- 2*qnorm(head(cumsum(softmax(c(gamma_use[j,],0))),nthres))


    }

    crit_samp <- as_tibble(crit_tmp,.name_repair = "minimal") %>%
      setNames(., paste0("cr",c(1:nthres))) %>%
      mutate(J_1 = c(1:nJ_1)) %>%
      mutate(it = i)

    grcrit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma[[i]][1,],0))),nthres))

    crit_samp2 <- as_tibble(t(grcrit_tmp),.name_repair = "minimal") %>%
      setNames(., paste0("cr",c(1:nthres))) %>%
      mutate(it = i)

    samp <- tibble(J_1 = J_1,
                   mu_r1=mu + mur1,
                   sigma_r1=coll[[1]]$X[,1] * pnorm(sigo + sigmar1),
                   nonattmutmp = coll[[1]]$X[,1] * pnorm(nonatt + nonattmur1),
                   isold = coll[[1]]$X[,1],
                   condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      mutate(it = i) %>%
      mutate(nonattmu_r1 = mu_r1 * nonattmutmp) %>%
      select(-nonattmutmp)

    samp2 <- tibble(J_1 = J_1,
                    mu_r1 = mu,
                    sigma_r1 = coll[[1]]$X[,1] * pnorm(sigo),
                    nonattmutmp = coll[[1]]$X[,1] * pnorm(nonatt),
                    isold = coll[[1]]$X[,1],
                    condition = tidydat$condition) %>%
      mutate(nonattmu_r1 = mu_r1 * nonattmutmp) %>%
      #filter(isold==1) %>%
      select(-J_1) %>%
      distinct() %>%
      mutate(it = i) %>%
      select(-nonattmutmp)

    r1s <- r1s %>% bind_rows(samp)
    grp <- grp %>% bind_rows(samp2)
    cr1s <- cr1s %>% bind_rows(crit_samp)
    crgrp <- crgrp %>% bind_rows(crit_samp2)
  }


  if(unique(tidydat$exp) == "GKH1999_e1"){

    r1s <- r1s %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
    grp <- grp %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
  }

  r1shdi <-bind_rows(as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1","nonattmu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%

                       group_by(J_1,condition,par) %>%
                       mean_hdi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5",
                              quantdesc2!="Mean0.8"),

                     as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1","nonattmu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%

                       group_by(J_1,condition,par) %>%
                       median_qi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5",
                              quantdesc2!="Median0.8"),
                     as_tibble(r1s) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1","nonattmu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(J_1,condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par,levels=c("mu_r1","sigma_r1","nonattmu_r1"),
                        labels=c("mu","sigma","nonattmu")))

  critr1s <-bind_rows(as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        mean_hdi(value,.width = c(.50,.80,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Mean","HDI","HDI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Mean0.5",
                               quantdesc2!="Mean0.8"),

                      as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        median_qi(value,.width = c(.50,.80,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Median","QI","QI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Median0.5",
                               quantdesc2!="Median0.8"),
                      as_tibble(cr1s) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(J_1,par) %>%

                        summarize(quants = sd(value)) %>%
                        mutate(.width=NA,
                               .point=NA,
                               .interval=NA,
                               quantdesc = "SD",
                               quantdesc2 = "SD"))%>%
    mutate(par = factor(par))


  grphdi <-bind_rows(as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1","nonattmu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       mean_hdi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5") %>%
                       filter(quantdesc2!="Mean0.8"),
                     as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1","nonattmu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       median_qi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5") %>%
                       filter(quantdesc2!="Median0.8"),
                     as_tibble(grp) %>%
                       pivot_longer(cols=c("mu_r1","sigma_r1","nonattmu_r1"),names_to="par",values_to = "value") %>%
                       mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par,levels=c("mu_r1","sigma_r1","nonattmu_r1"),
                        labels=c("mu","sigma","nonattmu")))


  critgrp <-bind_rows(as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        mean_hdi(value,.width = c(.50,.8,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Mean","HDI","HDI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Mean0.5") %>%
                        filter(quantdesc2!="Mean0.8"),
                      as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        median_qi(value,.width = c(.50,.8,.95)) %>%
                        pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                        mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                  labels=c("Median","QI","QI"))) %>%
                        mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                        filter(quantdesc2!="Median0.5") %>%
                        filter(quantdesc2!="Median0.8"),
                      as_tibble(crgrp) %>%
                        pivot_longer(cols=paste0("cr",c(1:nthres)),names_to="par",values_to = "value") %>%
                        group_by(par) %>%
                        summarize(quants = sd(value)) %>%
                        mutate(.width=NA,
                               .point=NA,
                               .interval=NA,
                               quantdesc = "SD",
                               quantdesc2 = "SD")) %>%
    mutate(par = factor(par))


  list(individual_dist = r1shdi,
       population_dist = grphdi,
       individual_crit = critr1s,
       population_crit = critgrp)
}
H2TM_distpar <- function(coll,tidydat){

  Y <- coll[[1]]$Y

  DO <- as.matrix(coll[[1]]$X_DO) %*% t(coll[[2]]$b_DO)
  DN <- as.matrix(coll[[1]]$X_DN) %*%  t(coll[[2]]$b_DN)
  go <- as.matrix(coll[[1]]$X_go) %*%  t(coll[[2]]$b_go)

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1

  lengthDO <- length(unique(J_1)) * coll[[1]]$G_DO
  lengthDN <- length(unique(J_1)) * coll[[1]]$G_DN
  lengthgo <- length(unique(J_1)) * coll[[1]]$G_go

  DOtmpr1 <- asplit(tmpr_1[,c(1:lengthDO)],1)
  DO_r_1 <- lapply(DOtmpr1,function(x) matrix(x,ncol=coll[[1]]$G_DO))

  DNtmpr1 <- asplit(tmpr_1[,c((lengthDO + 1) : (lengthDO + lengthDN))],1)
  DN_r_1 <- lapply(DNtmpr1,function(x) matrix(x,ncol=coll[[1]]$G_DN))

  gotmpr1 <- asplit(tmpr_1[,c((lengthDO + lengthDN + 1) : (lengthDO + lengthDN + lengthgo))],1)
  go_r_1 <- lapply(gotmpr1,function(x) matrix(x,ncol=coll[[1]]$G_go))


  lengthsam <- 6000

  r1s <- NULL
  grp <- NULL
  for (i in c(1:lengthsam)){


    Do = DO[,i]
    Dn = DN[,i]
    Go = go[,i]
    DOr1 <- rowSums(apply(DO_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_DO)
    DNr1 <- rowSums(apply(DN_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_DN)
    gor1 <- rowSums(apply(go_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_go)

    samp <- tibble(J_1 = J_1,
                   DO_r1=pnorm(Do + DOr1) * coll[[1]]$X_DO[,1],
                   DN_r1=pnorm(Dn + DNr1) * coll[[1]]$X_DN[,1],
                   go_r1=pnorm(Go + gor1) * coll[[1]]$X_go[,1],

                   isold = coll[[1]]$X_DO[,1],
                   condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      distinct() %>%
      mutate(it = i)

    samp2 <- tibble(J_1 = J_1,
                    DO_r1=pnorm(Do) * coll[[1]]$X_DO[,1],
                    DN_r1=pnorm(Dn) * coll[[1]]$X_DN[,1],
                    go_r1=pnorm(Go) * coll[[1]]$X_go[,1],
                    isold = coll[[1]]$X_DO[,1],
                    condition = tidydat$condition) %>%
      select(-J_1) %>%
      distinct() %>%
      mutate(it = i)

    r1s <- r1s %>% bind_rows(samp)
    grp <- grp %>% bind_rows(samp2)
  }

  # if(unique(tidydat$exp) == "GKH1999_e1"){
  #
  #   r1s <- r1s %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
  #   grp <- grp %>% mutate(condition = ifelse(condition == "LowF_New_0" & isold == 0, "HighF_New_0",condition))
  # }

  DOr1s <- r1s %>% select(DO_r1,J_1,condition,it,isold) %>%
    filter(isold==1) %>%
    mutate(par = "DO_r1") %>%
    rename("value" = "DO_r1")
  DNr1s <- r1s %>% select(DN_r1,J_1,condition,it,isold) %>%
    filter(isold==0) %>%
    mutate(par = "DN_r1") %>%
    rename("value" = "DN_r1")
  gor1s <- r1s %>% select(go_r1,J_1,condition,it) %>%
    mutate(condition = "None") %>%
    distinct() %>%
    mutate(par = "go_r1") %>%
    rename("value" = "go_r1")



  r1shdi <-bind_rows(as_tibble(bind_rows(DOr1s,DNr1s,gor1s)) %>%
                       #pivot_longer(cols=c("DO_r1","DN_r1","go_r1"),names_to="par",values_to = "value") %>%
                       #mutate(condition = ifelse(!grepl("_0|_1",condition),paste0(condition,"_",isold),condition)) %>%
                       group_by(J_1,condition,par) %>%
                       mean_hdi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5",
                              quantdesc2!="Mean0.8"),

                     as_tibble(bind_rows(DOr1s,DNr1s,gor1s)) %>%
                       group_by(J_1,condition,par) %>%
                       median_qi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5",
                              quantdesc2!="Median0.8"),
                     as_tibble(bind_rows(DOr1s,DNr1s,gor1s)) %>%
                       group_by(J_1,condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par))

  DOgrp <- grp %>% select(DO_r1,condition,it,isold) %>%
    filter(isold==1) %>%
    mutate(par = "DO_r1") %>%
    rename("value" = "DO_r1")
  DNgrp <- grp %>% select(DN_r1,condition,it,isold) %>%
    filter(isold==0) %>%
    mutate(par = "DN_r1") %>%
    rename("value" = "DN_r1")
  gogrp <- grp %>% select(go_r1,condition,it) %>%
    mutate(condition = "None") %>%
    distinct() %>%
    mutate(par = "go_r1") %>%
    rename("value" = "go_r1")

  grphdi <-bind_rows(as_tibble(bind_rows(DOgrp,DNgrp,gogrp)) %>%
                       group_by(condition,par) %>%
                       mean_hdi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5") %>%
                       filter(quantdesc2!="Mean0.8"),
                     as_tibble(bind_rows(DOgrp,DNgrp,gogrp)) %>%
                       group_by(condition,par) %>%
                       median_qi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5") %>%
                       filter(quantdesc2!="Median0.8"),
                     as_tibble(bind_rows(DOgrp,DNgrp,gogrp)) %>%
                       group_by(condition,par) %>%
                       summarize(quants = sd(value)) %>%
                       mutate(.width=NA,
                              .point=NA,
                              .interval=NA,
                              quantdesc = "SD",
                              quantdesc2 = "SD"))%>%
    mutate(par = factor(par))

  list(individual_distpar = r1shdi,
       population_distpar = grphdi)

}

summarizecorrelations <- function(coll,speccoll,specstring){

  if(specstring == "L_1"){
    ncoef = coll[[1]]$ncoef
  } else {
    ncoef = coll[[1]]$M2
  }
  if(dim(speccoll)[2] != 1){

    totL1 <- apply(speccoll,1,function(x) extendcor2(x,ncoef))
    corrs <- as_tibble(t(totL1)) %>% set_names(paste0("c",c(1:dim(totL1)[1]))) %>%
      pivot_longer(cols=paste0("c",c(1:dim(totL1)[1])),names_to="par",
                   values_to = "value")

    corrcomp <- bind_rows(corrs %>%
                            group_by(par) %>%
                            mean_hdi(value,.width = c(.50,.8,.95)) %>%
                            pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                            mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                      labels=c("Mean","HDI","HDI"))) %>%
                            mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                            filter(quantdesc2!="Mean0.5") %>%
                            filter(quantdesc2!="Mean0.8"),
                          corrs %>%
                            group_by(par) %>%
                            median_qi(value,.width = c(.50,.8,.95)) %>%
                            pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                            mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                      labels=c("Median","QI","QI"))) %>%
                            mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                            filter(quantdesc2!="Median0.5") %>%
                            filter(quantdesc2!="Median0.8"),
                          corrs %>%
                            group_by(par) %>%
                            summarize(quants = sd(value)) %>%
                            mutate(.width=NA,
                                   .point=NA,
                                   .interval=NA,
                                   quantdesc = "SD",
                                   quantdesc2 = "SD"))

    list(summary_corr = corrcomp,
         raw_corr = corrs)

  } else {
    list(summary_corr = NULL,
         raw_corr = NULL)
  }

}

extendcor2 <- function(singleL_1,ncoef){

  matrix(singleL_1,nrow=ncoef,byrow=F) %*% t(matrix(singleL_1,nrow=ncoef,byrow=F))

  #zeromat <- matrix(0,nthres,ncoef)

  #cbind(rbind(cor_coef,zeromat),rbind(t(zeromat),diag(nthres)))

}


#6,58,60,18,65,77,70,34,64

for (i in exps[c(4,10,12,14,15,33,41,63,68,79)]){

   model<-"2HTM"
   modeln<-"H2TM"
  fitfiles <- list.files(paste0(pathfit),full.names=T)


  if(length(fitfiles[grepl(i,fitfiles)])!=0){

    fitfiles <- list.files(pathfit,full.names=T)
    fitfiles <- fitfiles[grepl(model,fitfiles)]
    fit <- readRDS(fitfiles[grepl(i,fitfiles)])

    exp_pred <- d1pred[[which(explist_pred==i)]] %>%
      mutate(rating = as.numeric(as.character(rating)))

    tidydat <- makeFreqDataItem(exp_pred,model)

    tmp <- makeTmpDatItem(tidydat,model)
    coll <- prepTheta(tmp,fit,model)

   # distpar_all <- summarizedistpar(coll,tidydat,model)
  res <- summarizeparameters(i,model,coll)
  saveRDS(res,file=paste0("D:/HS01/ParameterSummary/ItemEffectsRandomCrit/",modeln,"/parametersummary_sdind_full_",model,"_",i,".rds"))

  itemeffects_mu <-summarizeitemeffects(coll,coll[[2]]$r_2_DO,tidydat) %>%
    mutate(specpar = "DO")
  saveRDS(itemeffects_mu,file=paste0("D:/HS01/ParameterSummary/ItemEffectsRandomCrit/",modeln,"/parametersummary_itemDO_full_",model,"_",i,".rds"))

  itemeffects_sigma <-  summarizeitemeffects(coll,coll[[2]]$r_2_DN,tidydat) %>%
    mutate(specpar = "DN")
  saveRDS(itemeffects_sigma,file=paste0("D:/HS01/ParameterSummary/ItemEffectsRandomCrit/",modeln,"/parametersummary_itemDN_full_",model,"_",i,".rds"))

  itemeffects_nonattmu <-  summarizeitemeffects(coll,coll[[2]]$r_2_go,tidydat) %>%
    mutate(specpar = "go")
  saveRDS(itemeffects_nonattmu,file=paste0("D:/HS01/ParameterSummary/ItemEffectsRandomCrit/",modeln,"/parametersummary_itemgo_full_",model,"_",i,".rds"))


  }
  }



