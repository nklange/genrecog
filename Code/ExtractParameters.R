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

pathfit <- "D:/HS01/Fits/NoItemEffects/Full/"
pathpar <- "D:/HS01/ParameterSummary/NoItemEffects/"




summarizeparameters <- function(expstring,model){


if(model == "2HTM"){
  modeln <- "H2TM"
} else {
  modeln <- model
}
  fitfiles <- list.files(paste0(pathfit,modeln,"/"),full.names=T)
  fit <- readRDS(fitfiles[grepl(expstring,fitfiles)])

  exp_pred <- d1pred[[which(explist_pred==expstring)]] %>%
    mutate(rating = as.numeric(as.character(rating))) %>%
    mutate(itemInd="Ind")

  tidydat <- makeFreqDataInd(exp_pred,model)

  tmp <- makeTmpDatInd(tidydat,model)
  coll <- prepTheta_Ind(tmp,fit,model)
  ncoef <- coll[[1]]$ncoef

  if(unique(tidydat$manipulation) == "Between"){
    condlabels<-"None"
  } else {
    condlabels <-unique(tidydat$condition)

  }

  if(model %in% c("EVSDT","Gumbel")){

    modelpar <- c("mu")

    grouplevel_sd <- summarizesds(coll) %>%
      mutate(par = as.character(par)) %>%
      mutate(condition = parse_number(par)) %>%
      mutate(condition = factor(condition,labels = condlabels)) %>%
      mutate(specpar = modelpar[[1]])

   grouplevel_crits <- summarizecrits(coll)


   population_weights <- summarizeweights(coll[[2]]$b,"b")%>%
     mutate(par = as.character(par)) %>%
     mutate(condition = parse_number(par)) %>%
     mutate(condition = as.numeric(condition)) %>%
     mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
     mutate(specpar = modelpar[[1]])

   distpar_all <- summarizedistpar(coll,tidydat,model)
   individualcrits <- summarizeindcrits(coll)

   correlations <- summarizecorrelations(coll)

   list(grouplevel_sd = grouplevel_sd,
        population_weights = population_weights,
        individual_distpar = distpar_all$individual_dist,
        population_distpar = distpar_all$population_dist,
        grouplevel_crits = grouplevel_crits,
        individual_crit = individualcrits,
        summary_correlations = correlations$summary_corr,
        raw_correlations = correlations$raw_corr,
        Exp = expstring)

  } else if (model %in% c("UVSDT","DPSDT","M0SDT")){

    modelpar <- c("mu","disc")

    seqb <-  ncoef/length(modelpar)
    bound <- seq(1,ncoef,seqb)
    #bound <- ifelse(bound[2]==2 | length(bound) == 1,1,bound[2])

    if(max(bound)==2){
      grouplevel_sd <- summarizesds(coll) %>%
        mutate(par = as.character(par)) %>%
        mutate(condition_tmp = parse_number(par)) %>%
        mutate(specpar = case_when(condition_tmp == 1 ~ modelpar[[1]],
                                   condition_tmp == 2 ~ modelpar[[2]],
                                   TRUE ~ "NA")) %>%
        mutate(condition = condlabels) %>%
        select(-condition_tmp)
    } else {

      grouplevel_sd <- summarizesds(coll) %>%
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

    grouplevel_crits <- summarizecrits(coll)

    population_weights <- bind_rows(summarizeweights(coll[[2]]$b,"b") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar[[1]]),
                                    summarizeweights(coll[[2]]$b_disc,"b_disc") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =  parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar[[2]]))

    distpar_all <- summarizedistpar(coll,tidydat,model)
    individualcrits <- summarizeindcrits(coll)
    correlations <- summarizecorrelations(coll)

    list(grouplevel_sd = grouplevel_sd,
         population_weights = population_weights,
         individual_distpar = distpar_all$individual_dist,
         population_distpar = distpar_all$population_dist,
         grouplevel_crits = grouplevel_crits,
         individual_crit = individualcrits,
         summary_correlations = correlations$summary_corr,
         raw_correlations = correlations$raw_corr,
         Exp = expstring)


  }else if (model %in% c("DPRMSDT","DPRM2SDT")){

    modelpar <- c("mu","disc")

    seqb <-  ncoef/length(modelpar)
    bound <- seq(1,ncoef,seqb)
    #bound <- ifelse(bound[2]==2 | length(bound) == 1,1,bound[2])

    if(max(bound)==2){
      grouplevel_sd <- summarizesds(coll) %>%
        mutate(par = as.character(par)) %>%
        mutate(condition_tmp = parse_number(par)) %>%
        mutate(specpar = case_when(condition_tmp == 1 ~ modelpar[[1]],
                                   condition_tmp == 2 ~ modelpar[[2]],
                                   TRUE ~ "NA")) %>%
        mutate(condition = condlabels) %>%
        select(-condition_tmp)
    } else {

      grouplevel_sd <- summarizesds(coll) %>%
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

    grouplevel_crits <- summarizecrits(coll)

    population_weights <- bind_rows(summarizeweights(coll[[2]]$b,"b") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar[[1]]),
                                    summarizeweights(coll[[2]]$b_disc,"b_disc") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =  parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar[[2]]))

    distpar_all <- summarizedistpar(coll,tidydat,model)
    individualcrits <- summarizeindcrits(coll)
    correlations <- summarizecorrelations(coll)

    rmprobs <- summarizermprobs(coll,coll[[2]]$s_m,"s_m")

    rmweights <- summarizerweights(coll[[2]]$mu_s,"mu_s")

    list(grouplevel_sd = grouplevel_sd,
         grouplevel_crits = grouplevel_crits,
         population_weights = population_weights,
         individual_distpar = distpar_all$individual_dist,
         population_distpar = distpar_all$population_dist,
         summary_correlations = correlations$summary_corr,
         raw_correlations = correlations$raw_corr,
         individual_crit = individualcrits,
         rm_weights = rmweights,
         rm_probs = rmprobs,
         Exp = expstring)

  } else if (model %in% c("MASDT")){

    modelpar <- c("mu","disc","nonattf")

    seqb <-  ncoef/length(modelpar)
    bound <- seq(1,ncoef,seqb)
    #bound <- ifelse(bound[2]==2,1,bound[2])

    if(max(bound) == 3){
      grouplevel_sd <- summarizesds(coll) %>%
        mutate(par = as.character(par)) %>%
        mutate(condition_tmp = parse_number(par)) %>%
        mutate(specpar = case_when(condition_tmp == 1 ~ modelpar[[1]],
                                   condition_tmp == 2 ~ modelpar[[2]],
                                   condition_tmp == 3 ~ modelpar[[3]],
                                   TRUE ~ "NA")) %>%
        mutate(condition = condlabels) %>%
        select(-condition_tmp)
    } else {

      grouplevel_sd <- summarizesds(coll) %>%
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

    grouplevel_crits <- summarizecrits(coll)

    population_weights <- bind_rows(summarizeweights(coll[[2]]$b,"b") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar[[1]]),
                                    summarizeweights(coll[[2]]$b_disc,"b_disc") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =  parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar[[2]]),
                                    summarizeweights(coll[[2]]$b_nonattmu,"b_nonattmu") %>%
                                      mutate(par = as.character(par)) %>%
                                      mutate(condition =  parse_number(par)) %>%
                                      mutate(condition = factor(condition,labels = unique(tidydat$condition))) %>%
                                      mutate(specpar = modelpar[[3]])
                                    )

    distpar_all <- summarizedistpar(coll,tidydat,model)
    individualcrits <- summarizeindcrits(coll)
    correlations <- summarizecorrelations(coll)

    list(grouplevel_sd = grouplevel_sd,
         grouplevel_crits = grouplevel_crits,
         population_weights = population_weights,
         individual_distpar = distpar_all$individual_dist,
         population_distpar = distpar_all$population_dist,
         individual_crit = individualcrits,
         summary_correlations = correlations$summary_corr,
         raw_correlations = correlations$raw_corr,
         Exp = expstring)

  } else if (model %in% c("2HTM")){
    modelpar <- c("DO","DN","go")

    modelpar2 <- c(rep("DO",coll[[1]]$G_DO),rep("DN",coll[[1]]$G_DN),"go")

    seqb <-  (ncoef-1)/2
    bound <- seq(1,ncoef,seqb)
    #bound <- ifelse(bound[2]==2,1,bound[2])

    if(!is.null(coll[[2]]$b_neut)){
      modelparn<-c("DO","DN","go","neut")
      seqb <-  (ncoef-2)/2
      bound <- seq(1,ncoef,seqb)

      if(max(bound) == 4){
        grouplevel_sd <- summarizesds(coll) %>%
          mutate(par = as.character(par)) %>%
          mutate(condition_tmp = parse_number(par)) %>%
          mutate(specpar = case_when(condition_tmp == 1 ~ modelparn[[1]],
                                     condition_tmp == 2 ~ modelparn[[2]],
                                     condition_tmp == 3 ~ modelparn[[3]],
                                     condition_tmp == 4 ~ modelparn[[4]],
                                     TRUE ~ "NA")) %>%
          mutate(condition = condlabels) %>%
          select(-condition_tmp)
      } else {

        grouplevel_sd <- summarizesds(coll) %>%
          mutate(par = as.character(par)) %>%
          mutate(condition_tmp = parse_number(par)) %>%
          mutate(specpar = case_when(condition_tmp <= coll[[1]]$G_DO ~ modelparn[[1]],
                                     condition_tmp <= (bound[[3]]-1) ~ modelparn[[2]],
                                     condition_tmp == (bound[[3]]) ~ modelparn[[3]],
                                     condition_tmp == bound[[3]]+1 ~ modelparn[[4]],
                                     TRUE ~ "NA")) %>%
          mutate(condition = case_when(condition_tmp <= coll[[1]]$G_DO ~ condition_tmp,
                                       condition_tmp <= (bound[[3]]-1) ~ condition_tmp - length(condlabels),
                                       condition_tmp == (bound[[3]]) ~ length(unique(tidydat$condition))+1,
                                       condition_tmp == bound[[3]]+1 ~ length(unique(tidydat$condition))+1,
                                       TRUE ~ 0)) %>%
          mutate(condition = factor(condition,labels = c(unique(tidydat$condition),"None"))) %>%
          select(-condition_tmp)
      }

    } else {
    if(max(bound) == 3){
      grouplevel_sd <- summarizesds(coll) %>%
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

      grouplevel_sd <- summarizesds(coll) %>%
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

    }


    rmprobs <- bind_rows(summarizermprobs(coll,coll[[2]]$s_m,"s_m"),
                         summarizermprobs(coll,coll[[2]]$s_m,"a_m_o"),
                         summarizermprobs(coll,coll[[2]]$s_m,"a_m_n"))

    rmweights <- bind_rows(summarizeweights(coll[[2]]$mu_s,"mu_s"),
                           summarizeweights(coll[[2]]$mu_a_o,"mu_a_o"),
                           summarizeweights(coll[[2]]$mu_a_n,"mu_a_n"))

    if(!is.null(coll[[2]]$b_neut)){

      population_weights <- bind_rows(summarizeweights(coll[[2]]$b_DO,"b_DO") %>%
                                        mutate(par = as.character(par)) %>%
                                        mutate(condition =parse_number(par)) %>%
                                        mutate(condition = factor(condition,labels = unique(tidydat %>% filter(isold==1) %>% .$condition))) %>%
                                        mutate(specpar = modelpar[[1]]),
                                      summarizeweights(coll[[2]]$b_DN,"b_DN") %>%
                                        mutate(par = as.character(par)) %>%
                                        mutate(condition =  parse_number(par)) %>%
                                        mutate(condition = factor(condition,labels = unique(tidydat %>% filter(isold==0) %>% .$condition))) %>%
                                        mutate(specpar = modelpar[[2]]),
                                      summarizeweights(coll[[2]]$b_go,"b_go") %>%

                                        mutate(condition =  "None") %>%

                                        mutate(specpar = modelpar[[3]]),
                                      summarizeweights(coll[[2]]$b_neut,"b_neut") %>%

                                        mutate(condition =  "None") %>%

                                        mutate(specpar = "neut")
      )
    } else{
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
    }


    distpar_all <- summarizedistpar(coll,tidydat,model)
    correlations <- summarizecorrelations(coll)

    list(grouplevel_sd = grouplevel_sd,
         population_weights = population_weights,
         rm_weights = rmweights,
         rm_probs = rmprobs,
         individual_distpar = distpar_all$individual_dist,
         population_distpar = distpar_all$population_dist,
         summary_correlations = correlations$summary_corr,
         raw_correlations = correlations$raw_corr,
         Exp = expstring)

  }


}

prepTheta_Ind <- function(tmpdat,fit, model){

  if (model %in% c("UVgamSDT","DPgamSDT")){
    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_gamma =  apply(fit$draws("b_gamma")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("crits")[,,],3,rbind),
                       L_1 = apply(fit$draws("L_1")[,,],3,rbind))

    # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1 - tmpdat$nthres

    extractdat <- list(J_1 = tmpdat$J_1,
                       nthres = tmpdat$nthres,
                       Y = tmpdat$Y,
                       X = tmpdat$X,
                       X_disc = tmpdat$X_disc,
                       X_gamma = tmpdat$X_gamma,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       sigma_Z_1 = tmpdat$Z_disc,
                       gamma_Z_1 = tmpdat$Z_gamma,
                       J_2 = tmpdat$J_2)

  } else if (model %in% c("UVSDT","DPSDT")){
    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("Intercept")[,,],3,rbind),
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
                       J_2 = tmpdat$J_2)

  } else if (model %in% c("DPRMSDT","DPRM2SDT")){
    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("Intercept")[,,],3,rbind),
                       mu_s = apply(fit$draws("mu_s")[,,],3,rbind),
                       s_m = apply(fit$draws("s_m")[,,],3,rbind),
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
                       X_nonattmu = tmpdat$X_nonattmu,
                       G = tmpdat$G,
                       G_disc = tmpdat$G_disc,
                       G_nonattmu = tmpdat$G_nonattmu,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       sigma_Z_1 = tmpdat$Z_disc,
                       nonattemu_Z_1 = tmpdat$Z_nonattmu,
                       lenrm = tmpdat$Chalf,
                       J_2 = tmpdat$J_2)

  } else if (model %in% c("EVSDT","Gumbel")){

    #tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    #Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1


    extractdat <- list(J_1 = tmpdat$J_1,
                       nthres = tmpdat$nthres,
                       Y = tmpdat$Y,
                       X = tmpdat$X,
                       # X_disc = tmpdat$X_disc[test,],
                       ncoef=ncoef,
                       mu_Z_1 = as.matrix(tmpdat$Z))#as.matrix(Z_1[,c(1:(ncoef/2))]),
    #sigma_Z_1 = as.matrix(Z_1[,c(((ncoef/2)+1):ncoef)]))

    if(ncoef == 1) {
      r_1 = apply(fit$draws("r_1_mu")[,,],3,rbind)
    } else {
      r_1 = apply(fit$draws("r_1_mu")[,,],3,rbind)
    }

    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       # b_disc =  b_disc <- apply(fit$draws("b_disc")[,,],3,rbind),

                       r_1 = r_1,
                       Intercept = apply(fit$draws("Intercept")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       L_1 = apply(fit$draws("L_1")[,,],3,rbind),
                       mu_cr = apply(fit$draws("mu_cr")[,,],3,rbind),
                       sigma_cr = apply(fit$draws("sigma_cr")[,,],3,rbind))

  } else if (model %in% c("M0SDT")){
    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       #b_nonattmu = apply(fit$draws("b_nonattmu")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("Intercept")[,,],3,rbind),

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
                       #X_nonattmu = tmpdat$X_nonattmu,
                       G = tmpdat$G,
                       G_disc = tmpdat$G_disc,
                      # G_nonattmu = tmpdat$G_nonattmu,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       sigma_Z_1 = tmpdat$Z_disc)
                      # nonattemu_Z_1 = tmpdat$Z_nonattmu)
  }else if (model %in% c("MASDT")){
    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_nonattmu = apply(fit$draws("b_nonattmu")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("Intercept")[,,],3,rbind),

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
                       X_nonattmu = tmpdat$X_nonattmu,
                       G = tmpdat$G,
                       G_disc = tmpdat$G_disc,
                       G_nonattmu = tmpdat$G_nonattmu,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       sigma_Z_1 = tmpdat$Z_disc,
                       nonattemu_Z_1 = tmpdat$Z_nonattmu,
                       J_2 = tmpdat$J_2)
  } else if (model %in% c("2HTM")){

    if (tmpdat$exp %in% c("HUW2015_e1","LBA2019_e1","LBA2019_e2","LM2020_e1")) {
      extractfit <- list(b_DO =  apply(fit$draws("b_DO")[,,],3,rbind),
                         b_DN =  apply(fit$draws("b_DN")[,,],3,rbind),
                         b_go = apply(fit$draws("b_go")[,,],3,rbind),
                         b_neut = apply(fit$draws("b_neut")[,,],3,rbind),
                         sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                         r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                         mu_s = apply(fit$draws("mu_s")[,,],3,rbind),
                         mu_a_o = apply(fit$draws("mu_a_o")[,,],3,rbind),
                         mu_a_n = apply(fit$draws("mu_a_n")[,,],3,rbind),
                         s_m = apply(fit$draws("s_m")[,,],3,rbind),
                         a_m_o = apply(fit$draws("a_m_o")[,,],3,rbind),
                         a_m_n = apply(fit$draws("a_m_n")[,,],3,rbind),

                         L_1 = apply(fit$draws("L_1")[,,],3,rbind))

      # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
      # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
      ncoef = tmpdat$M_1


      extractdat <- list(J_1 = tmpdat$J_1,
                         nthres = tmpdat$nthres,
                         Y = tmpdat$Y,
                         X_DO = tmpdat$X_DO,
                         X_DN = tmpdat$X_DN,
                         X_go = tmpdat$X_go,
                         X_neut = tmpdat$X_neut,
                         G_DO = tmpdat$G_DO,
                         G_DN = tmpdat$G_DN,
                         G_go = tmpdat$G_go,
                         G_neut = tmpdat$G_neut,
                         K_DO = tmpdat$K_DO,
                         K_DN = tmpdat$K_DN,
                         K_go = tmpdat$K_go,
                         K_neut = tmpdat$K_neut,
                         ncoef=ncoef,
                         lenrm = tmpdat$Chalf,
                         Z_DO = tmpdat$Z_DO,
                         Z_DN = tmpdat$Z_DN,
                         Z_go = tmpdat$Z_go,
                         Z_neut = tmpdat$Z_neut)
    } else {
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

                         L_1 = apply(fit$draws("L_1")[,,],3,rbind))

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
                         Z_DO = tmpdat$Z_DO,
                         Z_DN = tmpdat$Z_DN,
                         Z_go = tmpdat$Z_go)
    }

  }

  list(extractdat,extractfit)

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
makeTmpDatInd <- function(expdats,model){

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

  tmpdat$Y <- as.matrix(expdats[,7:dim(expdats)[2]])


  tmpdat$C <- dim(expdats)[2] - 6
  if(model=="DPRM2SDT"){
    tmpdat$Chalf <- 2
    } else {
      tmpdat$Chalf <- floor(tmpdat$C / 2)
    }
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


    tmpdat$K <- length(unique(expdats$condition))


    tmpdat$K_disc <- length(unique(expdats$condition))


    tmpdat$K_nonattmu <- length(unique(expdats$condition))
    tmpdat$X <- as.matrix(selectorth)
    tmpdat$X_disc <- as.matrix(selectorth)
    tmpdat$X_nonattmu <- as.matrix(selectorth) # this is MASDT specific
    tmpdat$X_gamma <- matrix(data = 1, nrow =  tmpdat$N , ncol = tmpdat$nthres)
    tmpdat$Z_gamma <- matrix(data = 1, nrow =  tmpdat$N , ncol = tmpdat$nthres)

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
    if(model %in% c("UVSDT","DPSDT","M0SDT","DPRMSDT","DPRM2SDT")){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc
    } else if(model == "UVgamSDT"){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc + tmpdat$nthres
    } else if(model == "MASDT"){
      tmpdat$M_1 <- tmpdat$G + tmpdat$G_disc + tmpdat$G_disc
    } else {
      tmpdat$M_1 <- tmpdat$G

    }
  }

  tmpdat$N_1 <- length(unique(expdats$id))

  tmpdat$J_1 <- as.numeric(factor(expdats$id))

  tmpdat$prior_only <- 0
  tmpdat$prior_alpha_mu <- c(10,3,1,0.1)[1:tmpdat$Chalf]
  tmpdat$prior_alpha_scale <- c(2,1,0.5,0.25)[1:tmpdat$Chalf]
  tmpdat$prior_alpha_mu_a <- c(1,2,3,5)[1:tmpdat$Chalf]
  tmpdat$prior_alpha_scale_a <- c(0.25,.5,1,1.5)[1:tmpdat$Chalf]
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

  return(sds)
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



  return(list(individual_crits = pargr))
}
summarizeindcrits <- function(coll){

  # criteria, group-level distributions
  nthres<-coll[[1]]$nthres

  allcrit <- NULL

  for(j in c(1:6000)){

  crit <- as_tibble(matrix(coll[[2]]$Intercept[j,],
         ncol = nthres,byrow=F)) %>%
    set_names(paste0("cr",c(1:nthres))) %>%
    mutate(J_1 = unique(coll[[1]]$J_1)) %>%
    mutate(it = j)

    allcrit <- allcrit %>% bind_rows(crit)
  }



  tes<-bind_rows(

   allcrit %>%
    pivot_longer(cols=c(paste0("cr",c(1:nthres))),
                 names_to="par",values_to="value") %>%
    group_by(J_1,par) %>%
    mean_hdi(value,.width = c(.50,0.8,.95))  %>%
    pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
    mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                              labels=c("Mean","HDI","HDI"))) %>%
    mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
    filter(quantdesc2!="Mean0.5") %>%
    filter(quantdesc2!="Mean0.8"),
   allcrit %>%
     pivot_longer(cols=c(paste0("cr",c(1:nthres))),
                  names_to="par",values_to="value") %>%
     group_by(J_1,par) %>%
    median_qi(value,.width = c(.50,0.8,.95))  %>%
    pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
    mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                              labels=c("Median","QI","QI"))) %>%
    mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
    filter(quantdesc2!="Median0.5") %>%
    filter(quantdesc2!="Median0.8"),
   allcrit %>%
     pivot_longer(cols=c(paste0("cr",c(1:nthres))),
                  names_to="par",values_to="value") %>%
     group_by(J_1,par) %>%
    summarize(quants = sd(value)) %>%
    mutate(.width=NA,
           .point=NA,
           .interval=NA,
           quantdesc = "SD",
           quantdesc2 = "SD"))


  pargr <- tes %>% distinct()




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
summarizedistpar <- function(coll,tidydat,model){

  if (model %in% c("EVSDT","Gumbel")){

    res <- EVSDT_distpar(coll,tidydat)
  } else if (model == "UVSDT"){

    res <- UVSDT_distpar(coll,tidydat)

  } else if (model == "DPSDT"){

    res <- DPSDT_distpar(coll,tidydat)
  } else if (model=="M0SDT"){
    res <- M0SDT_distpar(coll,tidydat)
  } else if (model=="MASDT"){
    res <- MASDT_distpar(coll,tidydat)
  } else if (model %in% c("DPRMSDT","DPRM2SDT")){
    res <- DPSDT_distpar(coll,tidydat)
  } else if (model %in% c("2HTM")){
    if(!is.null(coll[[2]]$b_neut)){
      res <- H2TM_distpar_neut(coll,tidydat)
    } else {
    res <- H2TM_distpar_std(coll,tidydat)
    }
  }

  res
}
EVSDT_distpar <- function(coll,tidydat){

  mu_o <- as.matrix(coll[[1]]$X)  %*% t(coll[[2]]$b)


  Y <- coll[[1]]$Y
  J_1 <- coll[[1]]$J_1
  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  mu_r_1 <- lapply(asplit(tmpr_1,1),function(x) matrix(x,ncol=ncoef))



  lengthsam <- 6000
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])

  r1s <- NULL
  grp <- NULL
  for (i in c(1:lengthsam)){

    mu = mu_o[,i]

    mur1 <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)


    samp <- tibble(J_1 = J_1,
                   mu_r1=mu + mur1,

                   isold = coll[[1]]$X[,1],
                   condition = tidydat$condition) %>%
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




  list(individual_dist = r1shdi,
       population_dist = grphdi)
}
UVSDT_distpar <- function(coll,tidydat){


  mu_o <- as.matrix(coll[[1]]$X)  %*% t(coll[[2]]$b)
  sigma_o <- as.matrix(coll[[1]]$X_disc) %*% t(coll[[2]]$b_disc)

  Y <- coll[[1]]$Y
  J_1 <- coll[[1]]$J_1
  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  mutmpr1 <- asplit(tmpr_1[,c(1:(dim(tmpr_1)[2]/2))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=ncoef/2))

  sigmatmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/2) + 1):dim(tmpr_1)[2])],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=ncoef/2))

  lengthsam <- 6000
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])

  r1s <- NULL
  grp <- NULL
  for (i in c(1:lengthsam)){

    mu = mu_o[,i]
    sigo = sigma_o[,i]
    mur1 <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1 <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)

    samp <- tibble(J_1 = J_1,
                   mu_r1=mu + mur1,
                   sigma_r1=1/exp(sigo + sigmar1),
                   isold = coll[[1]]$X[,1],
                   condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      mutate(it = i)

    samp2 <- tibble(J_1 = J_1,
                    mu_r1 = mu,
                    sigma_r1 = 1/exp(sigo),
                    isold = coll[[1]]$X[,1],
                    condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      select(-J_1) %>%
      distinct() %>%
      mutate(it = i)

    r1s <- r1s %>% bind_rows(samp)
    grp <- grp %>% bind_rows(samp2)
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


  list(individual_dist = r1shdi,
       population_dist = grphdi)
}
DPSDT_distpar <- function(coll,tidydat){



  mu_o <- as.matrix(coll[[1]]$X)  %*% t(coll[[2]]$b)
  sigma_o <- as.matrix(coll[[1]]$X_disc) %*% t(coll[[2]]$b_disc)

  Y <- coll[[1]]$Y
  J_1 <- coll[[1]]$J_1
  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  mutmpr1 <- asplit(tmpr_1[,c(1:(dim(tmpr_1)[2]/2))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=ncoef/2))

  sigmatmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/2) + 1):dim(tmpr_1)[2])],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=ncoef/2))

  lengthsam <- 6000
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])

  r1s <- NULL
  grp <- NULL
  for (i in c(1:lengthsam)){

    mu = mu_o[,i]
    sigo = sigma_o[,i]
    mur1 <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1 <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)

    samp <- tibble(J_1 = J_1,
                   mu_r1=mu + mur1,
                   sigma_r1=coll[[1]]$X[,1] *pnorm(sigo + sigmar1),
                   isold = coll[[1]]$X[,1],
                   condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      mutate(it = i)

    samp2 <- tibble(J_1 = J_1,
                    mu_r1 = mu,
                    sigma_r1 = coll[[1]]$X[,1] * pnorm(sigo),
                    isold = coll[[1]]$X[,1],
                    condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      select(-J_1) %>%
      distinct() %>%
      mutate(it = i)

    r1s <- r1s %>% bind_rows(samp)
    grp <- grp %>% bind_rows(samp2)
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


  list(individual_dist = r1shdi,
       population_dist = grphdi)
}
M0SDT_distpar <- function(coll,tidydat){


  mu_o <- as.matrix(coll[[1]]$X)  %*% t(coll[[2]]$b)
  sigma_o <- as.matrix(coll[[1]]$X_disc) %*% t(coll[[2]]$b_disc)

  Y <- coll[[1]]$Y
  J_1 <- coll[[1]]$J_1
  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  mutmpr1 <- asplit(tmpr_1[,c(1:(dim(tmpr_1)[2]/2))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=ncoef/2))

  sigmatmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/2) + 1):dim(tmpr_1)[2])],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=ncoef/2))

  lengthsam <- 6000
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])

  r1s <- NULL
  grp <- NULL
  for (i in c(1:lengthsam)){

    mu = mu_o[,i]
    sigo = sigma_o[,i]
    mur1 <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1 <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)

    samp <- tibble(J_1 = J_1,
                   mu_r1=mu + mur1,
                   sigma_r1=coll[[1]]$X[,1] * pnorm(sigo + sigmar1),
                   isold = coll[[1]]$X[,1],
                   condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      mutate(it = i)

    samp2 <- tibble(J_1 = J_1,
                    mu_r1 = mu,
                    sigma_r1 = coll[[1]]$X[,1] * pnorm(sigo),
                    isold = coll[[1]]$X[,1],
                    condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      select(-J_1) %>%
      distinct() %>%
      mutate(it = i)

    r1s <- r1s %>% bind_rows(samp)
    grp <- grp %>% bind_rows(samp2)
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


  list(individual_dist = r1shdi,
       population_dist = grphdi)
}
MASDT_distpar <- function(coll,tidydat){


  mu_o <- as.matrix(coll[[1]]$X)  %*% t(coll[[2]]$b)
  sigma_o <- as.matrix(coll[[1]]$X_disc) %*% t(coll[[2]]$b_disc)
  nonattmu_o <- as.matrix(coll[[1]]$X_nonattmu) %*% t(coll[[2]]$b_nonattmu)

  Y <- coll[[1]]$Y
  J_1 <- coll[[1]]$J_1
  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  mutmpr1 <- asplit(tmpr_1[,c(1:(dim(tmpr_1)[2]/3))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=ncoef/3))

  sigmatmpr1 <- asplit(tmpr_1[,c((dim(tmpr_1)[2]/3) + 1):((dim(tmpr_1)[2]/3) *2)],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=ncoef/3))

  nonattmutmpr1 <- asplit(tmpr_1[,c((((dim(tmpr_1)[2]/3) *2) + 1):dim(tmpr_1)[2])],1)
  nonattmu_r_1 <- lapply(nonattmutmpr1,function(x) matrix(x,ncol=ncoef/3))

  lengthsam <- 6000
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  nonattmur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])

  r1s <- NULL
  grp <- NULL
  for (i in c(1:lengthsam)){

    mu = mu_o[,i]
    sigo = sigma_o[,i]
    nonatt = nonattmu_o[,i]
    mur1 <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1 <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$sigma_Z_1)
    nonattmur1 <- rowSums(apply(nonattmu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$nonattemu_Z_1)

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
                        labels=c("mu","sigma_r1","mu_unattended")))


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
                        labels=c("mu","sigma_r1","mu_unattended")))


  list(individual_dist = r1shdi,
       population_dist = grphdi)
}
H2TM_distpar_std <- function(coll,tidydat){

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
H2TM_distpar_neut <- function(coll,tidydat){

  Y <- coll[[1]]$Y

  DO <- as.matrix(coll[[1]]$X_DO) %*% t(coll[[2]]$b_DO)
  DN <- as.matrix(coll[[1]]$X_DN) %*%  t(coll[[2]]$b_DN)
  go <- as.matrix(coll[[1]]$X_go) %*%  t(coll[[2]]$b_go)
  neut <- as.matrix(coll[[1]]$X_neut)%*%  t(coll[[2]]$b_neut)

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1

  lengthDO <- length(unique(J_1)) * coll[[1]]$G_DO
  lengthDN <- length(unique(J_1)) * coll[[1]]$G_DN
  lengthgo <- length(unique(J_1)) * coll[[1]]$G_go
  lengthneut <- length(unique(J_1)) * coll[[1]]$G_neut


  DOtmpr1 <- asplit(tmpr_1[,c(1:lengthDO)],1)
  DO_r_1 <- lapply(DOtmpr1,function(x) matrix(x,ncol=coll[[1]]$G_DO))

  DNtmpr1 <- asplit(tmpr_1[,c((lengthDO + 1) : (lengthDO + lengthDN))],1)
  DN_r_1 <- lapply(DNtmpr1,function(x) matrix(x,ncol=coll[[1]]$G_DN))

  gotmpr1 <- asplit(tmpr_1[,c((lengthDO + lengthDN + 1) : (lengthDO + lengthDN + lengthgo))],1)
  go_r_1 <- lapply(gotmpr1,function(x) matrix(x,ncol=coll[[1]]$G_go))

  neutmpr1 <- asplit(tmpr_1[,c((lengthDO + lengthDN + lengthgo + 1) : (lengthDO + lengthDN + lengthgo + lengthneut))],1)
  neut_r_1 <- lapply(neutmpr1,function(x) matrix(x,ncol=coll[[1]]$G_neut))

  lengthsam <- 6000

  r1s <- NULL
  grp <- NULL
  for (i in c(1:lengthsam)){

    Do = DO[,i]
    Dn = DN[,i]
    Go = go[,i]
    Neut = neut[,i]

    DOr1 <- rowSums(apply(DO_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_DO)
    DNr1 <- rowSums(apply(DN_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_DN)
    gor1 <- rowSums(apply(go_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_go)
    neutr1 <- rowSums(apply(neut_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_neut)

    samp <- tibble(J_1 = J_1,
                   DO_r1=pnorm(Do + DOr1) * coll[[1]]$X_DO[,1],
                   DN_r1=pnorm(Dn + DNr1) * coll[[1]]$X_DN[,1],
                   go_r1=pnorm(Go + gor1) * coll[[1]]$X_go[,1],
                   neut_r1=pnorm(Neut + neutr1) * coll[[1]]$X_neut[,1],
                   isold = coll[[1]]$X_DO[,1],
                   condition = tidydat$condition) %>%
      #filter(isold==1) %>%
      mutate(it = i)

    samp2 <- tibble(J_1 = J_1,
                    DO_r1=pnorm(Do) * coll[[1]]$X_DO[,1],
                    DN_r1=pnorm(Dn) * coll[[1]]$X_DN[,1],
                    go_r1=pnorm(Go) * coll[[1]]$X_go[,1],
                    neut_r1=pnorm(Neut) * coll[[1]]$X_neut[,1],
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
  neutr1s <- r1s %>% select(neut_r1,J_1,condition,it) %>%
    mutate(condition = "None") %>%
    distinct() %>%
    mutate(par = "neut_r1") %>%
    rename("value" = "neut_r1")


  r1shdi <-bind_rows(as_tibble(bind_rows(DOr1s,DNr1s,gor1s,neutr1s)) %>%
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

                     as_tibble(bind_rows(DOr1s,DNr1s,gor1s,neutr1s)) %>%
                       group_by(J_1,condition,par) %>%
                       median_qi(value,.width = c(.50,.80,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5",
                              quantdesc2!="Median0.8"),
                     as_tibble(bind_rows(DOr1s,DNr1s,gor1s,neutr1s)) %>%
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
  neutgrp <- grp %>% select(neut_r1,condition,it) %>%
    mutate(condition = "None") %>%
    distinct() %>%
    mutate(par = "neut_r1") %>%
    rename("value" = "neut_r1")

  grphdi <-bind_rows(as_tibble(bind_rows(DOgrp,DNgrp,gogrp,neutgrp)) %>%
                       group_by(condition,par) %>%
                       mean_hdi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Mean","HDI","HDI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Mean0.5") %>%
                       filter(quantdesc2!="Mean0.8"),
                     as_tibble(bind_rows(DOgrp,DNgrp,gogrp,neutgrp)) %>%
                       group_by(condition,par) %>%
                       median_qi(value,.width = c(.50,.8,.95)) %>%
                       pivot_longer(cols=c("value",".lower",".upper"),names_to="quantdesc",values_to = "quants") %>%
                       mutate(quantdesc = factor(quantdesc,levels=c("value",".lower",".upper"),
                                                 labels=c("Median","QI","QI"))) %>%
                       mutate(quantdesc2 = paste0(quantdesc,.width)) %>%
                       filter(quantdesc2!="Median0.5") %>%
                       filter(quantdesc2!="Median0.8"),
                     as_tibble(bind_rows(DOgrp,DNgrp,gogrp,neutgrp)) %>%
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


summarizecorrelations <- function(coll){

  if(dim(coll[[2]]$L_1)[2] != 1){

  totL1 <- apply(coll[[2]]$L_1,1,function(x) extendcor2(x,coll[[1]]$ncoef))
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

dprm2 <- c("AP2007_e1","AP2007_e2","AP2007_e3","BKSR2013_e1","NFR2013_e1","NFR2013_e3","GKH1999_e1",
           "JCD2012_e1a","JCM2019_e1","KAWY2013_e4","LBA2019_e1","LM2020_e1","TR2017_e1",
           "WKS2020_e1","KL2012_e2","KL2012_e3", "BTL2013_e1","CSR2015_e1","DR2012_e1b",
           "DW2020_e1","KK2015_e1","KUO2017_e2","MKG2018_e2","OZH2010_e1","PCM2006_e1",
           "PRM2010_e1","RSP2012_e1","SB2020_e1","SBT2018_e1","SCR2019_e1","JW2019_e1",
           "SD2014_e2","SHJ2005_e1","TMP2014_e1","WHD2018_e1")

dprm <- c("APH2016_e1","BSG2014_e1","D2007_e1a","D2007_e1b","FBH2013_e1",
          "FGR2019_e1","FO2016_e2a", "FO2016_e2b","FO2016_e3","GKH1999_e2","GKH1999_e3",
          "GKH1999_e4","HDM2006_e1","HDM2006_e2","HUW2015_e1","JCD2012_e1b","JWH2009_e1",
          "KAWY2013_e2","KAWY2013_e3","KFH2013_e1","KL2012_e1","KL2012_e4","KY2010_e1",
          "KY2011_e1","KY2016_e1","LBA2019_e2", "LP2006_e2","LP2013_e1","NFR2013_e2",
          "OBD2017_e1", "QGM2021_e1","RS2009_e1","RS2009_e2","RSM2009_e1","SB2020_e2",
          "SB2020_e3","SD2004_e2","SD2014_e1", "USS2015_e1","WKH2020_e1",
          "WMS2012_e1","ZMD2011_e1","ZOL2021_e1")


for(model in c("EVSDT","UVSDT","DPSDT","M0SDT","MASDT","2HTM")){
for (i in exps){


  if(model=="2HTM"){
    modeln <- "H2TM"
  } else {
    modeln<-model
  }
  fitfiles <- list.files(paste0(pathfit,modeln,"/"),full.names=T)
  #parfiles <- list.files(paste0(pathpar,model,"/"),full.names=T)
  if(length(fitfiles[grepl(i,fitfiles)])!=0){

   # old <- readRDS(parfiles[grepl(i,parfiles)])

    #test <- summarizeindcrits(coll)
    #li2 <- append(old,list("individual_crit"=test))

  res <- summarizeparameters(i,model)
  saveRDS(res,file=paste0("D:/HS01/ParameterSummary/NoItemEffects/",modeln,"/parametersummary_full_",model,"_",i,".rds"))
}
  }
}
