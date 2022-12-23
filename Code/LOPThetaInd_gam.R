# General information ------------------------

# No item effects, Random crit (i.e., criteria implemented as deviations from grand mean for individuals)
# for models EVgamSDT, UVgamSDT, DPgamSDT, M0gamSDT, MAgamSDT
# functions with 'Full' take individual posterior samples for prediction
# functions with 'LOP' take group-level distributions, wjhere necessary, to sample/approximate individual posterior samples

# for 'Full', the model prediction functions are identical to 'No item effects' set of modles

# tmpdat<-tmpdat_pred
# fit<-fit_f
# model<-"DPgamSDT"
# Top level prediction

# collect standata and posterior samples in a single object for all models -----

prepTheta_Ind <- function(tmpdat,fit, model){

  if (model %in% c("UVgamSDT")){
    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_gamma =  apply(fit$draws("b_gamma")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("crits")[,,],3,rbind),
                       L_1 = apply(fit$draws("L_1")[,,],3,rbind))

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
                       sigma_Z_1 = tmpdat$Z_disc,
                       gamma_Z_1 = tmpdat$Z_gamma,
                       J_2 = tmpdat$J_2)

  } else if (model %in% c("EVgamSDT")){
    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       #b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_gamma =  apply(fit$draws("b_gamma")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("crits")[,,],3,rbind),
                       L_1 = apply(fit$draws("L_1")[,,],3,rbind))

    # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1

    extractdat <- list(J_1 = tmpdat$J_1,
                       nthres = tmpdat$nthres,
                       Y = tmpdat$Y,
                       X = tmpdat$X,
                       #X_disc = tmpdat$X_disc,
                       X_gamma = tmpdat$X_gamma,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       #sigma_Z_1 = tmpdat$Z_disc,
                       gamma_Z_1 = tmpdat$Z_gamma)

  } else if (model %in% c("DPRMgamSDT")){
    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_gamma =  apply(fit$draws("b_gamma")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       mu_s = apply(fit$draws("mu_s")[,,],3,rbind),
                       s_m = apply(fit$draws("s_m")[,,],3,rbind),
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

  } else if (model %in% c("MAgamSDT")){

    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_nonattmu = apply(fit$draws("b_nonattmu")[,,],3,rbind),
                       b_gamma =  apply(fit$draws("b_gamma")[,,],3,rbind),

                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),


                       L_1 = apply(fit$draws("L_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("crits")[,,],3,rbind))

    # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1


    extractdat <- list(J_1 = tmpdat$J_1,
                       nthres = tmpdat$nthres,
                       Y = tmpdat$Y,
                       X = tmpdat$X,
                       X_disc = tmpdat$X_disc,
                       X_nonattmu = tmpdat$X_nonattmu,
                       X_gamma = tmpdat$X_gamma,

                       G = tmpdat$G,
                       G_disc = tmpdat$G_disc,
                       G_nonattmu = tmpdat$G_nonattmu,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       sigma_Z_1 = tmpdat$Z_disc,
                       gamma_Z_1 = tmpdat$Z_gamma,

                       nonattemu_Z_1 = tmpdat$Z_nonattmu
    )
  } else if (model %in% c("M0gamSDT","DPgamSDT")){

    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_gamma =  apply(fit$draws("b_gamma")[,,],3,rbind),
                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),
                       L_1 = apply(fit$draws("L_1")[,,],3,rbind),
                       Intercept = apply(fit$draws("crits")[,,],3,rbind))

    # tmpZ_1 = tmpdat[grepl("Z", names(tmpdat))]
    # Z_1 = matrix(unlist(tmpZ_1),ncol=length(tmpZ_1))
    ncoef = tmpdat$M_1


    extractdat <- list(J_1 = tmpdat$J_1,
                       nthres = tmpdat$nthres,
                       Y = tmpdat$Y,
                       X = tmpdat$X,
                       X_disc = tmpdat$X_disc,
                       X_gamma = tmpdat$X_gamma,

                       G = tmpdat$G,
                       G_disc = tmpdat$G_disc,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       sigma_Z_1 = tmpdat$Z_disc,
                       gamma_Z_1 = tmpdat$Z_gamma
    )
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

  } else if (model %in% c("DPRMSDT")){
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
                       X_nonattmu = tmpdat$X_nonattmu,
                       G = tmpdat$G,
                       G_disc = tmpdat$G_disc,
                       G_nonattmu = tmpdat$G_nonattmu,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       sigma_Z_1 = tmpdat$Z_disc,
                       nonattemu_Z_1 = tmpdat$Z_nonattmu,
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

  } else if (model %in% c("MASDT")){
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
                         Z_DO = tmpdat$Z_DO,
                         Z_DN = tmpdat$Z_DN,
                         Z_go = tmpdat$Z_go)
    }

  }

  list(extractdat,extractfit)

}

# Top-level prediction functions -------------------------------

predictionLOP_Ind <- function(tmpdat,fit,model){

  coll <- prepTheta_Ind(tmpdat,fit,model)

  if (model == "UVgamSDT"){

    out <- UVgamSDTpredictionLOP_Ind(coll)

  } else if (model == "EVgamSDT"){

    out <- EVgamSDTpredictionLOP_Ind(coll)

  } else if (model == "MAgamSDT"){

    out <- MAgamSDTpredictionLOP_Ind(coll)

  } else if (model == "M0gamSDT"){

    out <- M0gamSDTpredictionLOP_Ind(coll)

  } else if (model == "DPgamSDT"){

    out <- DPgamSDTpredictionLOP_Ind(coll)

  }

  out
}

predictionFULL_Ind <- function(tmpdat,fit,model){


  #test = rep(TRUE,tmpdat$N)

  coll <- prepTheta_Ind(tmpdat,fit,model)

  if (model %in% c("UVgamSDT")){

    out <- UVSDTpredictionFULL_Ind(coll)

  } else if (model %in% c("EVgamSDT")){

    out <- EVSDTpredictionFULL_Ind(coll)

  } else if (model == "DPgamSDT"){

    out <- DPSDTpredictionFull_Ind(coll)
  } else if (model %in% c("MAgamSDT")){

    out <- MASDTpredictionFull_Ind(coll)
  } else if (model == "M0gamSDT"){

    out <- M0SDTpredictionFull_Ind(coll)
  }

  out
}


# LOP Model prediction functions ---------------

UVgamSDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)

  gamma_use <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma_use[[m]] <- tmp
  }
  #r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres


  p_sd <- cbind(coll[[2]]$sd_1)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]))
  L_1 <- coll[[2]]$L_1
  totL1 <- apply(L_1,1,function(x) extendcor2(x,ncoef))

  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu# + r2
  lengthsam <- dim(mu)[2]

  res <- list()
  for (j in 1:lengthsam){

    r_1 <- sampleParticipantUVgamSDT(sds=p_sd[j,],
                                     mus = p_mu[j,],
                                     tmpL1=totL1[,j],
                                     ncoef = ncoef,
                                     nthres=nthres,
                                     nJ_1 = nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){



        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- disc_use[J_1==i,j] + as.numeric(r_1[[i]][["sigma_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        gamma_tmp <- gamma_use[[j]][i,] + r_1[[i]][["crit_r_1"]][k,]
        crit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma_tmp,0))),nthres))

        tmp <- makethetaUVSDT_LOP(as.numeric(mu_tmp),
                                  disc_tmp,
                                  crit_tmp)

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}
EVgamSDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  #disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)

  gamma_use <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma_use[[m]] <- tmp
  }
  #r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres


  p_sd <- cbind(coll[[2]]$sd_1)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]))
  L_1 <- coll[[2]]$L_1
  totL1 <- apply(L_1,1,function(x) extendcor2(x,ncoef))

  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu# + r2
  lengthsam <- dim(mu)[2]

  res <- list()
  for (j in 1:lengthsam){

    r_1 <- sampleParticipantEVgamSDT(sds=p_sd[j,],
                                     mus = p_mu[j,],
                                     tmpL1=totL1[,j],
                                     ncoef = ncoef,
                                     nthres=nthres,
                                     nJ_1 = nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){



        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        #disc_tmp <- disc_use[J_1==i,j] + as.numeric(r_1[[i]][["sigma_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        gamma_tmp <- gamma_use[[j]][i,] + r_1[[i]][["crit_r_1"]][k,]
        crit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma_tmp,0))),nthres))

        tmp <- makethetaUVSDT_LOP(as.numeric(mu_tmp),
                                  0,
                                  crit_tmp)

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

MAgamSDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)
  nonattenmu_use <- as.matrix(coll[[1]]$X_nonattmu) %*%  t(coll[[2]]$b_nonattmu)
  gamma_use <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma_use[[m]] <- tmp
  }
  #r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres


  p_sd <- cbind(coll[[2]]$sd_1)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]))


  L_1 <- coll[[2]]$L_1
  totL1 <- apply(L_1,1,function(x) extendcor2(x,ncoef))

  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu# + r2
  lengthsam <- dim(mu)[2]

  res <- list()
  for (j in 1:lengthsam){

    #s_m <- sampleParticipantdirichlet(p_mus[j,],nJ_1)

    r_1 <- sampleParticipantMAgamSDT(sds=p_sd[j,],
                                     mus = p_mu[j,],
                                     tmpL1=totL1[,j],
                                     ncoef = ncoef,
                                     nthres=nthres,
                                     nJ_1 = nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){



        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- as.matrix(coll[[1]]$X_disc[J_1==i,])[,1] * pnorm(disc_use[J_1==i,j] + as.numeric(r_1[[i]][["disc_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,])))

        attenu_tmp <- pnorm(nonattenmu_use[J_1==i,j] + as.numeric(r_1[[i]][["att_r_1"]][k,] %*% t(coll[[1]]$nonattemu_Z_1[J_1==i,])))

        for(m in seq_along(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1])){

          disc_tmp[m] <- ifelse(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1][m] == 0,1,disc_tmp[m])
          attenu_tmp[m] <- ifelse(as.matrix(coll[[1]]$X_nonattmu)[J_1 ==i,1][m] == 0,1,attenu_tmp[m])
        }


        attenmu_use = attenu_tmp * mu_tmp
        gamma_tmp <- gamma_use[[j]][i,] + r_1[[i]][["crit_r_1"]][k,]
        crit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma_tmp,0))),nthres))

        tmp <- makethetaMASDT_LOP(mu_tmp = as.numeric(mu_tmp),
                                  disc_tmp = as.numeric(disc_tmp),
                                  attenmu_use = attenmu_use,

                                  criteria = crit_tmp)


        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}
M0gamSDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)

  gamma_use <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma_use[[m]] <- tmp
  }
  #r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres


  p_sd <- cbind(coll[[2]]$sd_1)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]))


  L_1 <- coll[[2]]$L_1
  totL1 <- apply(L_1,1,function(x) extendcor2(x,ncoef))

  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu# + r2
  lengthsam <- dim(mu)[2]

  res <- list()
  for (j in 1:lengthsam){

    #s_m <- sampleParticipantdirichlet(p_mus[j,],nJ_1)

    r_1 <- sampleParticipantUVgamSDT(sds=p_sd[j,],
                                     mus = p_mu[j,],
                                     tmpL1=totL1[,j],
                                     ncoef = ncoef,
                                     nthres=nthres,
                                     nJ_1 = nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){



        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- as.matrix(coll[[1]]$X_disc[J_1==i,])[,1] * pnorm(disc_use[J_1==i,j] + as.numeric(r_1[[i]][["sigma_r_1"]][k,] %*% t(coll[[1]]$sigma_Z_1[J_1==i,])))

        for(m in seq_along(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1])){

          disc_tmp[m] <- ifelse(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1][m] == 0,1,disc_tmp[m])

        }



        gamma_tmp <- gamma_use[[j]][i,] + r_1[[i]][["crit_r_1"]][k,]
        crit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma_tmp,0))),nthres))

        tmp <- makethetaM0SDT_LOP(mu_tmp = as.numeric(mu_tmp),
                                  disc_tmp = as.numeric(disc_tmp),

                                  criteria = crit_tmp)


        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}
DPgamSDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)

  gamma_use <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma_use[[m]] <- tmp
  }
  #r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres


  p_sd <- cbind(coll[[2]]$sd_1)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]))
  #p_mus <- coll[[2]]$mu_s

  L_1 <- coll[[2]]$L_1
  totL1 <- apply(L_1,1,function(x) extendcor2(x,ncoef))

  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu# + r2
  lengthsam <- dim(mu)[2]

  res <- list()
  for (j in 1:lengthsam){

    #s_m <- sampleParticipantdirichlet(p_mus[j,],nJ_1)

    r_1 <- sampleParticipantUVgamSDT(sds=p_sd[j,],
                                     mus = p_mu[j,],
                                     tmpL1=totL1[,j],
                                     ncoef = ncoef,
                                     nthres=nthres,
                                     nJ_1 = nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){



        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- as.matrix(coll[[1]]$X_disc[J_1==i,])[,1] * pnorm(disc_use[J_1==i,j] + as.numeric(r_1[[i]][["sigma_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,])))
        #sm_tmp <- s_m[[i]]$sm_r_1[k,]
        gamma_tmp <- gamma_use[[j]][i,] + r_1[[i]][["crit_r_1"]][k,]
        crit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma_tmp,0))),nthres))

        tmp <- makethetaDPSDT_LOP(mu_tmp = as.numeric(mu_tmp),
                                  disc_tmp = as.numeric(disc_tmp),
                                  criteria = crit_tmp)

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

# Deviance model prediction funcions --------------------

UVSDTpredictionFULL_Ind <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  J_1 <- coll[[1]]$J_1

  uniqueid <- length(unique(J_1))
  tmpr_1 <- coll[[2]]$r_1


  mutmpr1 <- asplit(tmpr_1[,c(1:((ncoef-nthres)/2 * uniqueid))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  sigmatmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/2 * uniqueid)+ 1):((ncoef-nthres)/2 * uniqueid*2))],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])


  for (i in c(1:lengthsam)){



    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1[,i] <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
  }

  #prep
  mu_use = mu + mur1
  disc_use = disc + sigmar1


  res <- list()
  for (j in 1:lengthsam){


    res[[j]]  <- makethetaUVSDT(mu_tmp=as.numeric(mu_use[,j]),
                                disc_tmp=as.numeric(disc_use[,j]),
                                criteria =matrix(coll[[2]]$Intercept[j,],
                                                 ncol = nthres,byrow=F)[J_1,])

  }
  res
}


DPSDTpredictionFull_Ind <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  J_1 <- coll[[1]]$J_1
  uniqueid <- length(unique(J_1))
  tmpr_1 <- coll[[2]]$r_1


  mutmpr1 <- asplit(tmpr_1[,c(1:((ncoef-nthres)/2 * uniqueid))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  sigmatmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/2 * uniqueid)+ 1):((ncoef-nthres)/2 * uniqueid*2))],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))


  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])


  for (i in c(1:lengthsam)){

    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1[,i] <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$sigma_Z_1)
  }

  #prep
  mu_use = mu + mur1
  disc_use =  as.matrix(coll[[1]]$X_disc)[,1] * pnorm(disc + sigmar1)



  res <- list()
  for (j in 1:lengthsam){


    res[[j]]  <- makethetaDPSDT(mu_tmp=as.numeric(mu_use[,j]),
                                disc_tmp=as.numeric(disc_use[,j]),
                                criteria =matrix(coll[[2]]$Intercept[j,],
                                                 ncol = nthres,byrow=F)[J_1,])

  }
  res
}
M0SDTpredictionFull_Ind <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  J_1 <- coll[[1]]$J_1
  uniqueid <- length(unique(J_1))
  tmpr_1 <- coll[[2]]$r_1


  mutmpr1 <- asplit(tmpr_1[,c(1:((ncoef-nthres)/2 * uniqueid))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  sigmatmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/2 * uniqueid)+ 1):((ncoef-nthres)/2 * uniqueid*2))],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))


  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])


  for (i in c(1:lengthsam)){



    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1[,i] <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
  }

  #prep
  mu_use = mu + mur1
  disc_use =  pnorm(disc + sigmar1)

  for (i in c(1:lengthsam)){
    disc_use[,i] =  ifelse(as.matrix(coll[[1]]$X_disc)[,1] == 0,1,disc_use[,i])
  }

  res <- list()
  for (j in 1:lengthsam){



    res[[j]]  <- makethetaM0SDT(mu_tmp=as.numeric(mu_use[,j]),
                                disc_tmp=as.numeric(disc_use[,j]),
                                criteria =matrix(coll[[2]]$Intercept[j,],
                                                 ncol = nthres,byrow=F)[J_1,])

  }
  res
}


MASDTpredictionFull_Ind <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)
  nonattenmu <- as.matrix(coll[[1]]$X_nonattmu) %*%  t(coll[[2]]$b_nonattmu)

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  J_1 <- coll[[1]]$J_1
  uniqueid <- length(unique(J_1))
  tmpr_1 <- coll[[2]]$r_1


  mutmpr1 <- asplit(tmpr_1[,c(1:((ncoef-nthres)/3* uniqueid))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/3))

  sigmatmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/3 * uniqueid)+ 1):((ncoef-nthres)/3 * uniqueid*2))],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/3))


  nonattmutmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres)/3 * uniqueid *2)+ 1):((ncoef-nthres)/3 * uniqueid*3))],1)
  nonattmu_r_1 <- lapply(nonattmutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/3))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  nonattmur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])


  for (i in c(1:lengthsam)){



    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1[,i] <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    nonattmur1[,i]<- rowSums(apply(nonattmu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$nonattemu_Z_1)

  }

  #prep
  mu_use = mu + mur1
  disc_use =  pnorm(disc + sigmar1)
  attenmu_use_f = pnorm(nonattenmu + nonattmur1)


  for (i in c(1:lengthsam)){
    disc_use[,i] =  ifelse(as.matrix(coll[[1]]$X_disc)[,1] == 0,1,disc_use[,i])
    attenmu_use_f[,i] =  ifelse(as.matrix(coll[[1]]$X_disc)[,1] == 0,1,attenmu_use_f[,i])
  }

  attenmu_use = attenmu_use_f * mu_use

  res <- list()
  for (j in 1:lengthsam){



    res[[j]]  <- makethetaMASDT(mu_tmp=as.numeric(mu_use[,j]),
                                disc_tmp=as.numeric(disc_use[,j]),
                                criteria =matrix(coll[[2]]$Intercept[j,],
                                                 ncol = nthres,byrow=F)[J_1,],
                                munon_tmp =as.numeric(attenmu_use[,j]))

  }
  res
}

EVSDTpredictionFULL_Ind <- function(coll){


  Y <- coll[[1]]$Y



  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  #mur2 <- apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2]) #+ item effects
  #disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres


  J_1 <- coll[[1]]$J_1
  uniqueid <- length(unique(J_1))
  tmpr_1 <- coll[[2]]$r_1


  mutmpr1 <- asplit(tmpr_1[,c(1:((ncoef-nthres) * uniqueid))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  #sigmar1 <- matrix(ncol = 6000,nrow=length(Y))


  for (i in c(1:lengthsam)){


    # i<-1
    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1]) * coll[[1]]$mu_Z_1)
    #sigmar1[,i] <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
  }


  #prep
  mu_use = mu + mur1 #+ mur2



  res <- list()
  for (j in 1:lengthsam){

    res[[j]]  <- makethetaUVSDT(mu_tmp=as.numeric(mu_use[,j]),
                                disc_tmp=rep(0,length(coll[[1]]$J_1)),
                                criteria =matrix(coll[[2]]$Intercept[j,],
                                                 ncol = nthres,byrow=F)[J_1,])



  }
  res
}

# sample deviations from group-level for LOP ----------

# for all SDT models, sample distributional parameters + criteria from
# MVN(0,Sigma_alpha) where Sigma_alpha is the covariance matrix
# criteria then are transformed from the real-line via softmax to be ordered in normal space

# 2 parameters + thresholds
sampleParticipantUVgamSDT <- function(sds,mus,tmpL1,ncoef,nthres,nJ_1){

  # sds <- p_sd[j,]
  #  mus <- p_mu[j,]
  #  tmpL1 <- totL1[j,]
  sepncoef <- (ncoef-nthres)/2


  L_1 <- tmpL1
  #matrix(tmpL1,nrow=length(sds),byrow=F) %*% t(matrix(tmpL1,nrow=length(sds),byrow=F))

  # turn L_1 cholesky into correlation matrix by L_1 %*% t(L_1)
  # turn sigma and correlation matrix into covariance by S %*% Corr %*% S
  # Make full covmat for participant effects on mu/sigma and criteria for single sample

  ind <- NULL
  for(i in c(1:nJ_1)){

    r_1 <- mvtnorm::rmvnorm(1000,
                            mean = mus,
                            sigma= outer(sds,sds) * L_1) #equiv:diag(sds) %*% L_1 %*% diag(sds))

    # only use samples where criteria samples are ordered
    #tmpr_1_use <- r_1[apply(r_1, 1, function(x) !is.unsorted(tail(x,nthres))),]
    #r_1_use<-tmpr_1_use[1:min(1000,dim(tmpr_1_use)[1]),]

    ind[[i]] <- list(mu_r_1 = matrix(r_1[,c(1:sepncoef)],ncol=sepncoef),
                     sigma_r_1 = matrix(r_1[,c((sepncoef + 1):(ncoef-nthres))],ncol=sepncoef),
                     crit_r_1 = r_1[,c(((ncoef-nthres)+1):dim(r_1)[2])])
  }

  ind
}
# 1 parameter + thresholds
sampleParticipantEVgamSDT <- function(sds,mus,tmpL1,ncoef,nthres,nJ_1){

  # sds <- p_sd[j,]
  #  mus <- p_mu[j,]
  #  tmpL1 <- totL1[j,]
  sepncoef <- (ncoef-nthres)


  L_1 <- tmpL1
  #matrix(tmpL1,nrow=length(sds),byrow=F) %*% t(matrix(tmpL1,nrow=length(sds),byrow=F))

  # turn L_1 cholesky into correlation matrix by L_1 %*% t(L_1)
  # turn sigma and correlation matrix into covariance by S %*% Corr %*% S
  # Make full covmat for participant effects on mu/sigma and criteria for single sample

  ind <- NULL
  for(i in c(1:nJ_1)){

    r_1 <- mvtnorm::rmvnorm(1000,
                            mean = mus,
                            sigma= outer(sds,sds) * L_1) #equiv:diag(sds) %*% L_1 %*% diag(sds))

    # only use samples where criteria samples are ordered
    #tmpr_1_use <- r_1[apply(r_1, 1, function(x) !is.unsorted(tail(x,nthres))),]
    #r_1_use<-tmpr_1_use[1:min(1000,dim(tmpr_1_use)[1]),]

    ind[[i]] <- list(mu_r_1 = matrix(r_1[,c(1:sepncoef)],ncol=sepncoef),
                     crit_r_1 = r_1[,c((sepncoef+1):dim(r_1)[2])])
  }

  ind
}
# 3 parameters + thresholds
sampleParticipantMAgamSDT <- function(sds,mus,tmpL1,ncoef,nthres,nJ_1){

  # sds <- p_sd[j,]
  #  mus <- p_mu[j,]
  #  tmpL1 <- totL1[j,]
  sepncoef <- (ncoef-nthres)/3


  L_1 <- tmpL1
  #matrix(tmpL1,nrow=length(sds),byrow=F) %*% t(matrix(tmpL1,nrow=length(sds),byrow=F))

  # turn L_1 cholesky into correlation matrix by L_1 %*% t(L_1)
  # turn sigma and correlation matrix into covariance by S %*% Corr %*% S
  # Make full covmat for participant effects on mu/sigma and criteria for single sample

  ind <- NULL
  for(i in c(1:nJ_1)){

    r_1 <- mvtnorm::rmvnorm(1000,
                            mean = mus,
                            sigma= outer(sds,sds) * L_1) #equiv:diag(sds) %*% L_1 %*% diag(sds))

    # only use samples where criteria samples are ordered
    #tmpr_1_use <- r_1[apply(r_1, 1, function(x) !is.unsorted(tail(x,nthres))),]
    #r_1_use<-tmpr_1_use[1:min(1000,dim(tmpr_1_use)[1]),]

    ind[[i]] <- list(mu_r_1 = matrix(r_1[,c(1:sepncoef)],ncol=sepncoef),
                     disc_r_1 = matrix(r_1[,c((sepncoef + 1):(sepncoef + sepncoef))],ncol=sepncoef),
                     att_r_1 =  matrix(r_1[,c((sepncoef + sepncoef + 1):(ncoef - nthres))],ncol=sepncoef),
                     crit_r_1 = r_1[,c((ncoef - nthres +1):dim(r_1)[2])])
  }

  ind
}
# sample response mapping from Dirichlet
sampleParticipantdirichlet <- function(mus,nJ_1){


  ind <- NULL
  for(i in c(1:nJ_1)){
    sm_sample <- DirichletReg::rdirichlet(n = 1000,alpha=mus)


    ind[[i]] <- list(sm_r_1 =sm_sample)
  }

  ind

}





# Calculation of Theta and likelihoods -------------------------------
# Make Theta
makethetaUVSDT <- function(mu_tmp,disc_tmp,criteria){


  theta <- matrix(ncol=(dim(criteria)[[2]] + 1),nrow=length(mu_tmp))
  crits <- cbind(-Inf,criteria,Inf)

  for(Y in c(1:(dim(crits)[[2]]-1))){

    theta[,Y] <- pnorm(exp(disc_tmp) * (crits[,Y+1] - mu_tmp)) -
      pnorm(exp(disc_tmp) * (crits[,Y] - mu_tmp))

  }

  theta

}

makethetaDPSDT <- function(mu_tmp,disc_tmp,criteria){


  theta <- matrix(ncol=(dim(criteria)[[2]] + 1),nrow=length(mu_tmp))
  crits <- cbind(-Inf,criteria,Inf)

  for(Y in c(1:(dim(crits)[[2]]-1))){
    if(Y == (dim(crits)[[2]]-1)){
      theta[,Y] <- disc_tmp + (1 - disc_tmp) * (pnorm(crits[,Y+1] - mu_tmp) -
                                                  pnorm(crits[,Y] - mu_tmp))
    } else {
      theta[,Y] <- (1 - disc_tmp) * (pnorm(crits[,Y+1] - mu_tmp) -
                                       pnorm(crits[,Y] - mu_tmp))

    }

  }

  theta

}

makethetaMASDT <- function(mu_tmp,disc_tmp,criteria,munon_tmp){


  theta <- matrix(ncol=(dim(criteria)[[2]] + 1),nrow=length(mu_tmp))
  crits <- cbind(-Inf,criteria,Inf)

  for(Y in c(1:(dim(crits)[[2]]-1))){

    theta[,Y] <- disc_tmp * (pnorm(crits[,Y+1] - mu_tmp) - pnorm(crits[,Y] - mu_tmp)) +
      (1 - disc_tmp) * (pnorm(crits[,Y+1] - munon_tmp) - pnorm(crits[,Y] - munon_tmp))
  }

  theta

}

makethetaM0SDT <- function(mu_tmp,disc_tmp,criteria){


  theta <- matrix(ncol=(dim(criteria)[[2]] + 1),nrow=length(mu_tmp))
  crits <- cbind(-Inf,criteria,Inf)

  for(Y in c(1:(dim(crits)[[2]]-1))){

    theta[,Y] <- disc_tmp * (pnorm(crits[,Y+1] - mu_tmp) - pnorm(crits[,Y] - mu_tmp)) +
      (1 - disc_tmp) * (pnorm(crits[,Y+1]) - pnorm(crits[,Y]))
  }

  theta

}

makethetaDPSDT_LOP <- function(mu_tmp,disc_tmp,criteria){

  theta <- matrix(ncol=(length(criteria) + 1),nrow=length(mu_tmp))
  crits <- c(-Inf,criteria,Inf)

  for(Y in c(1:(length(crits)-1))){

    if(Y == (length(crits)-1)){
      theta[,Y] <- disc_tmp + (1 - disc_tmp) * (pnorm(crits[Y+1] - mu_tmp) -
                                                  pnorm(crits[Y] - mu_tmp))
    } else {
      theta[,Y] <- (1 - disc_tmp) * (pnorm(crits[Y+1] - mu_tmp) -
                                       pnorm(crits[Y] - mu_tmp))

    }


  }

  theta

}


makethetaUVSDT_LOP <- function(mu_tmp,disc_tmp,criteria){

  theta <- matrix(ncol=(length(criteria) + 1),nrow=length(mu_tmp))
  crits <- c(-Inf,criteria,Inf)

  for(Y in c(1:(length(crits)-1))){

    theta[,Y] <- pnorm(exp(disc_tmp) * (crits[Y+1] - mu_tmp)) -
      pnorm(exp(disc_tmp) * (crits[Y] - mu_tmp))

  }

  theta

}
makethetaMASDT_LOP <- function(mu_tmp,disc_tmp,attenmu_use,criteria){


  theta <- matrix(ncol=(length(criteria) + 1),nrow=length(mu_tmp))
  crits <- c(-Inf,criteria,Inf)

  for(Y in c(1:(length(crits)-1))){

    theta[,Y] <- disc_tmp * (pnorm(crits[Y+1] - mu_tmp) - pnorm(crits[Y] - mu_tmp)) +
      (1 - disc_tmp) * (pnorm(crits[Y+1] - attenmu_use) - pnorm(crits[Y] - attenmu_use))
  }

  theta

}
makethetaM0SDT_LOP <-  function(mu_tmp,disc_tmp,criteria){


  theta <- matrix(ncol=(length(criteria) + 1),nrow=length(mu_tmp))
  crits <- c(-Inf,criteria,Inf)

  for(Y in c(1:(length(crits)-1))){

    theta[,Y] <- disc_tmp * (pnorm(crits[Y+1] - mu_tmp) - pnorm(crits[Y] - mu_tmp)) +
      (1 - disc_tmp) * (pnorm(crits[Y+1]) - pnorm(crits[Y]))
  }

  theta

}





# Deviance from probabilities ----------------

getDeviancefromTheta_Ind <- function(obs,thetaList){

  # thetaList<-pred2[1]
  thetaList<-lapply(thetaList,function(x) replace(x,x==0,1e-8))
  -2 * (lapply(lapply(lapply(thetaList,FUN=log),FUN=function(x) x * obs), FUN=sum) %>% unlist())

  # test <-lapply(lapply(thetaList,FUN=log),FUN=function(x) x * obs)[[1]]
  # srat <-matrix(c(1:length(obs),obs),byrow=F,ncol=2)
  # probs <- lapply(thetaList, FUN = `[`, srat)
  # -2 * (lapply(lapply(probs,FUN=log),FUN = sum) %>% unlist())

}

getDeviancefromTheta_Ind_individual <- function(obs,thetaList,J_1){
  # obs<-tmpdat_pred$Y
  #J_1 <- ids$id
  # thetaList <- pred

  probs<-lapply(thetaList,function(x) replace(x,x==0,1e-8))


  lprobs <- lapply(lapply(lapply(probs,FUN=log),FUN=function(x) x * obs),function(x) rowSums(x))

  allprob <- NULL
  for(i in c(1:length(lprobs))){

    provs <- tibble(prob = lprobs[[i]],
                    id = J_1) %>%
      group_by(id) %>%
      summarize(sumprob = -2 * sum(prob))

    allprob <- allprob %>% bind_rows(provs)
  }

  return(allprob)
}

# point-wise LL ----------------------

getpwLL_Ind <- function(ids,thetaList){

  thetaList<-lapply(lapply(thetaList,function(x) replace(x,x==0,1e-8)),FUN=log)

  datas <- ids %>% ungroup() %>% select(-exp,-manipulation,-condition,-LOPFolds,-id,-isold)

  lists <- list()
  for(j in c(1:6000)){

    lists[[j]] <- rep(as.vector(as.matrix(thetaList[[j]])),as.vector(as.matrix(datas)))
  }

  lists
}

# posterior p-values -----------------

T1statX2 <- function(n,nhat){
  (n-nhat)^2/nhat
}

subsetdf_ind <- function(n,J_1,conditions){
  data.frame(n) %>% mutate(id = J_1,
                           cond = conditions) %>%
    group_by(id,cond) %>%
    summarize_all(sum) %>%
    select(-cond) %>% ungroup()
}

subsetdf_exp <- function(n,J_1,conditions){

  data.frame(n) %>% mutate(id = J_1,
                           cond = conditions) %>%
    group_by(id,cond) %>%
    summarize_all(sum) %>%
    ungroup() %>%
    select(-id) %>%
    group_by(cond) %>%
    summarize_all(mean) %>%
    select(-cond) %>% ungroup()

}


T1variouslevel <- function(theta,exp,model){


    prepdata <- selectModel(exp,model)
    expd <- prepdata$tmpdat
    expd2 <- makeFreqDataInd(exp,model)

    J_1 <- expd2$id # participant number per trial
    conditions <- paste0(expd2$condition,"_",expd2$isold) # condition per trial (_isold for within-subject/within-block conditions)
    Y <- expd$Y # ratings per trial

    theta[theta == 0] <- 1e-8

    nhat <- rowSums(Y) * theta
    npred <- extraDistr::rmnom(dim(theta)[1],rowSums(Y),theta)
    nobs <- Y


  # individual level T1:

  npred_ind <-subsetdf_ind(npred,J_1,conditions)%>% nest(npred = !id)
  nobs_ind <- subsetdf_ind(nobs,J_1,conditions) %>% nest(nobs = !id)
  nhat_ind <- subsetdf_ind(nhat,J_1,conditions)%>% nest(nhat = !id)

  ind_subj<-left_join(left_join(npred_ind,nobs_ind),nhat_ind) %>%
    group_by(id) %>%
    mutate(indT1preds = map2(npred,nhat, ~sum(T1statX2(.x,.y))),
           indT1obs = map2(nobs,nhat, ~sum(T1statX2(.x,.y)))) %>%
    select(-npred,-nobs,-nhat) %>%
    unnest(cols = c(indT1preds, indT1obs))


  #exp-level T1

  npred_exp<-subsetdf_exp(npred,J_1,conditions)
  nobs_exp <- subsetdf_exp(nobs,J_1,conditions)
  nhat_exp <- subsetdf_exp(nhat,J_1,conditions)

  expT1preds <- sum(T1statX2(npred_exp,nhat_exp))
  expT1obs <- sum(T1statX2(nobs_exp,nhat_exp))



    T1 <- ind_subj %>% bind_cols(expT1preds=expT1preds,
                                 expT1obs = expT1obs)


  return(list(T1 = T1,
              avpredfreq = data.frame(npred_exp)))

}




# Helper functions ------------

extendcor <- function(singleL_1,ncoef,nthres){

  cor_coef <- matrix(singleL_1,nrow=ncoef,byrow=F) %*% t(matrix(singleL_1,nrow=ncoef,byrow=F))

  zeromat <- matrix(0,nthres,ncoef)

  cbind(rbind(cor_coef,zeromat),rbind(t(zeromat),diag(nthres)))

}

extendcor2 <- function(singleL_1,ncoef){

  matrix(singleL_1,nrow=ncoef,byrow=F) %*% t(matrix(singleL_1,nrow=ncoef,byrow=F))

  #zeromat <- matrix(0,nthres,ncoef)

  #cbind(rbind(cor_coef,zeromat),rbind(t(zeromat),diag(nthres)))

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

