# General information ------------------------

# Item effects, criteria sampled as deviations from grand mean
# for models EVgamSDT, UVgamSDT, DPgamSDT, M0gamSDT, MAgamSDT, 2HTM
# functions with 'Full' take individual posterior samples for prediction
# functions with 'LOP' take group-level distributions, where necessary, to sample/approximate individual posterior samples


# tmpdat<-tmpdat_pred
# fit<-fit_f
# model<-modelpred
# Top level prediction

# collect standata and posterior samples in a single object for all models -----


prepTheta <- function(tmpdat,fit, model){

  if (model %in% c("UVgamSDT")){

    extractfit <- list(b =  apply(fit$draws("b")[,,],3,rbind),
                       b_disc =  apply(fit$draws("b_disc")[,,],3,rbind),
                       b_gamma =  apply(fit$draws("b_gamma")[,,],3,rbind),

                       sd_1 = apply(fit$draws("sd_1")[,,],3,rbind),
                       r_1 = apply(fit$draws("r_1")[,,],3,rbind),

                       Intercept = apply(fit$draws("crits")[,,],3,rbind),

                       L_1 = apply(fit$draws("L_1")[,,],3,rbind),

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
                       Z_2_mu = tmpdat$Z_2_mu,
                       Z_2_disc = tmpdat$Z_2_disc,

                       Z_2_nonattmu = tmpdat$Z_2_nonattmu)
  } else if (model %in% c("EVgamSDT")){

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

# Top-level prediction functions -------------------------------

# for out-of-sample LOP prediction
predictionLOP <- function(tmpdat,fit,model){


  coll <- prepTheta(tmpdat,fit,model)

  if (model == "UVgamSDT"){

    out <- UVSDTpredictionLOP(coll)

  } else if (model == "EVgamSDT"){

    out <- EVSDTpredictionLOP(coll)

  } else if (model == "DPSDT"){

    out <- DPSDTpredictionLOP(coll)

  } else if (model == "M0gamSDT"){

    out <- M0SDTpredictionLOP(coll)

  } else if (model == "2HTM"){
    out <- H2TMpredictionLOP(coll)
  }


  out
}
# for deviance and out-of-sample KFCV prediction
predictionFULL <- function(tmpdat,fit,model){


  # tmpdat<-tmpdat_pred
  # fit<-fit_f
  # model<-"Gumbel"

  coll <- prepTheta(tmpdat,fit,model)

  if (model == "UVgamSDT"){

    out <- UVSDTpredictionFull(coll)

  } else if (model == "EVgamSDT"){

    out <- EVSDTpredictionFull(coll)

  } else if (model == "Gumbel"){

    out <- GumbelpredictionFULL(coll)
  } else if (model == "DPgamSDT"){

    out <- DPSDTpredictionFull(coll)
  } else if (model == "M0gamSDT"){

    out <- M0SDTpredictionFull(coll)
  }else if (model == "MAgamSDT"){

    out <- MASDTpredictionFull(coll)
  }else if (model == "2HTM"){

    out <- H2TMpredictionFull(coll)
  }

  out
}



# LOP Model prediction functions ---------------
UVSDTpredictionLOP <- function(coll){


  r2 <-  apply(coll[[2]]$r_2_mu,1,function(x) x[coll[[1]]$J_2])


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

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  L_1 <- coll[[2]]$L_1

  p_sd <- cbind(coll[[2]]$sd_1)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]))
  L_1 <- coll[[2]]$L_1
  totL1 <- apply(L_1,1,function(x) extendcor2(x,ncoef))


  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu + r2
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
        disc_tmp <- disc_use[J_1==i,j] + as.numeric(r_1[[i]][["disc_r_1"]][k,] %*% t(coll[[1]]$disc_Z_1[J_1==i,]))
        gamma_tmp <- gamma_use[[j]][i,] + r_1[[i]][["crit_r_1"]][k,]
        crit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma_tmp,0))),nthres))

        tmp <- makethetaUVSDT_LOP(mu_tmp=as.numeric(mu_tmp),
                                  disc_tmp = disc_tmp,
                                  criteria = crit_tmp,
                                  Yactual = Y[J_1==i])

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp))
      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

EVSDTpredictionLOP <- function(coll){


  r2 <-  apply(coll[[2]]$r_2_mu,1,function(x) x[coll[[1]]$J_2])


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

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  L_1 <- coll[[2]]$L_1

  p_sd <- cbind(coll[[2]]$sd_1)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]))
  L_1 <- coll[[2]]$L_1
  totL1 <- apply(L_1,1,function(x) extendcor2(x,ncoef))


  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu + r2
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
        #disc_tmp <- disc_use[J_1==i,j] + as.numeric(r_1[[i]][["disc_r_1"]][k,] %*% t(coll[[1]]$disc_Z_1[J_1==i,]))
        gamma_tmp <- gamma_use[[j]][i,] + r_1[[i]][["crit_r_1"]][k,]
        crit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma_tmp,0))),nthres))

        tmp <- makethetaUVSDT_LOP(mu_tmp=as.numeric(mu_tmp),
                                  disc_tmp = as.numeric(0),
                                  criteria = crit_tmp,
                                  Yactual = Y[J_1==i])

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp))
      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

M0SDTpredictionLOP <- function(coll){



  Y <- coll[[1]]$Y

  mu_use <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b) +
    apply(coll[[2]]$r_2_mu,1,function(x) x[coll[[1]]$J_2]) #+ item effects
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)+
    apply(coll[[2]]$r_2_disc,1,function(x) x[coll[[1]]$J_2])


  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers



  gamma_use <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma_use[[m]] <- tmp
  }

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  L_1 <- coll[[2]]$L_1

  p_sd <- cbind(coll[[2]]$sd_1)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]))
  L_1 <- coll[[2]]$L_1
  totL1 <- apply(L_1,1,function(x) extendcor2(x,ncoef))

  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]


  lengthsam <- dim(mu_use)[2]

  res <- list()
  for (j in 1:lengthsam){

    r_1 <- sampleParticipantUVgamSDT(p_sd[j,],p_mu[j,],totL1[,j],ncoef,nthres,nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){



        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- as.matrix(coll[[1]]$X_disc)[J_1==i,1] * pnorm(disc_use[J_1==i,j] + as.numeric(r_1[[i]][["disc_r_1"]][k,] %*% t(coll[[1]]$disc_Z_1[J_1==i,])))
        gamma_tmp <- gamma_use[[j]][i,] + r_1[[i]][["crit_r_1"]][k,]
        crit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma_tmp,0))),nthres))

        for(m in seq_along(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1])){

          disc_tmp[m] <- ifelse(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1][m] == 0,1,disc_tmp[m])

        }


        tmp <- makethetaM0SDT_LOP(mu_tmp=as.numeric(mu_tmp),
                                  disc_tmp = as.numeric(disc_tmp),
                                  criteria =crit_tmp,
                                  Yactual = Y[J_1==i])

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp))
      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

DPSDTpredictionLOP <- function(coll){

  Y <- coll[[1]]$Y

  mu_use <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b) +
    apply(coll[[2]]$r_2_mu,1,function(x) x[coll[[1]]$J_2]) #+ item effects
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)+
    apply(coll[[2]]$r_2_disc,1,function(x) x[coll[[1]]$J_2])


  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  L_1 <- coll[[2]]$L_1

  p_sd <- cbind(coll[[2]]$sd_1,coll[[2]]$sigma_cr)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]),
                coll[[2]]$mu_cr)

  totL1 <- apply(L_1,1,function(x) extendcor(x,ncoef,nthres))

  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]


  lengthsam <- dim(mu_use)[2]

  res <- list()
  for (j in 1:lengthsam){

    r_1 <- sampleParticipantUVSDT(p_sd[j,],p_mu[j,],totL1[,j],ncoef,nthres,nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){



        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- as.matrix(coll[[1]]$X_disc)[J_1==i,1] * pnorm(disc_use[J_1==i,j] + as.numeric(r_1[[i]][["disc_r_1"]][k,] %*% t(coll[[1]]$disc_Z_1[J_1==i,])))

        tmp <- makethetaDPSDT_LOP(mu_tmp=as.numeric(mu_tmp),
                                  disc_tmp = as.numeric(disc_tmp),
                                  criteria = as.numeric(r_1[[i]][["crit_r_1"]][k,]),
                                  Yactual = Y[J_1==i])

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp))
      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

MASDTpredictionLOP <- function(coll){



  Y <- coll[[1]]$Y

  mu_use <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b) +
    apply(coll[[2]]$r_2_mu,1,function(x) x[coll[[1]]$J_2]) #+ item effects
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)+
    apply(coll[[2]]$r_2_disc,1,function(x) x[coll[[1]]$J_2])
  nonattenmu_use <- as.matrix(coll[[1]]$X_nonattmu) %*%  t(coll[[2]]$b_nonattmu) +
   apply(coll[[2]]$r_2_nonattmu,1,function(x) x[coll[[1]]$J_2])


  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers



  gamma_use <- vector(mode = "list", length = 6000)

  for(m in c(1:6000)){
    tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
    gamma_use[[m]] <- tmp
  }

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  L_1 <- coll[[2]]$L_1

  p_sd <- cbind(coll[[2]]$sd_1)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]))
  L_1 <- coll[[2]]$L_1
  totL1 <- apply(L_1,1,function(x) extendcor2(x,ncoef))

  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]


  lengthsam <- dim(mu_use)[2]

  res <- list()
  for (j in 1:lengthsam){

    r_1 <- sampleParticipantMAgamSDT(p_sd[j,],p_mu[j,],totL1[,j],ncoef,nthres,nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){



        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- as.matrix(coll[[1]]$X_disc)[J_1==i,1] * pnorm(disc_use[J_1==i,j] + as.numeric(r_1[[i]][["disc_r_1"]][k,] %*% t(coll[[1]]$disc_Z_1[J_1==i,])))
        attenu_tmp <- pnorm(nonattenmu_use[J_1==i,j] + as.numeric(r_1[[i]][["att_r_1"]][k,] %*% t(coll[[1]]$nonattemu_Z_1[J_1==i,])))
        gamma_tmp <- gamma_use[[j]][i,] + r_1[[i]][["crit_r_1"]][k,]
        crit_tmp <- 2*qnorm(head(cumsum(softmax(c(gamma_tmp,0))),nthres))

        for(m in seq_along(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1])){

          disc_tmp[m] <- ifelse(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1][m] == 0,1,disc_tmp[m])
          attenu_tmp[m] <- ifelse(as.matrix(coll[[1]]$X_nonattmu)[J_1 ==i,1][m] == 0,1,attenu_tmp[m])
        }


        attenmu_use = attenu_tmp * mu_tmp

        tmp <- makethetaMASDT_LOP(mu_tmp=as.numeric(mu_tmp),
                                  disc_tmp = as.numeric(disc_tmp),
                                  attenmu_use = attenmu_use,
                                  criteria =crit_tmp,
                                  Yactual = Y[J_1==i])

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp))
      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

H2TMpredictionLOP <- function(coll){

  r2DO <-  apply(coll[[2]]$r_2_DO,1,function(x) x[coll[[1]]$J_2])
  r2DN <-  apply(coll[[2]]$r_2_DN,1,function(x) x[coll[[1]]$J_2])
  r2go <-  apply(coll[[2]]$r_2_go,1,function(x) x[coll[[1]]$J_2])

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  DO <- as.matrix(coll[[1]]$X_DO) %*% t(coll[[2]]$b_DO)
  DN <- as.matrix(coll[[1]]$X_DN) %*%  t(coll[[2]]$b_DN)
  go <- as.matrix(coll[[1]]$X_go) %*%  t(coll[[2]]$b_go)

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres


  p_sd <- cbind(coll[[2]]$sd_1)
  p_mu <- matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1])
  L_1 <- coll[[2]]$L_1
  totL1 <- apply(L_1,1,function(x) extendcor(x,ncoef,0))

  mus <- coll[[2]]$mu_s
  mu_a_o <- coll[[2]]$mu_a_o
  mu_a_n <- coll[[2]]$mu_a_n
  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  DO_use = DO + r2DO
  DN_use = DN + r2DN
  go_use = go + r2go

  lengthsam <- dim(DO)[2]

  res <- list()
  for (j in 1:lengthsam){

    r_1 <- sampleParticipant2HTM(sds=p_sd[j,],
                                 mus=p_mu[j,],
                                 tmpL1 = totL1[,j],
                                 ncoef=ncoef,
                                 nthres=nthres,
                                 nJ_1=nJ_1,
                                 coll= coll)
    s_m <- sampleParticipantdirichlet(mus[j,],nJ_1)
    a_m_o <- sampleParticipantdirichlet(mu_a_o[j,],nJ_1)
    a_m_n <- sampleParticipantdirichlet(mu_a_n[j,],nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["DO_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["DO_r_1"]])[1])){

        DO_tmp <- pnorm(DO_use[J_1==i,j] + as.numeric(r_1[[i]][["DO_r_1"]][k,] %*% t(coll[[1]]$Z_DO[J_1==i,]))) * coll[[1]]$X_DO[J_1==i,1]
        DN_tmp <- pnorm(DN_use[J_1==i,j] + as.numeric(r_1[[i]][["DN_r_1"]][k,] %*% t(coll[[1]]$Z_DN[J_1==i,])))* coll[[1]]$X_DN[J_1==i,1]
        go_tmp <- pnorm(go_use[J_1==i,j] + as.numeric(r_1[[i]][["go_r_1"]][k,] %*% t(coll[[1]]$Z_go[J_1==i,])))




        tmp <- maketheta2HTM_LOP(DO_tmp,
                                 DN_tmp,
                                 go_tmp,
                                 s_m[[i]]$sm_r_1[k,],
                                 a_m_o[[i]]$sm_r_1[k,],
                                 a_m_n[[i]]$sm_r_1[k,],
                                 old = coll[[1]]$X_DO[J_1==i,1],
                                 Yactual = Y[J_1==i])

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp))

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}


# Deviance model prediction funcions --------------------


DPSDTpredictionFull <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b) +
    apply(coll[[2]]$r_2_mu,1,function(x) x[coll[[1]]$J_2]) #+ item effects
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)+
    apply(coll[[2]]$r_2_disc,1,function(x) x[coll[[1]]$J_2])


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(J_1))
  mutmpr1 <- asplit(tmpr_1[,c(1:(((ncoef-nthres)/2) * nJ_1))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  sigmatmpr1 <- asplit(tmpr_1[,c(((((ncoef-nthres)/2) * nJ_1) + 1):((ncoef-nthres)*nJ_1))],1)
  disc_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=length(Y))
  discr1 <- matrix(ncol = lengthsam,nrow=length(Y))

  for (i in c(1:lengthsam)){

    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1]) * coll[[1]]$mu_Z_1)
    discr1[,i] <- rowSums(apply(disc_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$disc_Z_1)

  }

  #prep
  mu_use = mu + mur1
  disc_use = as.matrix(coll[[1]]$X_disc)[,1] *  pnorm(disc + discr1)

  # recollection only applies to old items

  res <- list()
  for (j in 1:lengthsam){

    res[[j]]  <- makethetaDPSDT(mu_tmp=as.numeric(mu_use[,j]),
                                disc_tmp=as.numeric(disc_use[,j]),
                                criteria =matrix(coll[[2]]$Intercept[j,],
                                                 ncol = nthres,byrow=F)[J_1,])

  }
  res
}
H2TMpredictionFull <- function(coll){

  Y <- coll[[1]]$Y

  DO <- as.matrix(coll[[1]]$X_DO) %*% t(coll[[2]]$b_DO)
  DN <- as.matrix(coll[[1]]$X_DN) %*%  t(coll[[2]]$b_DN)
  go <- as.matrix(coll[[1]]$X_go) %*%  t(coll[[2]]$b_go)

  r2DO <-  apply(coll[[2]]$r_2_DO,1,function(x) x[coll[[1]]$J_2])
  r2DN <-  apply(coll[[2]]$r_2_DN,1,function(x) x[coll[[1]]$J_2])
  r2go <-  apply(coll[[2]]$r_2_go,1,function(x) x[coll[[1]]$J_2])

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

  lengthsam <- dim(DO)[2]
  DOr1 <- matrix(ncol = lengthsam,nrow=length(Y))
  DNr1 <- matrix(ncol = lengthsam,nrow=length(Y))
  gor1 <- matrix(ncol = lengthsam,nrow=length(Y))

  for (i in c(1:lengthsam)){

    DOr1[,i] <- rowSums(apply(DO_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_DO)
    DNr1[,i] <- rowSums(apply(DN_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_DN)
    gor1[,i] <- rowSums(apply(go_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_go)
  }

  #prep

  DO_use = apply(pnorm(DO + DOr1 + r2DO),2,function(x) x * coll[[1]]$X_DO[,1])
  DN_use = apply(pnorm(DN + DNr1 + r2DN),2,function(x) x * coll[[1]]$X_DN[,1])
  go_use = apply(pnorm(go + gor1 + r2go),2,function(x) x * coll[[1]]$X_go[,1])




  res <- list()
  for (j in 1:lengthsam){

    s_m_use <- matrix(coll[[2]]$s_m[j,],
                      ncol = (nthres + 1)/2,byrow=F)[J_1,] /
      rowSums(matrix(coll[[2]]$s_m[j,],
                     ncol = (nthres + 1)/2,byrow=F)[J_1,])


    a_m_o_use <- matrix(coll[[2]]$a_m_o[j,],
                        ncol = (nthres + 1)/2,byrow=F)[J_1,] /
      rowSums(matrix(coll[[2]]$a_m_o[j,],
                     ncol = (nthres + 1)/2,byrow=F)[J_1,])

    a_m_n_use <- matrix(coll[[2]]$a_m_n[j,],
                        ncol = (nthres + 1)/2,byrow=F)[J_1,] /
      rowSums(matrix(coll[[2]]$a_m_n[j,],
                     ncol = (nthres + 1)/2,byrow=F)[J_1,])


    res[[j]]  <- maketheta2HTM(DO_tmp=as.numeric(DO_use[,j]),
                               DN_tmp=as.numeric(DN_use[,j]),
                               go = as.numeric(go_use[,j]),
                               s_m = s_m_use,
                               a_m_o = a_m_o_use,
                               a_m_n = a_m_n_use,

                               old = coll[[1]]$X_DO[,1])

  }
  res
}

UVSDTpredictionFull <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b) +
    apply(coll[[2]]$r_2_mu,1,function(x) x[coll[[1]]$J_2]) #+ item effects
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(J_1))
  mutmpr1 <- asplit(tmpr_1[,c(1:(((ncoef-nthres)/2) * nJ_1))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  sigmatmpr1 <- asplit(tmpr_1[,c(((((ncoef-nthres)/2) * nJ_1) + 1):((ncoef-nthres)*nJ_1))],1)
  disc_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=length(Y))
  discr1 <- matrix(ncol = lengthsam,nrow=length(Y))


  for (i in c(1:lengthsam)){



    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    discr1[,i] <- rowSums(apply(disc_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$disc_Z_1)
  }

  #prep
  mu_use = mu + mur1
  disc_use = disc + discr1

  res <- list()
  for (j in 1:lengthsam){

    res[[j]]  <- makethetaUVSDT(mu_tmp=as.numeric(mu_use[,j]),
                                disc_tmp=as.numeric(disc_use[,j]),
                                criteria =matrix(coll[[2]]$Intercept[j,],
                                                 ncol = nthres,byrow=F)[J_1,])

  }
  res
}
EVSDTpredictionFull <- function(coll){

  J_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(J_1))

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b) +
    apply(coll[[2]]$r_2_mu,1,function(x) x[coll[[1]]$J_2]) #+ item effects
  #mu<- mu + muJ2
  # gamma<- vector(mode = "list", length = 6000)
  #
  # for(m in c(1:6000)){
  #   tmp <- as.matrix(rep.row( coll[[2]]$b_gamma[m,],n=nJ_1))
  #   gamma[[m]] <- tmp
  # }
  #


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  # tmpr_1 <- coll[[2]]$r_1


  tmpr_1 <- coll[[2]]$r_1


  mutmpr1 <- asplit(tmpr_1[,c(1:((ncoef-nthres) * nJ_1))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=length(Y))


  #gamtmpr1 <- asplit(tmpr_1[,c((((ncoef-nthres) * nJ_1) + 1):dim(tmpr_1)[[2]])],1)
  #gam_r_1 <- lapply(gamtmpr1,function(x) matrix(x,ncol=nthres))


  #gamma_use <- list()

  for (i in c(1:lengthsam)){

    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1]) * coll[[1]]$mu_Z_1)
    #gamma_use[[i]] <- gamma[[i]] + gam_r_1[[i]]
  }


  #prep
  mu_use = mu + mur1
  # recollection only applies to old items


  # gamma_tmp <- gamma_use[[1]]

  #2*qnorm(head(cumsum(softmax(cbind(gamma_tmp,0)[1,])),nthres))


  # already rejigged it back to criteria in stand generated quantiteis block. Just double-checking here whether it does
  # actually produce the same numbers when puzzling it back togetehr outside of stan (as that is what we are doing for LOP)
  # it does.


  res <- list()
  for (j in 1:lengthsam){

    res[[j]]  <- makethetaUVSDT(mu_tmp=as.numeric(mu_use[,j]),
                                disc_tmp=rep(0,length(coll[[1]]$J_1)),
                                criteria =matrix(coll[[2]]$Intercept[j,],
                                                 ncol = nthres,byrow=F)[J_1,])

  }
  res
}
M0SDTpredictionFull <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b) +
    apply(coll[[2]]$r_2_mu,1,function(x) x[coll[[1]]$J_2]) #+ item effects
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)+
    apply(coll[[2]]$r_2_disc,1,function(x) x[coll[[1]]$J_2])


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(J_1))

  mutmpr1 <- asplit(tmpr_1[,c(1:(((ncoef-nthres)/2) * nJ_1))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  sigmatmpr1 <- asplit(tmpr_1[,c(((((ncoef-nthres)/2) * nJ_1) + 1):((ncoef-nthres)*nJ_1))],1)
  disc_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/2))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=length(Y))
  discr1 <- matrix(ncol = lengthsam,nrow=length(Y))


  for (i in c(1:lengthsam)){



    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    discr1[,i] <- rowSums(apply(disc_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$disc_Z_1)
  }


  mu_use = mu + mur1
  disc_use =  pnorm(disc + discr1)

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
MASDTpredictionFull <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b) +
    apply(coll[[2]]$r_2_mu,1,function(x) x[coll[[1]]$J_2]) #+ item effects
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)+
    apply(coll[[2]]$r_2_disc,1,function(x) x[coll[[1]]$J_2])
  nonattmu <- as.matrix(coll[[1]]$X_nonattmu) %*%  t(coll[[2]]$b_nonattmu)+
    apply(coll[[2]]$r_2_nonattmu,1,function(x) x[coll[[1]]$J_2])

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(J_1))
  mutmpr1 <- asplit(tmpr_1[,c(1:(((ncoef-nthres)/3) * nJ_1))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/3))

  sigmatmpr1 <- asplit(tmpr_1[,c(  ((((ncoef-nthres)/3) * nJ_1) + 1):( (((ncoef-nthres)/3) * nJ_1) * 2))],1)
  disc_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/3))

  nonattmpr1 <- asplit(tmpr_1[,c(  (  ((((ncoef-nthres)/3) * nJ_1) * 2) + 1):( (((ncoef-nthres)/3) * nJ_1) * 3))],1)
  nonatt_r_1 <- lapply(nonattmpr1,function(x) matrix(x,ncol=(ncoef-nthres)/3))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=length(Y))
  discr1 <- matrix(ncol = lengthsam,nrow=length(Y))
  nonattr1 <- matrix(ncol = lengthsam,nrow=length(Y))

  for (i in c(1:lengthsam)){


    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    discr1[,i] <- rowSums(apply(disc_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$disc_Z_1)
    nonattr1[,i] <- rowSums(apply(nonatt_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$nonattemu_Z_1)
  }


  mu_use = mu + mur1
  disc_use =  pnorm(disc + discr1)
  nonatt_use =  pnorm(nonattmu + nonattr1)

  for (i in c(1:lengthsam)){
    disc_use[,i] =  ifelse(as.matrix(coll[[1]]$X_disc)[,1] == 0,1,disc_use[,i])
    nonatt_use[,i] =  ifelse(as.matrix(coll[[1]]$X_nonattmu)[,1] == 0,1,nonatt_use[,i])
  }

  attenmu_use = nonatt_use * mu_use

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


# sample deviations from group-level for LOP ----------

# for all SDT models, sample distributional parameters + criteria from
# MVN(0,Sigma_alpha) where Sigma_alpha is the covariance matrix
# criteria then are transformed from the real-line via softmax to be ordered in normal space
# functions are the same as in LOPThetaInd_gam.R

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
                     disc_r_1 = matrix(r_1[,c((sepncoef + 1):(sepncoef * 2))],ncol=sepncoef),
                     crit_r_1 = r_1[,c((sepncoef*2+1):dim(r_1)[2])])
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

# 2HTM sampling
sampleParticipant2HTM <- function(sds,mus,tmpL1,ncoef,nthres,nJ_1,coll){

  # sds <- p_sd[j,]
  #  mus <- p_mu[j,]
  #  tmpL1 <- totL1[,j]



  L_1 <- tmpL1
  #matrix(tmpL1,nrow=length(sds),byrow=F) %*% t(matrix(tmpL1,nrow=length(sds),byrow=F))

  # turn L_1 cholesky into correlation matrix by L_1 %*% t(L_1)
  # turn sigma and correlation matrix into covariance by S %*% Corr %*% S
  # Make full covmat for participant effects on mu/sigma and criteria for single sample

  ind <- NULL
  for(i in c(1:nJ_1)){

    r_1_use <- mvtnorm::rmvnorm(1000,
                                mean = mus,
                                sigma= outer(sds,sds) * L_1) #equiv:diag(sds) %*% L_1 %*% diag(sds))



    ind[[i]] <- list(DO_r_1 = matrix(r_1_use[,c(1:coll[[1]]$G_DO)],ncol=coll[[1]]$G_DO),
                     DN_r_1 = matrix(r_1_use[,c((coll[[1]]$G_DO + 1):(coll[[1]]$G_DO + coll[[1]]$G_DN))],ncol=coll[[1]]$G_DN),
                     go_r_1 = matrix(r_1_use[,c((coll[[1]]$G_DO + coll[[1]]$G_DN + 1):(coll[[1]]$G_DO + coll[[1]]$G_DN + coll[[1]]$G_go))],
                                     ncol=coll[[1]]$G_go))

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

makethetaM0SDT <- function(mu_tmp,disc_tmp,criteria){
  theta <- matrix(ncol=(dim(criteria)[[2]] + 1),nrow=length(mu_tmp))
  crits <- cbind(-Inf,criteria,Inf)

  for(Y in c(1:(dim(crits)[[2]]-1))){

    theta[,Y] <- disc_tmp * (pnorm(crits[,Y+1] - mu_tmp) - pnorm(crits[,Y] - mu_tmp)) +
      (1 - disc_tmp) * (pnorm(crits[,Y+1]) - pnorm(crits[,Y]))
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

maketheta2HTM <- function(DO_tmp,DN_tmp,go,s_m,a_m_o,a_m_n,old){


  theta <- matrix(ncol=(dim(a_m_o)[[2]])*2,nrow=length(DO_tmp))
  Cold <- dim(a_m_o)[[2]]
  C <- Cold * 2

  for(i in c(1:length(DO_tmp))){

    if(old[i] == 1){

      for(y in c(1:C))
        if(y > Cold){

          theta[i,y] <- DO_tmp[i] * s_m[i,C + 1 - y] + (1-DO_tmp[i]) * go[i] * a_m_o[i,C + 1 - y]
        } else {
          theta[i,y] <- (1-DO_tmp[i]) * (1-go[i]) * a_m_n[i,y]
        }


    } else {

      for(y in c(1:C))
        if(y <= Cold){

          theta[i,y] <- DN_tmp[i] * s_m[i,y] + (1-DN_tmp[i]) * (1-go[i])* a_m_n[i,y]
        } else {
          theta[i,y] <- (1-DN_tmp[i]) * go[i] * a_m_o[i,C + 1 - y]
        }


    }

  }
  theta

}



makethetaUVSDT_LOP <- function(mu_tmp,disc_tmp,criteria,Yactual){

  crits <- c(-Inf,criteria,Inf)

  theta <- pnorm(exp(disc_tmp) * (crits[(Yactual+1)] - mu_tmp)) -
    pnorm(exp(disc_tmp) * (crits[Yactual] - mu_tmp))

}

makethetaDPSDT_LOP <- function(mu_tmp,disc_tmp,criteria,Yactual){

  crits <- c(-Inf,criteria,Inf)


  theta <- (1 - disc_tmp) * (pnorm(crits[Yactual+1] - mu_tmp) -
                               pnorm(crits[Yactual] - mu_tmp))

  theta[Yactual == (length(crits)-1)] <- theta[Yactual == (length(crits)-1)] + disc_tmp[Yactual == (length(crits)-1)]

  return(theta)

}



makethetaM0SDT_LOP <- function(mu_tmp,disc_tmp,criteria,Yactual){

  crits <- c(-Inf,criteria,Inf)

  theta <- disc_tmp * (pnorm(crits[Yactual+1] - mu_tmp) - pnorm(crits[Yactual] - mu_tmp)) +
    (1 - disc_tmp) * (pnorm(crits[Yactual+1]) - pnorm(crits[Yactual]))

  return(theta)

}


makethetaMASDT_LOP <- function(mu_tmp,disc_tmp,att_tmp,criteria,Yactual){

  crits <- c(-Inf,criteria,Inf)

  theta <- disc_tmp * (pnorm(crits[Yactual+1] - mu_tmp) - pnorm(crits[Yactual] - mu_tmp)) +
    (1 - disc_tmp) * (pnorm(crits[Yactual+1]- att_tmp) - pnorm(crits[Yactual] - att_tmp))

  return(theta)

}


maketheta2HTM_LOP <- function(DO_tmp,DN_tmp,go,s_m,a_m_o,a_m_n,old,Yactual){


  theta <- vector(length=length(DO_tmp))
  Cold <- length(a_m_o)
  C <- Cold * 2

  for(i in c(1:length(DO_tmp))){

    if(old[i] == 1){

      #for(y in c(1:C))
      if(Yactual[i] > Cold){

        theta[i] <- DO_tmp[i] * s_m[C + 1 - Yactual[i]] + (1-DO_tmp[i]) * go[i] * a_m_o[C + 1 - Yactual[i]]
      } else {
        theta[i] <- (1-DO_tmp[i]) * (1-go[i]) * a_m_n[Yactual[i]]
      }


    } else {

      # for(y in c(1:C))
      if(Yactual[i] <= Cold){

        theta[i] <- DN_tmp[i] * s_m[Yactual[i]] + (1-DN_tmp[i]) * (1-go[i])* a_m_n[Yactual[i]]
      } else {
        theta[i] <- (1-DN_tmp[i]) * go[i] * a_m_o[C + 1 - Yactual[i]]
      }


    }

  }
  theta

}

# Deviance from probabilities ----------------


getDeviancefromTheta_individual <- function(obs,thetaList,J_1){
  # obs<-tmpdat_pred$Y
  #J_1 <- exp_pred %>% filter(Folds==fold) %>% .$id
  # thetaList <- pred
  srat <-matrix(c(1:length(obs),obs),byrow=F,ncol=2)
  probs <- lapply(thetaList, FUN = `[`, srat)



  probs<-lapply(probs,function(x) replace(x,x==0,1e-8))


  lprobs <- lapply(probs,FUN=log)

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

# Make Deviance from Theta

getDeviancefromTheta <- function(obs,thetaList){

  srat <-matrix(c(1:length(obs),obs),byrow=F,ncol=2)
  probs <- lapply(thetaList, FUN = `[`, srat)

  probs<-lapply(probs,function(x) replace(x,x==0,1e-8))
  -2 * (lapply(lapply(probs,FUN=log),FUN = sum) %>% unlist())

}

# point-wise LL ----------------------


getpwLL <- function(obs,thetaList){

  srat <-matrix(c(1:length(obs),obs),byrow=F,ncol=2)
  probs <- lapply(thetaList, FUN = `[`, srat)


  probs<-lapply(probs,function(x) replace(x,x==0,1e-8))

  lapply(probs,FUN=log)

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


  J_1 <- exp$id # participant number per trial
  conditions <- paste0(exp$condition,"_",exp$isold) # condition per trial (_isold for within-subject/within-block conditions)
  Y <- exp$rating # ratings per trial

  # ratings per trial in matrix of trial by rating category
  nobs <- matrix(0,nrow=dim(theta)[1],ncol=max(Y))

  for(i in seq_along(Y)){
    nobs[i,Y[i]] <- 1
  }

  # item level T1:

  nhat <- theta
  nhat[nhat == 0] <- 1e-8

  npred <- extraDistr::rmnom(dim(theta)[1],1,theta)

  T1itempred <- T1statX2(npred,nhat)
  T1itemobs <- T1statX2(nobs,nhat)

  item_subj <- bind_cols(itemT1preds = rowSums(T1itempred),
                         itemT1obs = rowSums(T1itemobs),
                         id = J_1) %>%
    group_by(id) %>%
    summarize_all(sum)



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


  T1 <- left_join(item_subj,ind_subj) %>% bind_cols(expT1preds=expT1preds,
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


