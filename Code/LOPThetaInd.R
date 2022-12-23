# General information ------------------------

# No item effects, criteria sampled from MVN(mu_cr,sigma_cr) for individuals
# for models EVSDT, Gumbel, UVSDT, DPSDT, DPRMSDT, DPRM2SDT, M0SDT, MASDT, 2HTM
# functions with 'Full' take individual posterior samples for prediction
# functions with 'LOP' take group-level distributions, wjhere necessary, to sample/approximate individual posterior samples


# tmpdat<-tmpdat_pred
# fit<-fit_f
# model<-modelpred
# Top level prediction

# collect standata and posterior samples in a single object for all models -----

prepTheta_Ind <- function(tmpdat,fit, model){

  if (model %in% c("UVSDT","DPSDT")){

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
                       J_2 = tmpdat$J_2)

  }  else if (model %in% c("MASDT")){
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
                       nonattemu_Z_1 = tmpdat$Z_nonattmu)
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
                       #G_nonattmu = tmpdat$G_nonattmu,
                       ncoef=ncoef,
                       mu_Z_1 = tmpdat$Z,
                       sigma_Z_1 = tmpdat$Z_disc)
  }else if (model %in% c("EVSDT","Gumbel")){


    # in early EVSDT fits, we used different convention for
    # Z matrices before switching to a more general model
    # so some of those early fits won't work with this template
    # of recovering the prediction

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
                       #L_1 = apply(fit$draws("L_1")[,,],3,rbind),
                       mu_cr = apply(fit$draws("mu_cr")[,,],3,rbind),
                       sigma_cr = apply(fit$draws("sigma_cr")[,,],3,rbind))

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

# for out-of-sample LOP prediction
predictionLOP_Ind <- function(tmpdat,fit,model){

  coll <- prepTheta_Ind(tmpdat,fit,model)

  if (model == "UVSDT"){

    out <- UVSDTpredictionLOP_Ind(coll)

  } else if (model == "EVSDT"){

    out <- EVSDTpredictionLOP_Ind(coll)

  } else if (model == "Gumbel"){

    out <- GumbelpredictionLOP_Ind(coll)

  } else if (model == "DPSDT"){

    out <- DPSDTpredictionLOP_Ind(coll)
  } else if (model == "DPRM2SDT"){

    out <- DPRM2SDTpredictionLOP_Ind(coll)
  } else if (model == "DPRMSDT"){

    out <- DPRMSDTpredictionLOP_Ind(coll)
  }  else if (model == "MASDT"){

    out <- MASDTpredictionLOP_Ind(coll)
  } else if (model == "M0SDT"){

    out <- M0SDTpredictionLOP_Ind(coll)
  }

  out
}

# for deviance and out-of-sample KFCV prediction
predictionFULL_Ind <- function(tmpdat,fit,model){


  #test = rep(TRUE,tmpdat$N)

  coll <- prepTheta_Ind(tmpdat,fit,model)

  if (model == "UVSDT"){

    out <- UVSDTpredictionFULL_Ind(coll)

  } else if (model == "EVSDT"){

    out <- EVSDTpredictionFULL_Ind(coll)

  } else if (model == "Gumbel"){

    out <- GumbelpredictionFULL_Ind(coll)


  } else if (model == "DPSDT"){

    out <- DPSDTpredictionFull_Ind(coll)
  } else if (model == "DPRM2SDT"){

    out <- DPRM2SDTpredictionFull_Ind(coll)
  } else if (model == "DPRMSDT"){

    out <- DPRMSDTpredictionFull_Ind(coll)
  } else if (model == "2HTM"){

    if (tmpdat$exp %in% c("HUW2015_e1","LBA2019_e1","LBA2019_e2","LM2020_e1")) {
      out <- H2TMneutralpredictionFull_Ind(coll)
    } else {
      out <- H2TMpredictionFull_Ind(coll)
    }
  } else if (model == "M0SDT"){
    out <- M0SDTpredictionFull_Ind(coll)
  } else if (model == "MASDT"){
    out <- MASDTpredictionFull_Ind(coll)
  }


  out
}

# LOP Model prediction functions ---------------
UVSDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)
  #r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  p_sd <- cbind(coll[[2]]$sd_1,coll[[2]]$sigma_cr)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]),
                coll[[2]]$mu_cr)
  L_1 <- coll[[2]]$L_1

  totL1 <- apply(L_1,1,function(x) extendcor(x,ncoef,nthres))

  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu# + r2
  lengthsam <- dim(mu)[2]

  res <- list()
  for (j in 1:lengthsam){

    r_1 <- sampleParticipantUVSDT(p_sd[j,],p_mu[j,],totL1[,j],ncoef,nthres,nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){



        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- disc_use[J_1==i,j] + as.numeric(r_1[[i]][["sigma_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))

        tmp <- makethetaUVSDT_LOP(as.numeric(mu_tmp),
                                  disc_tmp,
                                  as.numeric(r_1[[i]][["crit_r_1"]][k,]))

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

EVSDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  # disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)
  # r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(J_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  #L_1 <- coll[[2]]$L_1


  p_sd <- cbind(coll[[2]]$sd_1,coll[[2]]$sigma_cr)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]),
                coll[[2]]$mu_cr)
  # totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #  totL1[,1:ncoef] <- L_1[,1:ncoef]
  # totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu #+ r2
  lengthsam <-  dim(mu)[2]


  # generate for each participant that is generated out of sample:

  res <- list()
  for (j in 1:lengthsam){

    r_1 <- sampleParticipantEVSDT(p_sd[j,],p_mu[j,],ncoef,nthres,nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      #tmpres<- matrix(nrow=dim(Y[J_1==i,])[1],ncol=dim(r_1[["mu_r_1"]])[1])

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){


        mu_tmp <- mu_use[J_1 == i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        # disc_tmp <- disc_use[,j] + rowSums(r_1[["sigma_r_1"]][k,] * coll[[1]]$mu_Z_1)

        tmp <- makethetaUVSDT_LOP(as.numeric(mu_tmp),
                                  as.numeric(0),
                                  as.numeric(r_1[[i]][["crit_r_1"]][k,]))

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

GumbelpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  # disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)
  # r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  #L_1 <- coll[[2]]$L_1


  p_sd <- cbind(coll[[2]]$sd_1,coll[[2]]$sigma_cr)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]),
                coll[[2]]$mu_cr)
  # totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #  totL1[,1:ncoef] <- L_1[,1:ncoef]
  # totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu #+ r2
  lengthsam <-  dim(mu)[2]


  # generate for each participant that is generated out of sample:

  res <- list()
  for (j in 1:lengthsam){

    r_1 <- sampleParticipantEVSDT(p_sd[j,],p_mu[j,],ncoef,nthres,nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      #tmpres<- matrix(nrow=dim(Y[J_1==i,])[1],ncol=dim(r_1[["mu_r_1"]])[1])

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){


        mu_tmp <- mu_use[J_1 == i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        # disc_tmp <- disc_use[,j] + rowSums(r_1[["sigma_r_1"]][k,] * coll[[1]]$mu_Z_1)

        tmp <- makethetaGumbel_LOP(as.numeric(mu_tmp),
                                   #as.numeric(0),
                                   as.numeric(r_1[[i]][["crit_r_1"]][k,]))

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

DPSDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)
  #r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

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

  #prep without r_1s
  mu_use = mu# + r2
  lengthsam <- dim(mu)[2]

  res <- list()
  for (j in 1:lengthsam){



    r_1 <- sampleParticipantUVSDT(p_sd[j,],p_mu[j,],totL1[,j],ncoef,nthres,nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){


      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){


        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- as.matrix(coll[[1]]$X_disc[J_1==i,])[,1] * pnorm(disc_use[J_1==i,j] + as.numeric(r_1[[i]][["sigma_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,])))

        tmp <- makethetaDPSDT_LOP(as.numeric(mu_tmp),
                                  disc_tmp,
                                  as.numeric(r_1[[i]][["crit_r_1"]][k,]))

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}



DPRMSDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)
  #r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  L_1 <- coll[[2]]$L_1

  p_sd <- cbind(coll[[2]]$sd_1,coll[[2]]$sigma_cr)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]),
                coll[[2]]$mu_cr)

  p_mus <- coll[[2]]$mu_s

  totL1 <- apply(L_1,1,function(x) extendcor(x,ncoef,nthres))

  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu# + r2
  lengthsam <- dim(mu)[2]

  res <- list()
  for (j in 1:lengthsam){


    s_m <- sampleParticipantdirichlet(p_mus[j,],nJ_1)

    r_1 <- sampleParticipantUVSDT(p_sd[j,],p_mu[j,],totL1[,j],ncoef,nthres,nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){


        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- as.matrix(coll[[1]]$X_disc[J_1==i,])[,1] * pnorm(disc_use[J_1==i,j] + as.numeric(r_1[[i]][["sigma_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,])))
        sm_tmp <- s_m[[i]]$sm_r_1[k,]

        tmp <- makethetaDPRMSDT_LOP(mu_tmp = as.numeric(mu_tmp),
                                    disc_tmp = disc_tmp,
                                    criteria = as.numeric(r_1[[i]][["crit_r_1"]][k,]),
                                    s_m = sm_tmp)

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}
DPRM2SDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)
  #r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres
  L_1 <- coll[[2]]$L_1

  p_sd <- cbind(coll[[2]]$sd_1,coll[[2]]$sigma_cr)
  p_mu <- cbind(matrix(0,ncol=dim(coll[[2]]$sd_1)[2],nrow=dim(coll[[2]]$sd_1)[1]),
                coll[[2]]$mu_cr)

  p_mus <- coll[[2]]$mu_s

  totL1 <- apply(L_1,1,function(x) extendcor(x,ncoef,nthres))

  #
  #   totL1 <- rep.row(matrix(diag(dim(p_sd)[2]),nrow=1,byrow=T),dim(p_sd)[1])# Make cholesky-ish for criteria
  #   totL1[,1:ncoef] <- L_1[,1:ncoef]
  #   totL1[,(dim(p_sd)[2]+1):(dim(p_sd)[2] + ncoef)] <- L_1[,(ncoef+1):(dim(L_1)[2])]

  #prep without r_1s
  mu_use = mu# + r2
  lengthsam <- dim(mu)[2]

  res <- list()
  for (j in 1:lengthsam){


    s_m <- sampleParticipantdirichlet(p_mus[j,],nJ_1)

    r_1 <- sampleParticipantUVSDT(p_sd[j,],p_mu[j,],totL1[,j],ncoef,nthres,nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){


        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- as.matrix(coll[[1]]$X_disc[J_1==i,])[,1] * pnorm(disc_use[J_1==i,j] + as.numeric(r_1[[i]][["sigma_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,])))
        sm_tmp <- s_m[[i]]$sm_r_1[k,]

        tmp <- makethetaDPRM2SDT_LOP(mu_tmp = as.numeric(mu_tmp),
                                     disc_tmp = disc_tmp,
                                     criteria = as.numeric(r_1[[i]][["crit_r_1"]][k,]),
                                     s_m = sm_tmp)

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

M0SDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)
  #r2 <-  apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2])

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

  #prep without r_1s
  mu_use = mu# + r2
  lengthsam <- dim(mu)[2]

  res <- list()
  for (j in 1:lengthsam){



    r_1 <- sampleParticipantUVSDT(p_sd[j,],p_mu[j,],totL1[,j],ncoef,nthres,nJ_1)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){


      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){


        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- pnorm(disc_use[J_1==i,j] + as.numeric(r_1[[i]][["sigma_r_1"]][k,] %*% t(coll[[1]]$sigma_Z_1[J_1==i,])))


        for(m in seq_along(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1])){

          disc_tmp[m] <- ifelse(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1][m] == 0,1,disc_tmp[m])
        }



        tmp <- makethetaM0SDT_LOP(as.numeric(mu_tmp),
                                  disc_tmp,
                                  as.numeric(r_1[[i]][["crit_r_1"]][k,]))

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

MASDTpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc_use <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)
  nonattenmu_use <- as.matrix(coll[[1]]$X_nonattmu) %*%  t(coll[[2]]$b_nonattmu)

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

  #prep without r_1s
  mu_use = mu# + r2
  lengthsam <- dim(mu)[2]

  res <- list()
  for (j in 1:lengthsam){


    r_1 <- sampleParticipantMASDT(p_sd[j,],p_mu[j,],totL1[,j],ncoef,nthres,nJ_1,coll)

    forLOP <- vector(length = nJ_1)
    for(i in c(1:nJ_1)){

      tmpres<- vector(length=dim(r_1[[i]][["mu_r_1"]])[1])
      for (k in c(1:dim(r_1[[i]][["mu_r_1"]])[1])){



        mu_tmp <- mu_use[J_1==i,j] + as.numeric(r_1[[i]][["mu_r_1"]][k,] %*% t(coll[[1]]$mu_Z_1[J_1==i,]))
        disc_tmp <- pnorm(disc_use[J_1==i,j] + as.numeric(r_1[[i]][["sigma_r_1"]][k,] %*% t(coll[[1]]$sigma_Z_1[J_1==i,])))
        attenu_tmp <- pnorm(nonattenmu_use[J_1==i,j] + as.numeric(r_1[[i]][["nonatt_r_1"]][k,] %*% t(coll[[1]]$nonattemu_Z_1[J_1==i,])))

        for(m in seq_along(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1])){

          disc_tmp[m] <- ifelse(as.matrix(coll[[1]]$X_disc)[J_1 ==i,1][m] == 0,1,disc_tmp[m])
          attenu_tmp[m] <- ifelse(as.matrix(coll[[1]]$X_nonattmu)[J_1 ==i,1][m] == 0,1,attenu_tmp[m])
        }


        attenmu_use = attenu_tmp * mu_tmp

        tmp <- makethetaMASDT_LOP(mu_tmp = as.numeric(mu_tmp),
                                  disc_tmp = disc_tmp,
                                  attenmu_use = attenmu_use,
                                  criteria= as.numeric(r_1[[i]][["crit_r_1"]][k,]))

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

H2TMpredictionLOP_Ind <- function(coll){

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
  DO_use = DO# + r2
  DN_use = DN
  go_use = go

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
                                 old = coll[[1]]$X_DO[J_1==i,1])

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}
H2TMneutralpredictionLOP_Ind <- function(coll){

  Y <- coll[[1]]$Y
  tJ_1 <- coll[[1]]$J_1
  nJ_1 <- length(unique(tJ_1))

  J_1 <- rep(c(1:nJ_1),each=length(tJ_1)/nJ_1) # do not use original id integers


  DO <- as.matrix(coll[[1]]$X_DO) %*% t(coll[[2]]$b_DO)
  DN <- as.matrix(coll[[1]]$X_DN) %*%  t(coll[[2]]$b_DN)
  go <- as.matrix(coll[[1]]$X_go) %*%  t(coll[[2]]$b_go)
  neut <- as.matrix(coll[[1]]$X_neut) %*%  t(coll[[2]]$b_neut)

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
  DO_use = DO# + r2
  DN_use = DN
  go_use = go
  neut_use = neut

  lengthsam <- dim(DO)[2]

  res <- list()
  for (j in 1:lengthsam){


    r_1 <- sampleParticipant2HTM_neutral(p_sd[j,],p_mu[j,],totL1[,j],ncoef,nthres,nJ_1,coll)
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




        tmp <- maketheta2HTMneutral_LOP(DO_tmp,
                                        DN_tmp,
                                        go_tmp,
                                        s_m[[i]]$sm_r_1[k,],
                                        a_m_o[[i]]$sm_r_1[k,],
                                        a_m_n[[i]]$sm_r_1[k,],
                                        coll[[1]]$X_DO[J_1==i,1])

        tmp[tmp == 0] <- 1e-8
        tmpres[[k]]<-sum(log(tmp) * Y[J_1==i,])

      }

      forLOP[[i]] <- mean(tmpres)


    }

    res[[j]] <- forLOP

  }
  res
}

# Deviance Model prediction functions ---------------------

EVSDTpredictionFULL_Ind <- function(coll){


  Y <- coll[[1]]$Y



  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  #mur2 <- apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2]) #+ item effects
  #disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1

  mu_r_1 <- lapply(asplit(tmpr_1,1),function(x) matrix(x,ncol=ncoef))

  #sigmatmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/2) + 1):dim(tmpr_1)[2])],1)
  #sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=ncoef/2))


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

GumbelpredictionFULL_Ind <- function(coll){


  Y <- coll[[1]]$Y



  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  #mur2 <- apply(coll[[2]]$r_2,1,function(x) x[coll[[1]]$J_2]) #+ item effects
  #disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1

  mu_r_1 <- lapply(asplit(tmpr_1,1),function(x) matrix(x,ncol=ncoef))

  #sigmatmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/2) + 1):dim(tmpr_1)[2])],1)
  #sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=ncoef/2))


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

    res[[j]]  <- makethetaGumbel(mu_tmp=as.numeric(mu_use[,j]),
                                 #disc_tmp=rep(0,length(coll[[1]]$J_1)),
                                 criteria =matrix(coll[[2]]$Intercept[j,],
                                                  ncol = nthres,byrow=F)[J_1,])



  }
  res
}

UVSDTpredictionFULL_Ind <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1

  mutmpr1 <- asplit(tmpr_1[,c(1:(dim(tmpr_1)[2]/2))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=ncoef/2))

  sigmatmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/2) + 1):dim(tmpr_1)[2])],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=ncoef/2))

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

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1

  mutmpr1 <- asplit(tmpr_1[,c(1:(dim(tmpr_1)[2]/2))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=ncoef/2))

  sigmatmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/2) + 1):dim(tmpr_1)[2])],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=ncoef/2))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])


  for (i in c(1:lengthsam)){



    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1[,i] <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
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


DPRMSDTpredictionFull_Ind <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1

  mutmpr1 <- asplit(tmpr_1[,c(1:(dim(tmpr_1)[2]/2))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=ncoef/2))

  sigmatmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/2) + 1):dim(tmpr_1)[2])],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=ncoef/2))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])


  for (i in c(1:lengthsam)){



    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1[,i] <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
  }

  #prep
  mu_use = mu + mur1
  disc_use =  as.matrix(coll[[1]]$X_disc)[,1] * pnorm(disc + sigmar1)



  res <- list()
  for (j in 1:lengthsam){


    # rounding error in stan output -> is defined as simplex in the model
    # rowSums != 1 in the 7th decimal place
    # rejig it here so that multinomial probabilities = 1

    s_m_use <- matrix(coll[[2]]$s_m[j,],
                      ncol = (nthres + 1)/2,byrow=F)[J_1,] /
      rowSums(matrix(coll[[2]]$s_m[j,],
                     ncol = (nthres + 1)/2,byrow=F)[J_1,])


    res[[j]]  <- makethetaDPRMSDT(mu_tmp=as.numeric(mu_use[,j]),
                                  disc_tmp=as.numeric(disc_use[,j]),
                                  criteria =matrix(coll[[2]]$Intercept[j,],
                                                   ncol = nthres,byrow=F)[J_1,],
                                  s_m = s_m_use)

  }
  res
}

DPRM2SDTpredictionFull_Ind <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)


  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1

  mutmpr1 <- asplit(tmpr_1[,c(1:(dim(tmpr_1)[2]/2))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=ncoef/2))

  sigmatmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/2) + 1):dim(tmpr_1)[2])],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=ncoef/2))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])


  for (i in c(1:lengthsam)){



    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1[,i] <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
  }

  #prep
  mu_use = mu + mur1
  disc_use =  as.matrix(coll[[1]]$X_disc)[,1] * pnorm(disc + sigmar1)



  res <- list()
  for (j in 1:lengthsam){


    # rounding error in stan output -> is defined as simplex in the model
    # rowSums != 1 in the 7th decimal place
    # rejig it here so that multinomial probabilities = 1

    s_m_use <- matrix(coll[[2]]$s_m[j,],
                      ncol = 2,byrow=F)[J_1,] /
      rowSums(matrix(coll[[2]]$s_m[j,],
                     ncol = 2,byrow=F)[J_1,])


    res[[j]]  <- makethetaDPRM2SDT(mu_tmp=as.numeric(mu_use[,j]),
                                   disc_tmp=as.numeric(disc_use[,j]),
                                   criteria =matrix(coll[[2]]$Intercept[j,],
                                                    ncol = nthres,byrow=F)[J_1,],
                                   s_m = s_m_use)

  }
  res
}
H2TMpredictionFull_Ind <- function(coll){

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

  lengthsam <- dim(DO)[2]
  DOr1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  DNr1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  gor1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])

  for (i in c(1:lengthsam)){

    DOr1[,i] <- rowSums(apply(DO_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_DO)
    DNr1[,i] <- rowSums(apply(DN_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_DN)
    gor1[,i] <- rowSums(apply(go_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_go)
  }

  #prep

  DO_use = apply(pnorm(DO + DOr1),2,function(x) x * coll[[1]]$X_DO[,1])
  DN_use = apply(pnorm(DN + DNr1),2,function(x) x * coll[[1]]$X_DN[,1])
  go_use = apply(pnorm(go + gor1),2,function(x) x * coll[[1]]$X_go[,1])




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
H2TMneutralpredictionFull_Ind <- function(coll){

  Y <- coll[[1]]$Y

  DO <- as.matrix(coll[[1]]$X_DO) %*% t(coll[[2]]$b_DO)
  DN <- as.matrix(coll[[1]]$X_DN) %*%  t(coll[[2]]$b_DN)
  go <- as.matrix(coll[[1]]$X_go) %*%  t(coll[[2]]$b_go)
  neut <- as.matrix(coll[[1]]$X_neut) %*% t(coll[[2]]$b_neut)

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

  neuttmpr1 <- asplit(tmpr_1[,c((lengthDO + lengthDN + lengthgo + 1) : (lengthDO + lengthDN + lengthgo + lengthneut))],1)
  neut_r_1 <- lapply(neuttmpr1,function(x) matrix(x,ncol=coll[[1]]$G_neut))

  lengthsam <- dim(DO)[2]
  DOr1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  DNr1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  gor1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  neutr1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])

  for (i in c(1:lengthsam)){

    DOr1[,i] <- rowSums(apply(DO_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_DO)
    DNr1[,i] <- rowSums(apply(DN_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_DN)
    gor1[,i] <- rowSums(apply(go_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_go)
    neutr1[,i] <- rowSums(apply(neut_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$Z_neut)
  }

  #prep

  DO_use = apply(pnorm(DO + DOr1),2,function(x) x * coll[[1]]$X_DO[,1])
  DN_use = apply(pnorm(DN + DNr1),2,function(x) x * coll[[1]]$X_DN[,1])
  go_use = apply(pnorm(go + gor1),2,function(x) x * coll[[1]]$X_go[,1])
  neut_use = apply(pnorm(neut + neutr1),2,function(x) x * coll[[1]]$X_neut[,1])



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


    res[[j]]  <- maketheta2HTMneutral(DO_tmp=as.numeric(DO_use[,j]),
                                      DN_tmp=as.numeric(DN_use[,j]),
                                      go = as.numeric(go_use[,j]),
                                      neut = as.numeric(neut_use[,j]),
                                      s_m = s_m_use,
                                      a_m_o = a_m_o_use,
                                      a_m_n = a_m_n_use,

                                      old = coll[[1]]$X_DO[,1])

  }
  res
}

M0SDTpredictionFull_Ind <- function(coll){

  Y <- coll[[1]]$Y

  mu <- as.matrix(coll[[1]]$X) %*% t(coll[[2]]$b)
  disc <- as.matrix(coll[[1]]$X_disc) %*%  t(coll[[2]]$b_disc)

  ncoef <- coll[[1]]$ncoef
  nthres <- coll[[1]]$nthres

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1

  mutmpr1 <- asplit(tmpr_1[,c(1:(dim(tmpr_1)[2]/2))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=ncoef/2))

  disctmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/2) + 1):(dim(tmpr_1)[2]))],1)
  disc_r_1 <- lapply(disctmpr1,function(x) matrix(x,ncol=ncoef/2))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  discr1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])

  for (i in c(1:lengthsam)){

    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1]) * coll[[1]]$mu_Z_1)
    discr1[,i] <- rowSums(apply(disc_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$sigma_Z_1)

  }

  #prep
  mu_use = mu + mur1
  disc_use =  pnorm(disc + discr1)

  # awkward round-about to set lambda_tmp = 1 for new items to enter default-branch for all new items,
  # and split it with some percentage for old items when calculating response probabilities to allow
  # easy extension to setting a non-attended-old-items mu in the secondary branch that won't be accessed
  # by new items

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

  tmpr_1 <- coll[[2]]$r_1
  J_1 <- coll[[1]]$J_1

  mutmpr1 <- asplit(tmpr_1[,c(1:(dim(tmpr_1)[2]/3))],1)
  mu_r_1 <- lapply(mutmpr1,function(x) matrix(x,ncol=ncoef/3))

  sigmatmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/3) + 1):(dim(tmpr_1)[2]/3 * 2))],1)
  sigma_r_1 <- lapply(sigmatmpr1,function(x) matrix(x,ncol=ncoef/3))

  nonattmutmpr1 <- asplit(tmpr_1[,c(((dim(tmpr_1)[2]/3 * 2) + 1):(dim(tmpr_1)[2]))],1)
  nonattmu_r_1 <- lapply(nonattmutmpr1,function(x) matrix(x,ncol=ncoef/3))

  lengthsam <- dim(mu)[2]
  mur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  sigmar1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])
  nonattmur1 <- matrix(ncol = lengthsam,nrow=dim(Y)[1])

  for (i in c(1:lengthsam)){



    mur1[,i] <- rowSums(apply(mu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$mu_Z_1)
    sigmar1[,i] <- rowSums(apply(sigma_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$sigma_Z_1)
    nonattmur1[,i]<- rowSums(apply(nonattmu_r_1[[i]],2,function(x) x[J_1])  * coll[[1]]$nonattemu_Z_1)
  }

  #prep
  mu_use = mu + mur1
  disc_use =  pnorm(disc + sigmar1)
  attenmu_use_f = pnorm(nonattenmu + nonattmur1)

  for(m in seq_along(as.matrix(coll[[1]]$X_disc)[,1])){

    disc_use[m] <- ifelse(as.matrix(coll[[1]]$X_disc)[m,1] == 0, 1, disc_use[m])
    attenmu_use_f[m] <- ifelse(as.matrix(coll[[1]]$X_disc)[m,1] == 0, 1, attenmu_use_f[m])
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

# sample deviations from group-level for LOP ----------

# for all SDT models, sample distributional parameters from
# N(0,sigma_alpha)
# for criteria, sample from N(mu_cr_i, sigma_cr_i) for i in 1:Nthresholds
# practically as a multivariate normal with mu = (rep(0, Npara),mu_cr) and
# covariance matrix where distributional parameters had correlations etc,
# while criteria were uncorrelated to them an to one another

# As criteria had to be ordered, we use rejection sampling, i.e.,
# generate many samples from the multivariate normal and choose the first 1000
# where criteria are ordered.

# sample 2 group-level parameters
sampleParticipantUVSDT <- function(sds,mus,tmpL1,ncoef,nthres,nJ_1){

  # sds <- p_sd[j,]
  #  mus <- p_mu[j,]
  #  tmpL1 <- totL1[j,]
  sepncoef <- ncoef/2


  L_1 <- tmpL1
  #matrix(tmpL1,nrow=length(sds),byrow=F) %*% t(matrix(tmpL1,nrow=length(sds),byrow=F))

  # turn L_1 cholesky into correlation matrix by L_1 %*% t(L_1)
  # turn sigma and correlation matrix into covariance by S %*% Corr %*% S
  # Make full covmat for participant effects on mu/sigma and criteria for single sample

  ind <- NULL
  for(i in c(1:nJ_1)){

    r_1 <- mvtnorm::rmvnorm(15000,
                            mean = mus,
                            sigma= outer(sds,sds) * L_1) #equiv:diag(sds) %*% L_1 %*% diag(sds))

    # only use samples where criteria samples are ordered
    tmpr_1_use <- r_1[apply(r_1, 1, function(x) !is.unsorted(tail(x,nthres))),]
    r_1_use<-tmpr_1_use[1:min(1000,dim(tmpr_1_use)[1]),]

    ind[[i]] <- list(mu_r_1 = matrix(r_1_use[,c(1:sepncoef)],ncol=sepncoef),
                     sigma_r_1 = matrix(r_1_use[,c((sepncoef + 1):ncoef)],ncol=sepncoef),
                     crit_r_1 = r_1_use[,c(c(ncoef+1):dim(r_1_use)[2])])
  }

  ind
}

# sample 3 group-level parameters
sampleParticipantMASDT <- function(sds,mus,tmpL1,ncoef,nthres,nJ_1,coll){

  # sds <- p_sd[j,]
  #  mus <- p_mu[j,]
  #  tmpL1 <- totL1[,j]
  #sepncoef <- ncoef/2


  L_1 <- tmpL1
  #matrix(tmpL1,nrow=length(sds),byrow=F) %*% t(matrix(tmpL1,nrow=length(sds),byrow=F))

  # turn L_1 cholesky into correlation matrix by L_1 %*% t(L_1)
  # turn sigma and correlation matrix into covariance by S %*% Corr %*% S
  # Make full covmat for participant effects on mu/sigma and criteria for single sample

  lengthmu <- dim(coll[[1]]$mu_Z_1)[2]
  lengthdisc <- dim(coll[[1]]$sigma_Z_1)[2]
  lengthnon <- dim(coll[[1]]$nonattemu_Z_1)[2]

  ind <- NULL
  for(i in c(1:nJ_1)){

    r_1 <- mvtnorm::rmvnorm(15000,
                            mean = mus,
                            sigma= outer(sds,sds) * L_1) #equiv:diag(sds) %*% L_1 %*% diag(sds))

    # only use samples where criteria samples are ordered
    tmpr_1_use <- r_1[apply(r_1, 1, function(x) !is.unsorted(tail(x,nthres))),]
    r_1_use<-tmpr_1_use[1:min(1000,dim(tmpr_1_use)[1]),]

    ind[[i]] <- list(mu_r_1 = matrix(r_1_use[,c(1:lengthmu)],ncol=lengthmu),
                     sigma_r_1 = matrix(r_1_use[,c((lengthmu + 1):(lengthmu+lengthdisc))],ncol=lengthdisc),
                     nonatt_r_1 = matrix(r_1_use[,c((lengthmu + lengthdisc+1):(lengthmu+lengthdisc+lengthnon))],ncol=lengthnon),
                     crit_r_1 = r_1_use[,c(c(lengthmu+lengthdisc+lengthnon+1):dim(r_1_use)[2])])
  }

  ind
}

# sample 1 group-level parameter
sampleParticipantEVSDT <- function(sds,mus,ncoef,nthres,nJ_1){


  #nJ_1<-length(unique(tmpdat$J_1))


  #  j<-1
  # sds <- p_sd[j,]
  #  mus <- p_mu[j,]
  # tmpL1 <- totL1[j,]
  sepncoef <- ncoef


  # L_1 <- matrix(tmpL1,nrow=length(sds),byrow=F) %*% t(matrix(tmpL1,nrow=length(sds),byrow=F))

  # turn L_1 cholesky into correlation matrix by L_1 %*% t(L_1)
  # turn sigma and correlation matrix into covariance by S %*% Corr %*% S
  # Make full covmat for participant effects on mu/sigma and criteria for single sample

  ind <- NULL
  for(i in c(1:nJ_1)){
    r_1 <- mvtnorm::rmvnorm(15000,
                            mean = mus,
                            sigma= diag(sds^2))

    # only use samples where criteria samples are ordered
    tmpr_1_use <- r_1[apply(r_1, 1, function(x) !is.unsorted(tail(x,nthres))),]
    r_1_use<-tmpr_1_use[1:min(1000,dim(tmpr_1_use)[1]),]



    ind[[i]] <- list(mu_r_1 = matrix(r_1_use[,c(1:sepncoef)],ncol=sepncoef),
                     #  sigma_r_1 = matrix(r_1_use[,c((sepncoef + 1):ncoef)],ncol=sepncoef),
                     crit_r_1 = r_1_use[,c(c(ncoef+1):dim(r_1_use)[2])])
  }

  ind

}

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
sampleParticipant2HTM_neutral <- function(sds,mus,tmpL1,ncoef,nthres,nJ_1,coll){

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
                                     ncol=coll[[1]]$G_go),
                     neut_r_1 = matrix(r_1_use[,c((coll[[1]]$G_DO + coll[[1]]$G_DN + coll[[1]]$G_go + 1):(coll[[1]]$G_DO +
                                                                                                            coll[[1]]$G_DN +
                                                                                                            coll[[1]]$G_go +
                                                                                                            coll[[1]]$G_neut))],
                                       ncol=coll[[1]]$G_neut))

  }
  ind
}

# sample response mapping parameters from dirichlet
sampleParticipantdirichlet <- function(mus,nJ_1){


  ind <- NULL
  for(i in c(1:nJ_1)){
    sm_sample <- DirichletReg::rdirichlet(n = 1000,alpha=mus)


    ind[[i]] <- list(sm_r_1 =sm_sample)
  }

  ind

}




# Calculation of Theta and likelihoods -------------------------------
# for full fits
makethetaUVSDT <- function(mu_tmp,disc_tmp,criteria){


  theta <- matrix(ncol=(dim(criteria)[[2]] + 1),nrow=length(mu_tmp))
  crits <- cbind(-Inf,criteria,Inf)

  for(Y in c(1:(dim(crits)[[2]]-1))){

    theta[,Y] <- pnorm(exp(disc_tmp) * (crits[,Y+1] - mu_tmp)) -
      pnorm(exp(disc_tmp) * (crits[,Y] - mu_tmp))

  }

  theta

}

makethetaGumbel <- function(mu_tmp,criteria){


  theta <- matrix(ncol=(dim(criteria)[[2]] + 1),nrow=length(mu_tmp))
  crits <- cbind(-Inf,criteria,Inf)

  for(Y in c(1:(dim(crits)[[2]]-1))){

    theta[,Y] <- pgumbelmin(xloc=crits[,Y+1],location=mu_tmp,scale=1) -
      pgumbelmin(xloc=crits[,Y],location=mu_tmp,scale=1)

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
maketheta2HTMneutral <- function(DO_tmp,DN_tmp,go,neut,s_m,a_m_o,a_m_n,old){


  theta <- matrix(ncol=(dim(a_m_o)[[2]])*2 + 1,nrow=length(DO_tmp))

  C <- dim(a_m_o)[[2]] * 2 + 1
  Cold <- ceiling(C/2)

  for(i in c(1:length(DO_tmp))){

    if(old[i] == 1){

      for(y in c(1:C))
        if(y > Cold){

          theta[i,y] <- DO_tmp[i] * s_m[i,C + 1 - y] + (1-DO_tmp[i]) * (1-neut[i])* go[i] * a_m_o[i,C + 1 - y]
        } else if(y == Cold) {

          theta[i,y] <- (1-DO_tmp[i]) * neut[i]
        } else {

          theta[i,y] <- (1-DO_tmp[i]) * (1-neut[i]) * (1-go[i]) * a_m_n[i,y]
        }


    } else {

      for(y in c(1:C))
        if(y < Cold){

          theta[i,y] <- DN_tmp[i] * s_m[i,y] + (1-DN_tmp[i]) * (1-neut[i]) * (1-go[i])* a_m_n[i,y]
        } else if(y == Cold){

          theta[i,y] <- (1-DN_tmp[i]) * neut[i]
        } else {
          theta[i,y] <- (1-DN_tmp[i]) * (1-neut[i]) * go[i] * a_m_o[i,C + 1 - y]
        }


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
makethetaDPRMSDT <- function(mu_tmp,disc_tmp,criteria,s_m){


  theta <- matrix(ncol=(dim(criteria)[[2]] + 1),nrow=length(mu_tmp))
  crits <- cbind(-Inf,criteria,Inf)

  for(Y in c(1:(dim(crits)[[2]]-1))){

    if(Y >  (ceiling((dim(crits)[[2]] -1)/2))){

      theta[,Y] <- disc_tmp * s_m[,dim(crits)[[2]] - Y]  + (1 - disc_tmp) * (pnorm(crits[,Y+1] - mu_tmp) -
                                                                               pnorm(crits[,Y] - mu_tmp))
    } else {
      theta[,Y] <- (1 - disc_tmp) * (pnorm(crits[,Y+1] - mu_tmp) -
                                       pnorm(crits[,Y] - mu_tmp))

    }

  }

  theta

}
makethetaDPRM2SDT <- function(mu_tmp,disc_tmp,criteria,s_m){


  theta <- matrix(ncol=(dim(criteria)[[2]] + 1),nrow=length(mu_tmp))
  crits <- cbind(-Inf,criteria,Inf)

  for(Y in c(1:(dim(crits)[[2]]-1))){

    if(Y == (dim(crits)[[2]] -1)){

      theta[,Y] <- disc_tmp * s_m[,1]  + (1 - disc_tmp) * (pnorm(crits[,Y+1] - mu_tmp) -
                                                             pnorm(crits[,Y] - mu_tmp))
    } else if(Y == (dim(crits)[[2]] -2)){

      theta[,Y] <- disc_tmp * s_m[,2]  + (1 - disc_tmp) * (pnorm(crits[,Y+1] - mu_tmp) -
                                                             pnorm(crits[,Y] - mu_tmp))
    } else {
      theta[,Y] <- (1 - disc_tmp) * (pnorm(crits[,Y+1] - mu_tmp) -
                                       pnorm(crits[,Y] - mu_tmp))

    }

  }

  theta


}

makethetaGumbel_LOP <- function(mu_tmp,criteria){


  theta <- matrix(ncol=(length(criteria) + 1),nrow=length(mu_tmp))
  crits <- c(-Inf,criteria,Inf)

  for(Y in c(1:(length(crits)-1))){

    theta[,Y] <- pgumbelmin(xloc=crits[Y+1],location=mu_tmp,scale=1) -
      pgumbelmin(xloc=crits[Y],location=mu_tmp,scale=1)

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
makethetaDPRMSDT_LOP <- function(mu_tmp,disc_tmp,criteria,s_m){


  theta <- matrix(ncol=(length(criteria) + 1),nrow=length(mu_tmp))
  crits <- c(-Inf,criteria,Inf)

  for(Y in c(1:(length(crits)-1))){

    if(Y >  (ceiling((length(crits) -1)/2))){

      theta[,Y] <- disc_tmp * s_m[length(crits) - Y]  + (1 - disc_tmp) * (pnorm(crits[Y+1] - mu_tmp) -
                                                                            pnorm(crits[Y] - mu_tmp))
    } else {
      theta[,Y] <- (1 - disc_tmp) * (pnorm(crits[Y+1] - mu_tmp) -
                                       pnorm(crits[Y] - mu_tmp))

    }

  }

  theta

}
makethetaDPRM2SDT_LOP <- function(mu_tmp,disc_tmp,criteria,s_m){


  theta <- matrix(ncol=(length(criteria) + 1),nrow=length(mu_tmp))
  crits <- c(-Inf,criteria,Inf)

  for(Y in c(1:(length(crits)-1))){

    if(Y == (length(crits) -1)){

      theta[,Y] <- disc_tmp * s_m[1]  + (1 - disc_tmp) * (pnorm(crits[Y+1] - mu_tmp) -
                                                            pnorm(crits[Y] - mu_tmp))
    } else if(Y == (length(crits) -2)){

      theta[,Y] <- disc_tmp * s_m[2]  + (1 - disc_tmp) * (pnorm(crits[Y+1] - mu_tmp) -
                                                            pnorm(crits[Y] - mu_tmp))
    } else {
      theta[,Y] <- (1 - disc_tmp) * (pnorm(crits[Y+1] - mu_tmp) -
                                       pnorm(crits[Y] - mu_tmp))

    }

  }

  theta


}


makethetaM0SDT_LOP <- function(mu_tmp,disc_tmp,criteria){


  theta <- matrix(ncol=(length(criteria) + 1),nrow=length(mu_tmp))
  crits <- c(-Inf,criteria,Inf)

  for(Y in c(1:(length(crits)-1))){

    theta[,Y] <- disc_tmp * (pnorm(crits[Y+1] - mu_tmp) - pnorm(crits[Y] - mu_tmp)) +
      (1 - disc_tmp) * (pnorm(crits[Y+1]) - pnorm(crits[Y]))
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

maketheta2HTM_LOP <- function(DO_tmp,DN_tmp,go,s_m,a_m_o,a_m_n,old){


  theta <- matrix(ncol=length(a_m_o)*2,nrow=length(DO_tmp))
  Cold <- length(a_m_o)
  C <- Cold * 2

  for(i in c(1:length(DO_tmp))){

    if(old[i] == 1){

      for(y in c(1:C)){
        if(y > Cold){

          theta[i,y] <- DO_tmp[i] * s_m[C + 1 - y] + (1-DO_tmp[i]) * go[i] * a_m_o[C + 1 - y]
        } else {
          theta[i,y] <- (1-DO_tmp[i]) * (1-go[i]) * a_m_n[y]
        }
      }

    } else {

      for(y in c(1:C)){
        if(y <= Cold){

          theta[i,y] <- DN_tmp[i] * s_m[y] + (1-DN_tmp[i]) * (1-go[i])* a_m_n[y]
        } else {
          theta[i,y] <- (1-DN_tmp[i]) * go[i] * a_m_o[C + 1 - y]
        }
      }

    }

  }
  theta

}
maketheta2HTMneutral_LOP <- function(DO_tmp,DN_tmp,go,neut,s_m,a_m_o,a_m_n,old){


  theta <- matrix(ncol=length(a_m_o)*2 + 1,nrow=length(DO_tmp))

  C <- length(a_m_o)*2 + 1
  Cold <- ceiling(C/2)

  for(i in c(1:length(DO_tmp))){

    if(old[i] == 1){

      for(y in c(1:C)){
        if(y > Cold){

          theta[i,y] <- DO_tmp[i] * s_m[C + 1 - y] + (1-DO_tmp[i]) * (1-neut[i]) * go[i] * a_m_o[C + 1 - y]
        } else if(y == Cold) {
          theta[i,y] <- (1-DO_tmp[i]) * neut[i]
        } else {
          theta[i,y] <- (1-DO_tmp[i]) * (1-neut[i]) * (1-go[i]) * a_m_n[y]
        }
      }

    } else {

      for(y in c(1:C)){
        if(y < Cold){

          theta[i,y] <- DN_tmp[i] * s_m[y] + (1-DN_tmp[i]) * (1-neut[i]) * (1-go[i])* a_m_n[y]
        } else if (y == Cold){

          theta[i,y] <- (1-DN_tmp[i]) * neut[i]

        } else {
          theta[i,y] <- (1-DN_tmp[i]) * (1-neut[i]) * go[i] * a_m_o[C + 1 - y]
        }

      }
    }

  }
  theta

}






# Deviance from probabilities ----------------

getDeviancefromTheta_Ind <- function(obs,thetaList){

  thetaList<-lapply(thetaList,function(x) replace(x,x==0,1e-8))
  -2 * (lapply(lapply(lapply(thetaList,FUN=log),FUN=function(x) x * obs), FUN=sum) %>% unlist())

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



pgumbelmin <- function(xloc,location,scale) {

  1- exp(-exp(-(-xloc-location)/scale))

}

extendcor <- function(singleL_1,ncoef,nthres){

  cor_coef <- matrix(singleL_1,nrow=ncoef,byrow=F) %*% t(matrix(singleL_1,nrow=ncoef,byrow=F))

  zeromat <- matrix(0,nthres,ncoef)

  cbind(rbind(cor_coef,zeromat),rbind(t(zeromat),diag(nthres)))

}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
