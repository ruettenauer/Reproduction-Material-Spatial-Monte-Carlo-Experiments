# rm(list = ls())


library(ggmap)
library(RgoogleMaps)
library(foreign)
library(sp)
library(GISTools)
library(sp)
library(spdep)
library(texreg)

library(rgdal)
library(raster)


### Impacs for SLX ###

impacts.SLX <- function(obj, ...) {
  stopifnot(!is.null(attr(obj, "mixedImps")))
  n <- nrow(obj$model)
  k <- obj$qr$rank
  impactsWX(attr(obj, "mixedImps"), n, k, type="SLX")
}


impactsWX <- function(obj, n, k, type="SLX") {
  imps <- lapply(obj, function(x) x[, 1])
  names(imps) <- c("direct", "indirect", "total")
  attr(imps, "bnames") <- rownames(obj[[1]])
  ses <- lapply(obj, function(x) x[, 2])
  names(ses) <- c("direct", "indirect", "total")
  attr(ses, "bnames") <- rownames(obj[[1]])
  res <- list(impacts=imps, se=ses)
  attr(res, "n") <- n
  attr(res, "k") <- k
  attr(res, "type") <- type
  attr(res, "method") <- "estimable"
  attr(res, "bnames") <- rownames(obj[[1]])
  class(res) <- "WXImpact"
  res
}


print.WXImpact <- function(x, ...) {
  mat <- lagImpactMat(x$impacts)
  cat("Impact measures (", attr(x, "type"), ", ",
      attr(x, "method"), "):\n", sep="")
  print(mat, ...)
  cat("\n")
  invisible(x)
}

print.summary.WXImpact <- function(x, ...) {
  mat <- x$mat
  cat("Impact measures (", attr(x, "type"), ", ",
      attr(x, "method"), "):\n", sep="")
  print(mat, ...)
  cat("========================================================\n")
  mat <- x$semat
  cat("Standard errors:\n", sep="")
  print(mat, ...)
  cat("========================================================\n")
  cat("Z-values:\n")
  mat <- x$zmat
  rownames(mat) <- attr(x, "bnames")
  print(mat, ...)
  cat("\np-values:\n")
  xx <- apply(x$pzmat, 2, format.pval)
  # 100928 Eelke Folmer
  if (length(attr(x, "bnames")) == 1L) {
    xx <- matrix(xx, ncol=3)
    colnames(xx) <- c("Direct", "Indirect", "Total")
  }
  rownames(xx) <- attr(x, "bnames")
  print(xx, ..., quote=FALSE)
  cat("\n")
  invisible(x)
}

summary.WXImpact <- function(object, ..., adjust_k=FALSE) {
  stopifnot(is.logical(adjust_k))
  stopifnot(length(adjust_k) == 1L)
  object$mat <- lagImpactMat(object$impacts)
  object$semat <- lagImpactMat(object$se)
  if (adjust_k) {
    object$semat <- object$semat * (attr(object, "n")/
                                      (attr(object, "n") - attr(object, "k")))
    attr(object, "method") <- paste(attr(object, "method"),
                                    ", n-k", sep="")
  }
  object$zmat <- object$mat/object$semat
  object$pzmat <- 2*(1-pnorm(abs(object$zmat)))
  class(object) <- c("summary.WXImpact", class(object))
  object
}



#### Start Code ####

# setwd("C:/Users/tor44za/Uni W�rzburg/Forschung/Monte Carlo Spatial/02_Daten")
# 
# 
# ### Parameters
# 
# N<-900
# 
# R<-30
# 
# seed=123579
# 
# # Autocorrelation in y
# ay<-0.2
# 
# # Betas
# bx1<-1.2
# bwx1<-0.5
# 
# bx2<-1.2
# bwx2<-0.5
# 
# bu<-1
# bwu<-0
# 
# # Omitted variable bias
# omv1<-0
# omv2<-0
# # Autocor of omv
# womv1<-0
# womv2<-0


spsim<-function(N=900, R=30, seed=123579, tol=1e-10,
                ay=0, ax1=0, ax2=0, au=0,
                bx1=0.4, bwx1=0, bx2=0.6, bwx2=0,
                bu=1, bwu=0, sdu=0.5,
                omv1=0, womv1=0, omv2=0, womv2=0){
  
  ### Save parameters
  parameters<-c(ay=ay, ax1=ax1, ax2=ax2, au=au, bx1=bx1, bwx1=bwx1, bx2=bx2, bwx2=bwx2,
                bu=bu, bwu=bwu, omv1=omv1, womv1=womv1, 
                N=N, R=R, seed=seed)
  
  
  
  ###################################
  ### Set up raster and listw obj ###
  ###################################
  
  n2<-sqrt(N)
  
  df.sp<-raster(nrows=n2, ncols=n2)
  #df.sp <- extent(df.sp)
  df.sp <- as(df.sp, "SpatialPolygonsDataFrame") 
  
  ### Set up IDs
  df.sp$id<-1:nrow(df.sp)
  spChFIDs(df.sp)<-df.sp$id
  
  
  # ### Create spatial listw object (Queens Nbs)
  # df.nb<-poly2nb(df.sp, queen=T)
  # df.listw<-nb2listw(df.nb, style="W")
  # trW<-trW(W=as(df.listw, "CsparseMatrix"), type="mult")
  
  ### Creat idw listw object with 100 nbs
  coords <- coordinates(df.sp)
  df.nbk <- knn2nb(knearneigh(coords, k = 100))
  dsts <- nbdists(df.nbk, coords)
  idw <- lapply(dsts, function(x) 1/((x/10))) 
  df.listw <- nb2listw(df.nbk, glist = idw, style = "minmax")
  #trW<-trW(W=as(df.listw, "CsparseMatrix"), type="mult")
  

  df.listw$style<-"W" # Fake row-standardization to avoid lagged intercept

  
  ###########################
  ### Compute real values ###
  ###########################
  
  ### Contruct Multiplyer mat
  W<-listw2mat(df.listw)
  
  M<-solve(diag(nrow(df.sp))-ay*W)
  
  Mx1<-solve(diag(nrow(df.sp))-ax1*W)
  Mx2<-solve(diag(nrow(df.sp))-ax2*W)
  
  Mu<-solve(diag(nrow(df.sp))-au*W)
  

  ### Calculate Impacts
  spMat1<- (M*bx1 + M %*% (W*bwx1))
  deff1<-mean(diag(spMat1))
  diag(spMat1)<-NA
  speff1<-mean(rowSums(spMat1, na.rm=T))
  
  spMat2<- (M*bx2 + M %*% (W*bwx2))
  deff2<-mean(diag(spMat2))
  diag(spMat2)<-NA
  speff2<-mean(rowSums(spMat2, na.rm=T))
  
  ### Save values
  theta<-matrix(NA, ncol=4, nrow=4)
  colnames(theta)<-c("", "Coef", "Direct", "Indirect")
  theta[,1]<-c("rho", "lambda", "beta1", "beta2")
  theta[,2]<-c(ay, au, bx1, bx2)
  theta[3:4,3]<-c(deff1, deff2)
  theta[3:4,4]<-c(speff1, speff2)
  
  
  ##############################
  ### Set up result matrices ###
  ##############################
  
  b.coef<-data.frame(matrix(NA, nrow=R, ncol=6))
  colnames(b.coef)<-c("rho", "lambda", "beta1", "beta2", "Wbeta1", "Wbeta2")
  sd.coef<-data.frame(matrix(NA, nrow=R, ncol=6))
  colnames(sd.coef)<-c("rho", "lambda", "beta1", "beta2", "Wbeta1", "Wbeta2")
  
  b.direct<-data.frame(matrix(NA, nrow=R, ncol=2))
  colnames(b.direct)<-c("beta1", "beta2")
  sd.direct<-data.frame(matrix(NA, nrow=R, ncol=2))
  colnames(sd.direct)<-c("beta1", "beta2")
  
  b.indirect<-data.frame(matrix(NA, nrow=R, ncol=2))
  colnames(b.indirect)<-c("beta1", "beta2")
  sd.indirect<-data.frame(matrix(NA, nrow=R, ncol=2))
  colnames(sd.indirect)<-c("beta1", "beta2")
  
  lm.res<-list(b.coef=b.coef, sd.coef=sd.coef, b.direct=b.direct, sd.direct=sd.direct, 
               b.indirect=b.indirect, sd.indirect=sd.indirect)
  slx.res<-list(b.coef=b.coef, sd.coef=sd.coef, b.direct=b.direct, sd.direct=sd.direct, 
                b.indirect=b.indirect, sd.indirect=sd.indirect)
  sar.res<-list(b.coef=b.coef, sd.coef=sd.coef, b.direct=b.direct, sd.direct=sd.direct, 
                b.indirect=b.indirect, sd.indirect=sd.indirect)
  sem.res<-list(b.coef=b.coef, sd.coef=sd.coef, b.direct=b.direct, sd.direct=sd.direct, 
                b.indirect=b.indirect, sd.indirect=sd.indirect)
  sac.res<-list(b.coef=b.coef, sd.coef=sd.coef, b.direct=b.direct, sd.direct=sd.direct, 
                b.indirect=b.indirect, sd.indirect=sd.indirect)
  sdm.res<-list(b.coef=b.coef, sd.coef=sd.coef, b.direct=b.direct, sd.direct=sd.direct, 
                b.indirect=b.indirect, sd.indirect=sd.indirect)
  sdem.res<-list(b.coef=b.coef, sd.coef=sd.coef, b.direct=b.direct, sd.direct=sd.direct, 
                 b.indirect=b.indirect, sd.indirect=sd.indirect)
  
  
  
  ########################
  ### Start Simulation ###
  ########################
  
  set.seed(seed)
  #maxint<-.Machine$integer.max
  maxint<-2147483647
  seeds<-runif(R, min=-maxint, max=maxint)
  
  # print(seeds)
  
  for(i in 1:R){
    
    set.seed(seeds[i])
    
    ### Create an X var
    df.sp$X1<-Mx1 %*% rnorm(nrow(df.sp), mean=0, sd=1)
                       
    df.sp$W_X1<-lag.listw(df.listw, df.sp$X1)
    
    ### Create an X2 var
    df.sp$X2<-Mx2 %*% rnorm(nrow(df.sp), mean=0, sd=1)
                       
    df.sp$W_X2<-lag.listw(df.listw, df.sp$X2)
  
    
    ### Create the disturbances
    df.sp$u<-Mu %*% (rnorm(nrow(df.sp), mean=0, sd=sdu)
                     + omv1*df.sp$X1 
                     + omv2*df.sp$X2)
    
    df.sp$W_u<-lag.listw(df.listw, df.sp$u)
    

    ### Create Y var 
    
    df.sp$Y<-M %*% (  df.sp$X1*bx1 +
                      df.sp$W_X1*bwx1 +
                      df.sp$X2*bx2 +
                      df.sp$W_X2*bwx2 +
                      df.sp$u*bu +
                      df.sp$W_u*bwu
    )
    
    
    
    
    
    
    ################
    ### Analysis ###
    ################
    lm.mod<-NA
    slx.mod<-NA
    sar.mod<-NA
    sem.mod<-NA
    sac.mod<-NA
    sdm.mod<-NA
    sdem.mod<-NA
    
    ### Linear OLS
    lm.mod<-lm(Y ~ X1 + X2, data=data.frame(df.sp) )
    #summary(lm.mod)
    
    ### Spatial SLX model
    slx.mod<-lmSLX(Y ~ X1 + X2, data=data.frame(df.sp), 
                   listw=df.listw)
    #summary(slx.mod)
    
    ### SAR
    sar.mod<-lagsarlm(Y ~ X1 + X2, data=data.frame(df.sp), 
                      listw=df.listw, type="lag", tol.solve=tol)
    #summary(sar.mod)
    
    ### SEM
    sem.mod<-errorsarlm(Y ~ X1 + X2, data=data.frame(df.sp), 
                        listw=df.listw, etype="error", tol.solve=tol)
    #summary(sem.mod)
    
    ### SAC
    sac.mod<-sacsarlm(Y ~ X1 + X2, data=data.frame(df.sp), 
                      listw=df.listw, type="sac", tol.solve=tol)
    
    #summary(sac.mod)
    
    ### SDM
    sdm.mod<-lagsarlm(Y ~ X1 + X2, data=data.frame(df.sp), 
                      listw=df.listw, type="mixed", tol.solve=tol)
    
    #summary(sdm.mod)
    
    ### SDEM
    sdem.mod<-errorsarlm(Y ~ X1 + X2, data=data.frame(df.sp), 
                         listw=df.listw, etype="emixed", tol.solve=tol)
    #summary(sdem.mod)
    
    
    
    ### Calculate Impacts ###
    sar.imp<-NA
    sdm.imp<-NA
    sac.imp<-NA
    # sar.imp<-summary(impacts(sar.mod, tr=trW, R=300))
    # sdm.imp<-summary(impacts(sdm.mod, tr=trW, R=300))
    # sac.imp<-summary(impacts(sac.mod, tr=trW, R=300))
    
    sar.imp<-impacts(sar.mod, listw=df.listw)
    sdm.imp<-impacts(sdm.mod, listw=df.listw)
    sac.imp<-impacts(sac.mod, listw=df.listw)
    
    # slx.imp<-summary(impacts.SLX(slx.mod, tr=trW, R=300))
    
    # sar.imp<-impacts(sar.mod, listw=df.listw)
    # sdm.imp<-impacts(sdm.mod, listw=df.listw)
    # sac.imp<-impacts(sac.mod, listw=df.listw)
    # slx.imp<-impacts.SLX(slx.mod, listw=df.listw)
    
    
    ### Paste into results matrix ###
    
    ### Coefs
    lm.res$b.coef[i,3:4]<-lm.mod$coefficients[-1]
    lm.res$sd.coef[i,3:4]<-sqrt(diag(vcov(lm.mod)))[-1]
    
    slx.res$b.coef[i,3:6]<-slx.mod$coefficients[-1]
    slx.res$sd.coef[i,3:6]<-sqrt(diag(vcov(slx.mod)))[-1]
    
    sar.res$b.coef[i,3:4]<-sar.mod$coefficients[-1]
    sar.res$sd.coef[i,3:4]<-sar.mod$rest.se[-1]
    sar.res$b.coef[i,1]<-sar.mod$rho
    sar.res$sd.coef[i,1]<-sar.mod$rho.se
    
    sem.res$b.coef[i,3:4]<-sem.mod$coefficients[-1]
    sem.res$sd.coef[i,3:4]<-sem.mod$rest.se[-1]
    sem.res$b.coef[i,2]<-sem.mod$lambda
    sem.res$sd.coef[i,2]<-sem.mod$lambda.se
    
    sac.res$b.coef[i,3:4]<-sac.mod$coefficients[-1]
    sac.res$sd.coef[i,3:4]<-sac.mod$rest.se[-1]
    sac.res$b.coef[i,2]<-sac.mod$lambda
    sac.res$sd.coef[i,2]<-sac.mod$lambda.se
    sac.res$b.coef[i,1]<-sac.mod$rho
    sac.res$sd.coef[i,1]<-sac.mod$rho.se
    
    sdm.res$b.coef[i,3:6]<-sdm.mod$coefficients[-1]
    sdm.res$sd.coef[i,3:6]<-sdm.mod$rest.se[-1]
    sdm.res$b.coef[i,1]<-sdm.mod$rho
    sdm.res$sd.coef[i,1]<-sdm.mod$rho.se
    
    sdem.res$b.coef[i,3:6]<-sdem.mod$coefficients[-1]
    sdem.res$sd.coef[i,3:6]<-sdem.mod$rest.se[-1]
    sdem.res$b.coef[i,2]<-sdem.mod$lambda
    sdem.res$sd.coef[i,2]<-sdem.mod$lambda.se
    
    ### Direct Impacts
    
    lm.res$b.direct[i,]<-lm.mod$coefficients[-1]
    lm.res$sd.direct[i,]<-sqrt(diag(vcov(lm.mod)))[-1]
    
    slx.res$b.direct[i,]<-slx.mod$coefficients[2:3]
    slx.res$sd.direct[i,]<-sqrt(diag(vcov(slx.mod)))[2:3]
    
    sar.res$b.direct[i,]<-sar.imp[[1]]
    sar.res$sd.direct[i,]<-NA
    
    sem.res$b.direct[i,]<-sem.mod$coefficients[-1]
    sem.res$sd.direct[i,]<-sem.mod$rest.se[-1]
    
    sac.res$b.direct[i,]<-sac.imp[[1]]
    sac.res$sd.direct[i,]<-NA
    
    sdm.res$b.direct[i,]<-sdm.imp[[1]]
    sdm.res$sd.direct[i,]<-NA
    
    sdem.res$b.direct[i,]<-sdem.mod$coefficients[2:3]
    sdem.res$sd.direct[i,]<-sdem.mod$rest.se[2:3]
    
    ### Direct Impacts
    
    slx.res$b.indirect[i,]<-slx.mod$coefficients[4:5]
    slx.res$sd.indirect[i,]<-sqrt(diag(vcov(slx.mod)))[4:5]
    
    sar.res$b.indirect[i,]<-sar.imp[[2]]
    sar.res$sd.indirect[i,]<-NA
    
    sac.res$b.indirect[i,]<-sac.imp[[2]]
    sac.res$sd.indirect[i,]<-NA
    
    sdm.res$b.indirect[i,]<-sdm.imp[[2]]
    sdm.res$sd.indirect[i,]<-NA
    
    sdem.res$b.indirect[i,]<-sdem.mod$coefficients[4:5]
    sdem.res$sd.indirect[i,]<-sdem.mod$rest.se[4:5]
    
    
    # print(seeds[i])
    if(i==1){
      cat(paste("Simulations completed (by 10):"," "))
    } 
    if(i %% 50==0){
      cat(paste(" ", i," "))
    }else if(i %% 10==0){
      cat(paste("+"))
    }
    if(i==R){
      cat("\n")
    }
    
  }
  
  
  ### Create output element ###

  result<-list(parameters = parameters,
               theta      = theta,
               seeds      = seeds,
               lm.res     = lm.res,
               slx.res    = slx.res,
               sar.res    = sar.res,
               sem.res    = sem.res,
               sac.res    = sac.res,
               sdm.res    = sdm.res,
               sdem.res   = sdem.res)
  
}







######################
#### Summary Sims ####
######################

mse<-function(sim=NULL, true=NULL){
  se<-(sim-true)^2
  mse<-mean(se)
  return(mse)
}

rmse<-function(sim=NULL, true=NULL){
  se<-(sim-true)^2
  mse<-mean(se)
  rmse<-sqrt(mse)
  return(rmse)
}

fn <- function(x) setNames(data.frame(.=paste("", rownames(x)), x, check.names=F, row.names=NULL),c(" ", colnames(x)))



summary.sim<-function(sim=NULL){
  
  ### Set up results matrix
  mat<-matrix(nrow=21, ncol=2)
  colnames(mat)<-c("Bias", "RMSE")
  rownames(mat)<-c("OLS", "x1", "x2",
                   "SLX", "x1", "x2",
                   "SAR", "x1", "x2",
                   "SEM", "x1", "x2",
                   "SAC", "x1", "x2",
                   "SDM", "x1", "x2",
                   "SDEM", "x1", "x2")
  
  mat_d<-mat
  mat_sp<-mat
  
  mat_sd<-mat
  colnames(mat_sd)<-c("Direct", "Indirect")
  
  ### Get true parameters
  dx1<-as.numeric(sim$theta[3,3])
  dx2<-as.numeric(sim$theta[4,3])
  
  spx1<-as.numeric(sim$theta[3,4])
  spx2<-as.numeric(sim$theta[4,4])
  
  ### Direct values
  
  # LM
  mat_d[2,1]<-mean(sim$lm.res$b.direct[,1])-dx1
  mat_d[2,2]<-rmse(sim$lm.res$b.direct[,1], dx1)
  mat_d[3,1]<-mean(sim$lm.res$b.direct[,2])-dx2
  mat_d[3,2]<-rmse(sim$lm.res$b.direct[,2], dx2)
  
  # SLX
  mat_d[5,1]<-mean(sim$slx.res$b.direct[,1])-dx1
  mat_d[5,2]<-rmse(sim$slx.res$b.direct[,1], dx1)
  mat_d[6,1]<-mean(sim$slx.res$b.direct[,2])-dx2
  mat_d[6,2]<-rmse(sim$slx.res$b.direct[,2], dx2)
  
  # SAR
  mat_d[8,1]<-mean(sim$sar.res$b.direct[,1])-dx1
  mat_d[8,2]<-rmse(sim$sar.res$b.direct[,1], dx1)
  mat_d[9,1]<-mean(sim$sar.res$b.direct[,2])-dx2
  mat_d[9,2]<-rmse(sim$sar.res$b.direct[,2], dx2)
  
  # SEM
  mat_d[11,1]<-mean(sim$sem.res$b.direct[,1])-dx1
  mat_d[11,2]<-rmse(sim$sem.res$b.direct[,1], dx1)
  mat_d[12,1]<-mean(sim$sem.res$b.direct[,2])-dx2
  mat_d[12,2]<-rmse(sim$sem.res$b.direct[,2], dx2)
  
  # SAC
  mat_d[14,1]<-mean(sim$sac.res$b.direct[,1])-dx1
  mat_d[14,2]<-rmse(sim$sac.res$b.direct[,1], dx1)
  mat_d[15,1]<-mean(sim$sac.res$b.direct[,2])-dx2
  mat_d[15,2]<-rmse(sim$sac.res$b.direct[,2], dx2)
  
  # SDM
  mat_d[17,1]<-mean(sim$sdm.res$b.direct[,1])-dx1
  mat_d[17,2]<-rmse(sim$sdm.res$b.direct[,1], dx1)
  mat_d[18,1]<-mean(sim$sdm.res$b.direct[,2])-dx2
  mat_d[18,2]<-rmse(sim$sdm.res$b.direct[,2], dx2)
  
  # SDEM
  mat_d[20,1]<-mean(sim$sdem.res$b.direct[,1])-dx1
  mat_d[20,2]<-rmse(sim$sdem.res$b.direct[,1], dx1)
  mat_d[21,1]<-mean(sim$sdem.res$b.direct[,2])-dx2
  mat_d[21,2]<-rmse(sim$sdem.res$b.direct[,2], dx2)
  
  
  ### Indirect values
  
  # LM
  mat_sp[2,1]<-mean(sim$lm.res$b.indirect[,1])-spx1
  mat_sp[2,2]<-rmse(sim$lm.res$b.indirect[,1], spx1)
  mat_sp[3,1]<-mean(sim$lm.res$b.indirect[,2])-spx2
  mat_sp[3,2]<-rmse(sim$lm.res$b.indirect[,2], spx2)
  
  # SLX
  mat_sp[5,1]<-mean(sim$slx.res$b.indirect[,1])-spx1
  mat_sp[5,2]<-rmse(sim$slx.res$b.indirect[,1], spx1)
  mat_sp[6,1]<-mean(sim$slx.res$b.indirect[,2])-spx2
  mat_sp[6,2]<-rmse(sim$slx.res$b.indirect[,2], spx2)
  
  # SAR
  mat_sp[8,1]<-mean(sim$sar.res$b.indirect[,1])-spx1
  mat_sp[8,2]<-rmse(sim$sar.res$b.indirect[,1], spx1)
  mat_sp[9,1]<-mean(sim$sar.res$b.indirect[,2])-spx2
  mat_sp[9,2]<-rmse(sim$sar.res$b.indirect[,2], spx2)
  
  # SEM
  mat_sp[11,1]<-mean(sim$sem.res$b.indirect[,1])-spx1
  mat_sp[11,2]<-rmse(sim$sem.res$b.indirect[,1], spx1)
  mat_sp[12,1]<-mean(sim$sem.res$b.indirect[,2])-spx2
  mat_sp[12,2]<-rmse(sim$sem.res$b.indirect[,2], spx2)
  
  # SAC
  mat_sp[14,1]<-mean(sim$sac.res$b.indirect[,1])-spx1
  mat_sp[14,2]<-rmse(sim$sac.res$b.indirect[,1], spx1)
  mat_sp[15,1]<-mean(sim$sac.res$b.indirect[,2])-spx2
  mat_sp[15,2]<-rmse(sim$sac.res$b.indirect[,2], spx2)
  
  # SDM
  mat_sp[17,1]<-mean(sim$sdm.res$b.indirect[,1])-spx1
  mat_sp[17,2]<-rmse(sim$sdm.res$b.indirect[,1], spx1)
  mat_sp[18,1]<-mean(sim$sdm.res$b.indirect[,2])-spx2
  mat_sp[18,2]<-rmse(sim$sdm.res$b.indirect[,2], spx2)
  
  # SDEM
  mat_sp[20,1]<-mean(sim$sdem.res$b.indirect[,1])-spx1
  mat_sp[20,2]<-rmse(sim$sdem.res$b.indirect[,1], spx1)
  mat_sp[21,1]<-mean(sim$sdem.res$b.indirect[,2])-spx2
  mat_sp[21,2]<-rmse(sim$sdem.res$b.indirect[,2], spx2)
  
  
  #### Standard deviations ####
  
  # LM
  mat_sd[2,1]<-sd(sim$lm.res$b.direct[,1])
  mat_sd[2,2]<-sd(sim$lm.res$b.indirect[,1])
  mat_sd[3,1]<-sd(sim$lm.res$b.direct[,2])
  mat_sd[3,2]<-sd(sim$lm.res$b.indirect[,2])
  
  # SLX
  mat_sd[5,1]<-sd(sim$slx.res$b.direct[,1])
  mat_sd[5,2]<-sd(sim$slx.res$b.indirect[,1])
  mat_sd[6,1]<-sd(sim$slx.res$b.direct[,2])
  mat_sd[6,2]<-sd(sim$slx.res$b.indirect[,2])
  
  # SAR
  mat_sd[8,1]<-sd(sim$sar.res$b.direct[,1])
  mat_sd[8,2]<-sd(sim$sar.res$b.indirect[,1])
  mat_sd[9,1]<-sd(sim$sar.res$b.direct[,2])
  mat_sd[9,2]<-sd(sim$sar.res$b.indirect[,2])
  
  # SEM
  mat_sd[11,1]<-sd(sim$sem.res$b.direct[,1])
  mat_sd[11,2]<-sd(sim$sem.res$b.indirect[,1])
  mat_sd[12,1]<-sd(sim$sem.res$b.direct[,2])
  mat_sd[12,2]<-sd(sim$sem.res$b.indirect[,2])
  
  # SAC
  mat_sd[14,1]<-sd(sim$sac.res$b.direct[,1])
  mat_sd[14,2]<-sd(sim$sac.res$b.indirect[,1])
  mat_sd[15,1]<-sd(sim$sac.res$b.direct[,2])
  mat_sd[15,2]<-sd(sim$sac.res$b.indirect[,2])
  
  # SDM
  mat_sd[17,1]<-sd(sim$sdm.res$b.direct[,1])
  mat_sd[17,2]<-sd(sim$sdm.res$b.indirect[,1])
  mat_sd[18,1]<-sd(sim$sdm.res$b.direct[,2])
  mat_sd[18,2]<-sd(sim$sdm.res$b.indirect[,2])
  
  # SDEM
  mat_sd[20,1]<-sd(sim$sdem.res$b.direct[,1])
  mat_sd[20,2]<-sd(sim$sdem.res$b.indirect[,1])
  mat_sd[21,1]<-sd(sim$sdem.res$b.direct[,2])
  mat_sd[21,2]<-sd(sim$sdem.res$b.indirect[,2])
  
  
  
  
  
  res<-list(mat_d=mat_d, mat_sp=mat_sp, mat_sd=mat_sd,
            theta=sim$theta, parameters=sim$parameters)
  
  class(res) <- c("summary.sim")
  
  res
}




covr.sim<-function(sim=NULL, confint=0.95){
  
  ### Set up results matrix
  mat<-matrix(nrow=21, ncol=1)
  colnames(mat)<-c("Cov Rate")
  rownames(mat)<-c("OLS", "x1", "x2",
                   "SLX", "x1", "x2",
                   "SAR", "x1", "x2",
                   "SEM", "x1", "x2",
                   "SAC", "x1", "x2",
                   "SDM", "x1", "x2",
                   "SDEM", "x1", "x2")
  
  mat_d<-mat
  mat_sp<-mat
  
  interval <- -qnorm((1-confint)/2)
  
  
  ### Get true parameters
  dx1<-as.numeric(sim$theta[3,3])
  dx2<-as.numeric(sim$theta[4,3])
  
  spx1<-as.numeric(sim$theta[3,4])
  spx2<-as.numeric(sim$theta[4,4])
  
  N<-nrow(sim$lm.res$b.coef)
  
  
  ### Calculate confidence intervals
  
  ci.df<-matrix(NA, ncol=8, nrow=N)
  colnames(ci.df)<-c("lb_dir1", "ub_dir1", "lb_dir2", "ub_dir2",
                     "lb_ind1", "ub_ind1", "lb_ind2", "ub_ind2")
  
  # LM
  ci_lm.df<-data.frame(ci.df)
  ci_lm.df$lb_dir1<-sim$lm.res$b.direct[,1]-(sim$lm.res$sd.direct[,1]*interval)
  ci_lm.df$ub_dir1<-sim$lm.res$b.direct[,1]+(sim$lm.res$sd.direct[,1]*interval)
  ci_lm.df$lb_dir2<-sim$lm.res$b.direct[,2]-(sim$lm.res$sd.direct[,2]*interval)
  ci_lm.df$ub_dir2<-sim$lm.res$b.direct[,2]+(sim$lm.res$sd.direct[,2]*interval)
  
  ci_lm.df$lb_ind1<-sim$lm.res$b.indirect[,1]-(sim$lm.res$sd.indirect[,1]*interval)
  ci_lm.df$ub_ind1<-sim$lm.res$b.indirect[,1]+(sim$lm.res$sd.indirect[,1]*interval)
  ci_lm.df$lb_ind2<-sim$lm.res$b.indirect[,2]-(sim$lm.res$sd.indirect[,2]*interval)
  ci_lm.df$ub_ind2<-sim$lm.res$b.indirect[,2]+(sim$lm.res$sd.indirect[,2]*interval)
  
  # SLX
  ci_slx.df<-data.frame(ci.df)
  ci_slx.df$lb_dir1<-sim$slx.res$b.direct[,1]-(sim$slx.res$sd.direct[,1]*interval)
  ci_slx.df$ub_dir1<-sim$slx.res$b.direct[,1]+(sim$slx.res$sd.direct[,1]*interval)
  ci_slx.df$lb_dir2<-sim$slx.res$b.direct[,2]-(sim$slx.res$sd.direct[,2]*interval)
  ci_slx.df$ub_dir2<-sim$slx.res$b.direct[,2]+(sim$slx.res$sd.direct[,2]*interval)
  
  ci_slx.df$lb_ind1<-sim$slx.res$b.indirect[,1]-(sim$slx.res$sd.indirect[,1]*interval)
  ci_slx.df$ub_ind1<-sim$slx.res$b.indirect[,1]+(sim$slx.res$sd.indirect[,1]*interval)
  ci_slx.df$lb_ind2<-sim$slx.res$b.indirect[,2]-(sim$slx.res$sd.indirect[,2]*interval)
  ci_slx.df$ub_ind2<-sim$slx.res$b.indirect[,2]+(sim$slx.res$sd.indirect[,2]*interval)
  
  # SAR
  ci_sar.df<-data.frame(ci.df)
  ci_sar.df$lb_dir1<-sim$sar.res$b.direct[,1]-(sim$sar.res$sd.direct[,1]*interval)
  ci_sar.df$ub_dir1<-sim$sar.res$b.direct[,1]+(sim$sar.res$sd.direct[,1]*interval)
  ci_sar.df$lb_dir2<-sim$sar.res$b.direct[,2]-(sim$sar.res$sd.direct[,2]*interval)
  ci_sar.df$ub_dir2<-sim$sar.res$b.direct[,2]+(sim$sar.res$sd.direct[,2]*interval)
  
  ci_sar.df$lb_ind1<-sim$sar.res$b.indirect[,1]-(sim$sar.res$sd.indirect[,1]*interval)
  ci_sar.df$ub_ind1<-sim$sar.res$b.indirect[,1]+(sim$sar.res$sd.indirect[,1]*interval)
  ci_sar.df$lb_ind2<-sim$sar.res$b.indirect[,2]-(sim$sar.res$sd.indirect[,2]*interval)
  ci_sar.df$ub_ind2<-sim$sar.res$b.indirect[,2]+(sim$sar.res$sd.indirect[,2]*interval)
  
  # SEM
  ci_sem.df<-data.frame(ci.df)
  ci_sem.df$lb_dir1<-sim$sem.res$b.direct[,1]-(sim$sem.res$sd.direct[,1]*interval)
  ci_sem.df$ub_dir1<-sim$sem.res$b.direct[,1]+(sim$sem.res$sd.direct[,1]*interval)
  ci_sem.df$lb_dir2<-sim$sem.res$b.direct[,2]-(sim$sem.res$sd.direct[,2]*interval)
  ci_sem.df$ub_dir2<-sim$sem.res$b.direct[,2]+(sim$sem.res$sd.direct[,2]*interval)
  
  ci_sem.df$lb_ind1<-sim$sem.res$b.indirect[,1]-(sim$sem.res$sd.indirect[,1]*interval)
  ci_sem.df$ub_ind1<-sim$sem.res$b.indirect[,1]+(sim$sem.res$sd.indirect[,1]*interval)
  ci_sem.df$lb_ind2<-sim$sem.res$b.indirect[,2]-(sim$sem.res$sd.indirect[,2]*interval)
  ci_sem.df$ub_ind2<-sim$sem.res$b.indirect[,2]+(sim$sem.res$sd.indirect[,2]*interval)
  
  
  # SAC
  ci_sac.df<-data.frame(ci.df)
  ci_sac.df$lb_dir1<-sim$sac.res$b.direct[,1]-(sim$sac.res$sd.direct[,1]*interval)
  ci_sac.df$ub_dir1<-sim$sac.res$b.direct[,1]+(sim$sac.res$sd.direct[,1]*interval)
  ci_sac.df$lb_dir2<-sim$sac.res$b.direct[,2]-(sim$sac.res$sd.direct[,2]*interval)
  ci_sac.df$ub_dir2<-sim$sac.res$b.direct[,2]+(sim$sac.res$sd.direct[,2]*interval)
  
  ci_sac.df$lb_ind1<-sim$sac.res$b.indirect[,1]-(sim$sac.res$sd.indirect[,1]*interval)
  ci_sac.df$ub_ind1<-sim$sac.res$b.indirect[,1]+(sim$sac.res$sd.indirect[,1]*interval)
  ci_sac.df$lb_ind2<-sim$sac.res$b.indirect[,2]-(sim$sac.res$sd.indirect[,2]*interval)
  ci_sac.df$ub_ind2<-sim$sac.res$b.indirect[,2]+(sim$sac.res$sd.indirect[,2]*interval)
  
  # SDM
  ci_sdm.df<-data.frame(ci.df)
  ci_sdm.df$lb_dir1<-sim$sdm.res$b.direct[,1]-(sim$sdm.res$sd.direct[,1]*interval)
  ci_sdm.df$ub_dir1<-sim$sdm.res$b.direct[,1]+(sim$sdm.res$sd.direct[,1]*interval)
  ci_sdm.df$lb_dir2<-sim$sdm.res$b.direct[,2]-(sim$sdm.res$sd.direct[,2]*interval)
  ci_sdm.df$ub_dir2<-sim$sdm.res$b.direct[,2]+(sim$sdm.res$sd.direct[,2]*interval)
  
  ci_sdm.df$lb_ind1<-sim$sdm.res$b.indirect[,1]-(sim$sdm.res$sd.indirect[,1]*interval)
  ci_sdm.df$ub_ind1<-sim$sdm.res$b.indirect[,1]+(sim$sdm.res$sd.indirect[,1]*interval)
  ci_sdm.df$lb_ind2<-sim$sdm.res$b.indirect[,2]-(sim$sdm.res$sd.indirect[,2]*interval)
  ci_sdm.df$ub_ind2<-sim$sdm.res$b.indirect[,2]+(sim$sdm.res$sd.indirect[,2]*interval)
  
  # SDEM
  ci_sdem.df<-data.frame(ci.df)
  ci_sdem.df$lb_dir1<-sim$sdem.res$b.direct[,1]-(sim$sdem.res$sd.direct[,1]*interval)
  ci_sdem.df$ub_dir1<-sim$sdem.res$b.direct[,1]+(sim$sdem.res$sd.direct[,1]*interval)
  ci_sdem.df$lb_dir2<-sim$sdem.res$b.direct[,2]-(sim$sdem.res$sd.direct[,2]*interval)
  ci_sdem.df$ub_dir2<-sim$sdem.res$b.direct[,2]+(sim$sdem.res$sd.direct[,2]*interval)
  
  ci_sdem.df$lb_ind1<-sim$sdem.res$b.indirect[,1]-(sim$sdem.res$sd.indirect[,1]*interval)
  ci_sdem.df$ub_ind1<-sim$sdem.res$b.indirect[,1]+(sim$sdem.res$sd.indirect[,1]*interval)
  ci_sdem.df$lb_ind2<-sim$sdem.res$b.indirect[,2]-(sim$sdem.res$sd.indirect[,2]*interval)
  ci_sdem.df$ub_ind2<-sim$sdem.res$b.indirect[,2]+(sim$sdem.res$sd.indirect[,2]*interval)
  
  
  ### Coverage rates
  
  ## Direct
  # LM
  mat_d[2,1]<-length(which(ci_lm.df$lb_dir1<=dx1 & ci_lm.df$ub_dir1>=dx1))/N
  mat_d[3,1]<-length(which(ci_lm.df$lb_dir2<=dx2 & ci_lm.df$ub_dir2>=dx2))/N
  
  # SLX
  mat_d[5,1]<-length(which(ci_slx.df$lb_dir1<=dx1 & ci_slx.df$ub_dir1>=dx1))/N
  mat_d[6,1]<-length(which(ci_slx.df$lb_dir2<=dx2 & ci_slx.df$ub_dir2>=dx2))/N
  
  # SAR
  mat_d[8,1]<-length(which(ci_sar.df$lb_dir1<=dx1 & ci_sar.df$ub_dir1>=dx1))/N
  mat_d[9,1]<-length(which(ci_sar.df$lb_dir2<=dx2 & ci_sar.df$ub_dir2>=dx2))/N
  
  # SEM
  mat_d[11,1]<-length(which(ci_sem.df$lb_dir1<=dx1 & ci_sem.df$ub_dir1>=dx1))/N
  mat_d[12,1]<-length(which(ci_sem.df$lb_dir2<=dx2 & ci_sem.df$ub_dir2>=dx2))/N
  
  # SAC
  mat_d[14,1]<-length(which(ci_sac.df$lb_dir1<=dx1 & ci_sac.df$ub_dir1>=dx1))/N
  mat_d[15,1]<-length(which(ci_sac.df$lb_dir2<=dx2 & ci_sac.df$ub_dir2>=dx2))/N
  
  # SDM
  mat_d[17,1]<-length(which(ci_sdm.df$lb_dir1<=dx1 & ci_sdm.df$ub_dir1>=dx1))/N
  mat_d[18,1]<-length(which(ci_sdm.df$lb_dir2<=dx2 & ci_sdm.df$ub_dir2>=dx2))/N
  
  # SDEM
  mat_d[20,1]<-length(which(ci_sdem.df$lb_dir1<=dx1 & ci_sdem.df$ub_dir1>=dx1))/N
  mat_d[21,1]<-length(which(ci_sdem.df$lb_dir2<=dx2 & ci_sdem.df$ub_dir2>=dx2))/N
  
  
  ## Indirect
  # # LM
  # mat_sp[2,1]<-length(which(ci_lm.df$lb_ind1<=spx1 & ci_lm.df$ub_ind1>=spx1))/N
  # mat_sp[3,1]<-length(which(ci_lm.df$lb_ind2<=spx2 & ci_lm.df$ub_ind2>=spx2))/N
  
  # SLX
  mat_sp[5,1]<-length(which(ci_slx.df$lb_ind1<=spx1 & ci_slx.df$ub_ind1>=spx1))/N
  mat_sp[6,1]<-length(which(ci_slx.df$lb_ind2<=spx2 & ci_slx.df$ub_ind2>=spx2))/N
  
  # SAR
  mat_sp[8,1]<-length(which(ci_sar.df$lb_ind1<=spx1 & ci_sar.df$ub_ind1>=spx1))/N
  mat_sp[9,1]<-length(which(ci_sar.df$lb_ind2<=spx2 & ci_sar.df$ub_ind2>=spx2))/N
  
  # # SEM
  # mat_sp[11,1]<-length(which(ci_sem.df$lb_ind1<=spx1 & ci_sem.df$ub_ind1>=spx1))/N
  # mat_sp[12,1]<-length(which(ci_sem.df$lb_ind2<=spx2 & ci_sem.df$ub_ind2>=spx2))/N
  
  # SAC
  mat_sp[14,1]<-length(which(ci_sac.df$lb_ind1<=spx1 & ci_sac.df$ub_ind1>=spx1))/N
  mat_sp[15,1]<-length(which(ci_sac.df$lb_ind2<=spx2 & ci_sac.df$ub_ind2>=spx2))/N
  
  # SDM
  mat_sp[17,1]<-length(which(ci_sdm.df$lb_ind1<=spx1 & ci_sdm.df$ub_ind1>=spx1))/N
  mat_sp[18,1]<-length(which(ci_sdm.df$lb_ind2<=spx2 & ci_sdm.df$ub_ind2>=spx2))/N
  
  # SDEM
  mat_sp[20,1]<-length(which(ci_sdem.df$lb_ind1<=spx1 & ci_sdem.df$ub_ind1>=spx1))/N
  mat_sp[21,1]<-length(which(ci_sdem.df$lb_ind2<=spx2 & ci_sdem.df$ub_ind2>=spx2))/N
  
  
  
  
  
  res<-list(mat_d=mat_d, mat_sp=mat_sp,
            theta=sim$theta, parameters=sim$parameters)

  
  res
}





print.summary.sim <- function(x, digits = 5, width=getOption("width"), stars=T, ...){


  ## Make printable mat direct
  mat_d<-x$mat_d
  mat_d2<-round(mat_d, digits)
  mat_d2<-apply(mat_d2, 2, FUN=function(x) as.character(format(x, digits=digits, scientific=F)))
  mat_d2[seq(1, nrow(mat_d), 3),]<-c("","")
  rownames(mat_d2)<-rownames(mat_d)
  
  ## Add star for best result
  if(stars){
    p11<-which(abs(mat_d[,1])==min(abs(mat_d[seq(2, nrow(mat_d), 3),1])))
    p12<-which(abs(mat_d[,1])==min(abs(mat_d[seq(3, nrow(mat_d), 3),1])))
    
    p21<-which(abs(mat_d[,2])==min(abs(mat_d[seq(2, nrow(mat_d), 3),2])))
    p22<-which(abs(mat_d[,2])==min(abs(mat_d[seq(3, nrow(mat_d), 3),2])))  
    
    mat_d2[c(p11, p12),1]<-paste(mat_d2[c(p11, p12),1], "*", sep="")
    mat_d2[c(p21, p22),2]<-paste(mat_d2[c(p21, p22),2], "*", sep="")
    
  }
  
  ## Make printable mat indirect
  mat_sp<-x$mat_sp
  mat_sp2<-round(mat_sp, digits)
  mat_sp2<-apply(mat_sp2, 2, FUN=function(x) as.character(format(x, digits=digits, scientific=F)))
  mat_sp2[seq(1, nrow(mat_sp), 3),]<-c("","")
  rownames(mat_sp2)<-rownames(mat_sp)
  
  ## Add star for best result
  if(stars){
    p11<-which(abs(mat_sp[,1])==min(abs(mat_sp[seq(2, nrow(mat_sp), 3),1]), na.rm=T))
    p12<-which(abs(mat_sp[,1])==min(abs(mat_sp[seq(3, nrow(mat_sp), 3),1]), na.rm=T))
    
    p21<-which(abs(mat_sp[,2])==min(abs(mat_sp[seq(2, nrow(mat_sp), 3),2]), na.rm=T))
    p22<-which(abs(mat_sp[,2])==min(abs(mat_sp[seq(3, nrow(mat_sp), 3),2]), na.rm=T))  
    
    mat_sp2[c(p11, p12),1]<-paste(mat_sp2[c(p11, p12),1], "*", sep="")
    mat_sp2[c(p21, p22),2]<-paste(mat_sp2[c(p21, p22),2], "*", sep="")
    
  }

  cat("\n")
  cat("\nCall:\n")
  print(x$parameters)
  cat("\n")
  cat("\nImpacts:\n")
  print(x$theta)
  cat("\n")
  
  matrix.list <- list(mat_d2, mat_sp2)
  matrix.chain <- do.call(cbind, lapply(matrix.list, fn))
  cat("\nResults direct:\n")
  cat(" ", paste0(c("Direct", "Indirect"), collapse = "                       "), "\n"); print(matrix.chain, row.names = FALSE)
  cat("\n")
  
  x$mat_d2<-mat_d2
  x$mat_sp2<-mat_sp2
  
  invisible(x)
}


plot.sim<-function(sim=NULL, xlim=NULL,..){
  
  x<-summary.sim(sim)
  n<-(nrow(x$mat_d)/3)
  
  ### Get true parameters
  dx1<-as.numeric(sim$theta[3,3])
  dx2<-as.numeric(sim$theta[4,3])
  
  spx1<-as.numeric(sim$theta[3,4])
  spx2<-as.numeric(sim$theta[4,4])
  
  ### Set up model Frame
  mf<-data.frame(matrix(NA, ncol=4, nrow=n*4))
  colnames(mf)<-c("Model", "Coefficient", "SE", "modelName")
  names<-rownames(x$mat_d)[which(is.na(x$mat_d[,1]))]
  mf$Model<-rep(names,4)
  mf$modelName<-c(rep("Direct x1",n), rep("Direct x2",n),
                  rep("Indirect x1",n), rep("Indirect x2",n))
  
  ### Paste bias values
  mat_d<-x$mat_d
  mat_sp<-x$mat_sp
  mat_sd<-x$mat_sd
  
  ### Bias
  # Direct
  mf[1:n,2]<-mat_d[seq(2, nrow(mat_d), by=3),1]
  mf[(n+1):(2*n),2]<-mat_d[seq(3, nrow(mat_d), by=3),1]
  
  # Indirect
  mf[(2*n+1):(3*n),2]<-mat_sp[seq(2, nrow(mat_sp), by=3),1]
  mf[(3*n+1):(4*n),2]<-mat_sp[seq(3, nrow(mat_sp), by=3),1]
  
  ### Standard-deviation
  # Direct
  mf[1:n,3]<-mat_sd[seq(2, nrow(mat_sd), by=3),1]
  mf[(n+1):(2*n),3]<-mat_sd[seq(3, nrow(mat_sd), by=3),1]
  
  # Indirect
  mf[(2*n+1):(3*n),3]<-mat_sd[seq(2, nrow(mat_sd), by=3),2]
  mf[(3*n+1):(4*n),3]<-mat_sd[seq(3, nrow(mat_sd), by=3),2]
  
  
  # Reorder levels
  mf$Model<-as.factor(mf$Model)
  mf$modelName<-as.factor(mf$modelName)
  
  mf$Model<-factor(mf$Model, 
                    levels(mf$Model)[rev(c(1, 7, 3, 6, 2, 5, 4))])
  mf$modelName<-factor(mf$modelName, 
                        levels(mf$modelName)[rev(c(1,2,3,4))])
  
  
  # Confidence intervals
  interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier
  
  
  # Plot
  zp1 <- ggplot(mf, aes(colour = modelName, shape = modelName))
  zp1 <- zp1 + geom_hline(yintercept = 0, colour = "black", lty = 2)
  zp1 <- zp1 + scale_shape_manual(values=c(15,16,17,18)) + scale_size_manual(values=0.5)
  zp1 <- zp1 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                   ymax = Coefficient + SE*interval2, 0.2),
                               lwd = 1, position = position_dodge(width = 1/2),
                               fill ="black")
  zp1 <- zp1 + coord_flip() + theme_bw()
  if(!is.null(xlim)){
    zp1 <- zp1 + ylim(xlim[1], xlim[2])
  }
  zp1 <- zp1 + theme(legend.title = element_blank())
  zp1 <- zp1 + labs(y="Bias") 
  zp1 <- zp1 + theme(legend.text = element_text(size = 16),
                     legend.position="bottom",
                     legend.key=element_blank(),
                     axis.text.x = element_text(size=16),
                     axis.text.y = element_text(size=16),
                     axis.title.x = element_text(size=16),
                     axis.title.y = element_text(size=16),
                     plot.title = element_text(size=16))
  zp1 <- zp1 + ggtitle(element_blank()) 
  zp1 <- zp1 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                      shape = guide_legend(reverse = T))
  print(zp1)
  
  invisible(zp1)
  
}




# ##################
# #### Analysis ####
# ##################
# 
# 
# sim1<-spsim(N=900, R=1000, seed=123789,
#             ay=0, ax1=0.5, ax2=0, au=0,
#             bx1=0.7, bwx1=0, bx2=1.5, bwx2=0,
#             bu=1, bwu=0,
#             omv1=0, womv1=0)
# save(sim1, file="Simulation_1.RData")
# summary.sim(sim1)
# plot.sim(sim1)
# 
# 
# sim2<-spsim(N=900, R=1000, seed=123789,
#             ay=0.2, ax1=0, ax2=0, au=0,
#             bx1=0.7, bwx1=0, bx2=1.5, bwx2=0,
#             bu=1, bwu=0,
#             omv1=0, womv1=0)
# save(sim2, file="Simulation_2.RData")
# summary.sim(sim2)
# plot.sim(sim2)
# 
# sim3<-spsim(N=900, R=1000, seed=123789,
#             ay=0.2, ax1=0.3, ax2=0, au=0,
#             bx1=0.7, bwx1=0, bx2=1.5, bwx2=0,
#             bu=1, bwu=0,
#             omv1=0, womv1=0)
# save(sim3, file="Simulation_3.RData")
# summary.sim(sim3)
# plot.sim(sim3)
# 
# 
# sim4<-spsim(N=900, R=1000, seed=123789,
#             ay=0.2, ax1=0.3, ax2=0.3, au=0,
#             bx1=0.7, bwx1=0, bx2=1.5, bwx2=0,
#             bu=1, bwu=0,
#             omv1=0, womv1=0)
# save(sim4, file="Simulation_4.RData")
# summary.sim(sim4)
# plot.sim(sim4)
# 
# 
# sim5<-spsim(N=900, R=1000, seed=123789,
#             ay=0.2, ax1=0.3, ax2=0.6, au=0,
#             bx1=0.7, bwx1=0, bx2=1.5, bwx2=0,
#             bu=1, bwu=0,
#             omv1=0, womv1=0)
# save(sim5, file="Simulation_5.RData")
# summary.sim(sim5)
# plot.sim(sim5)
# 
# sim6<-spsim(N=900, R=1000, seed=123789,
#             ay=0.2, ax1=0.3, ax2=0.6, au=0.2,
#             bx1=0.7, bwx1=0, bx2=1.5, bwx2=0,
#             bu=1, bwu=0,
#             omv1=0, womv1=0)
# save(sim6, file="Simulation_6.RData")
# summary.sim(sim6)
# plot.sim(sim6)
# 
# sim7<-spsim(N=900, R=1000, seed=123789,
#             ay=0.2, ax1=0, ax2=0, au=0.2,
#             bx1=0.7, bwx1=0, bx2=1.5, bwx2=0,
#             bu=1, bwu=0,
#             omv1=0, womv1=0)
# save(sim7, file="Simulation_7.RData")
# summary.sim(sim7)
# plot.sim(sim7)
# 
# sim8<-spsim(N=900, R=1000, seed=123789,
#             ay=0, ax1=0.3, ax2=0.6, au=0.2,
#             bx1=0.7, bwx1=0, bx2=1.5, bwx2=0,
#             bu=1, bwu=0,
#             omv1=0, womv1=0)
# save(sim8, file="Simulation_8.RData")
# summary.sim(sim8)
# plot.sim(sim8)
# 
# sim9<-spsim(N=900, R=1000, seed=123789,
#             ay=0, ax1=0, ax2=0, au=0.2,
#             bx1=0.7, bwx1=0, bx2=1.5, bwx2=0,
#             bu=1, bwu=0,
#             omv1=0, womv1=0)
# save(sim9, file="Simulation_9.RData")
# summary.sim(sim9)
# plot.sim(sim9)
# 
# 
# 
# 
# ### Local Spillovers
# 
# 
# sim10<-spsim(N=900, R=1000, seed=123789,
#             ay=0, ax1=0, ax2=0, au=0,
#             bx1=0.7, bwx1=0.2, bx2=1.5, bwx2=0.2,
#             bu=1, bwu=0,
#             omv1=0, womv1=0)
# save(sim10, file="Simulation_10.RData")
# summary.sim(sim10)
# plot.sim(sim10)
# 
# sim11<-spsim(N=900, R=1000, seed=123789,
#              ay=0, ax1=0, ax2=0, au=0,
#              bx1=0.7, bwx1=0.2, bx2=1.5, bwx2=0.4,
#              bu=1, bwu=0,
#              omv1=0, womv1=0)
# save(sim11, file="Simulation_11.RData")
# summary.sim(sim11)
# plot.sim(sim11)
# 
# sim12<-spsim(N=900, R=1000, seed=123789,
#              ay=0, ax1=0, ax2=0, au=0.2,
#              bx1=0.7, bwx1=0.2, bx2=1.5, bwx2=0.4,
#              bu=1, bwu=0,
#              omv1=0, womv1=0)
# save(sim12, file="Simulation_12.RData")
# summary.sim(sim12)
# plot.sim(sim12)
# 
# 
# sim13<-spsim(N=900, R=1000, seed=123789,
#              ay=0, ax1=0, ax2=0.5, au=0,
#              bx1=0.7, bwx1=0.2, bx2=1.5, bwx2=0.4,
#              bu=1, bwu=0,
#              omv1=0, womv1=0)
# save(sim13, file="Simulation_13.RData")
# summary.sim(sim13)
# plot.sim(sim13)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# test<-spsim(N=900, R=30, seed=123789,
#             ay=0.2, ax1=0.3, ax2=0.6, au=0.2,
#             bx1=-0.7, bwx1=0, bx2=1.5, bwx2=0,
#             bu=1, bwu=0,
#             omv1=0, womv1=0)
# summary.sim(test)
# plot.sim(test)









