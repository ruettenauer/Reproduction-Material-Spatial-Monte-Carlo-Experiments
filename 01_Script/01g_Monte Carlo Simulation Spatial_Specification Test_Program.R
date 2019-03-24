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

# setwd("C:/Users/tor44za/Uni Würzburg/Forschung/Monte Carlo Spatial/02_Daten")
# 
# 
### Parameters

# N=900
# R=30
# seed=123579
# tol=1e-10
# ay=0
# ax1=0
# ax2=0
# au=0
# bx1=0.4
# bwx1=0
# bx2=0.6
# bwx2=0
# bu=1
# bwu=0
# sdu=0.5
# omv1=0
# womv1=0
# omv2=0
# womv2=0



spsim_lmtest<-function(N=900, R=30, seed=123579, tol=1e-10,
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
  
  
  ### Create spatial listw object (Queens Nbs)
  df.nb<-poly2nb(df.sp, queen=T)
  df.listw<-nb2listw(df.nb, style="W")
  # trW<-trW(W=as(df.listw, "CsparseMatrix"), type="mult")
  
  
  
  ###########################
  ### Compute real values ###
  ###########################
  
  ### Contruct Multiplyer mat
  W<-listw2mat(df.listw)
  
  M<-solve(diag(nrow(df.sp))-ay*W)
  
  Mx1<-solve(diag(nrow(df.sp))-ax1*W)
  Mx2<-solve(diag(nrow(df.sp))-ax2*W)
  
  Mu<-solve(diag(nrow(df.sp))-au*W)
  

  # ### Calculate Impacts
  # spMat1<- (M*bx1 + M %*% (W*bwx1))
  # deff1<-mean(diag(spMat1))
  # diag(spMat1)<-NA
  # speff1<-mean(rowSums(spMat1, na.rm=T))
  # 
  # spMat2<- (M*bx2 + M %*% (W*bwx2))
  # deff2<-mean(diag(spMat2))
  # diag(spMat2)<-NA
  # speff2<-mean(rowSums(spMat2, na.rm=T))
  # 
  # ### Save values
  # theta<-matrix(NA, ncol=4, nrow=4)
  # colnames(theta)<-c("", "Coef", "Direct", "Indirect")
  # theta[,1]<-c("rho", "lambda", "beta1", "beta2")
  # theta[,2]<-c(ay, au, bx1, bx2)
  # theta[3:4,3]<-c(deff1, deff2)
  # theta[3:4,4]<-c(speff1, speff2)
  
  
  ##############################
  ### Set up result matrices ###
  ##############################
  

  test.res<-vector("list", R)
  
  
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
    
    ### Linear OLS
    lm.mod<-lm(Y ~ X1 + X2, data=data.frame(df.sp) )
    
    
    ### LM tests
    LMtest.res<-lm.LMtests(lm.mod, df.listw, test=c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA"))
    
    
    ### Paste into results object
    test.res[[i]]<-LMtest.res
    
    
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
               seeds      = seeds,
               LMtest.res = test.res)
  
  class(result) <- c("spsim_lmtest")
  
  return(result)
  
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

comb<-function(x){
  lapply(seq_along(x[[1]]), function(i) 
    do.call(Map, c(f = rbind, lapply(x, `[[`, i))))
}


summary.spsim_lmtest<-function(sim=NULL, digits=5, ...){
  
  test.res<-sim$LMtest.res
  
  ### Set up results matrix
  mat<-matrix(nrow=2, ncol=5)
  colnames(mat)<-c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA")
  rownames(mat)<-c("Statistic", "p<=0.05")
  
  
  ### Get statistic
  mat[1,1]<-mean(sapply(sapply(test.res, "[[", 1, simplify = F), "[[", 1))
  mat[1,2]<-mean(sapply(sapply(test.res, "[[", 2, simplify = F), "[[", 1))
  mat[1,3]<-mean(sapply(sapply(test.res, "[[", 3, simplify = F), "[[", 1))
  mat[1,4]<-mean(sapply(sapply(test.res, "[[", 4, simplify = F), "[[", 1))
  mat[1,5]<-mean(sapply(sapply(test.res, "[[", 5, simplify = F), "[[", 1))
  
  ### Get p values
  p1<-sapply(sapply(test.res, "[[", 1, simplify = F), "[[", 3)
  p2<-sapply(sapply(test.res, "[[", 2, simplify = F), "[[", 3)
  p3<-sapply(sapply(test.res, "[[", 3, simplify = F), "[[", 3)
  p4<-sapply(sapply(test.res, "[[", 4, simplify = F), "[[", 3)
  p5<-sapply(sapply(test.res, "[[", 5, simplify = F), "[[", 3)
  
  mat[2,1]<-length(which(p1<=0.05))/length(p1)
  mat[2,2]<-length(which(p2<=0.05))/length(p2)
  mat[2,3]<-length(which(p3<=0.05))/length(p3)
  mat[2,4]<-length(which(p4<=0.05))/length(p4)
  mat[2,5]<-length(which(p5<=0.05))/length(p5)
  
  
  res<-list(mat=mat, parameters=sim$parameters, digits=digits)
  
  class(res) <- c("summary.spsim_lmtest")
  
  res
}

print.summary.spsim_lmtest <- function(x, digits = 5, width=getOption("width")){
  
  digits<-x$digits
  mat<-x$mat
  names<-rownames(mat)
  ## Make printable mat direct

  mat<-apply(mat, 2, FUN=function(x) as.character(format(x, digits=digits, scientific=F)))
  rownames(mat)<-names
  

  cat("\n")
  cat("\nCall:\n")
  print(x$parameters)
  cat("\n")
  print(mat, row.names = TRUE, quote=FALSE)
  cat("\n")
  
  
  invisible(x)
}





