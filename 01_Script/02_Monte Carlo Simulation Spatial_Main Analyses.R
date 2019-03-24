
###########################################
#### Rproduction File Spatial Analyses ####
####        Tobias Ruettenauer         ####
####          Main Analyses            ####
###########################################

rm(list = ls())

#### Specify directory ####

dir <- "C:/Users/tor44za/Uni Würzburg/Forschung/Monte Carlo Spatial/09_Reproduction-Material"

#### Set up directories ####

dd<-paste(dir, "/01_Script", sep="")
wd<-paste(dir, "/02_Data", sep="")
od<-paste(dir, "/03_Output/", sep="")


#### Start Code ####

setwd(wd)

#### Load packages ####

library(ggmap)
library(RgoogleMaps)
library(foreign)
library(sp)
library(GISTools)
library(sp)
library(spdep)
library(texreg)
library(xtable)

library(rgdal)
library(raster)
library(extrafont)
loadfonts()

### Load program
source(paste(dd, "01_Monte Carlo Simulation Spatial_Program.R", sep="/"))



######################################
######################################
#### 1) Main simulation results   ####
######################################
######################################

######################################
#### Define possible combinations ####
######################################

beta<-c(0.2, 0.5)

rho<-c(0, 0.5)
delta<-c("0, 0", "0.4, 0.7")
lambda<-c(0, 0.5)
gamma<-c(0, 0.3)
theta<-c("0, 0", "0.1, 0.8")


### Set up matrix of all combinations ###

nc<-2^5

l<-list(rho, delta, lambda, gamma, theta)

comb_mat4<-expand.grid(l)

colnames(comb_mat4)<-c("rho", "delta", "lambda", "gamma", "theta")
comb_mat4$delta<-as.character(comb_mat4$delta)
comb_mat4$theta<-as.character(comb_mat4$theta)

save(comb_mat4, file="Combinations4.RData")





###########################################
#### Start simulations of combinations ####
###########################################

# Set up manual loops (to spare capacity)


k<-c(1:32)

### Simumlation Run
for(i in k){

  sim4<-NA # avoid using the last sim (in case of error)

  par<-comb_mat4[i,]
  deltai<-as.numeric(c(unlist(strsplit(par$delta, split=", "))))
  thetai<-as.numeric(c(unlist(strsplit(par$theta, split=", "))))

  sim4<-spsim(N=900, R=1000, seed=1573262181,
              ay=par$rho, ax1=deltai[1], ax2=deltai[2], au=par$lambda,
              bx1=beta[1], bwx1=thetai[1], bx2=beta[2], bwx2=thetai[2],
              bu=1, bwu=0, sdu=1,
              omv1=par$gamma, womv1=0)

  assign(paste0("sim4_", i, sep=""), sim4)


  save(list=c(paste0("sim4_", i, sep="")), file=paste0("Simulation4_", i, ".RData", sep=""))

}


# In case you reload results:
# for(i in k){
#   load(paste("Simulation4_", i, ".RData", sep=""))
# }




#######################################
#### Make plots of each Simulation ####
#######################################


### Jpeg ###

for(i in k){
  cm<-comb_mat4[i,]
  
  tmp.sim<-get(paste("sim4_", i, sep=""))
  
  p<-plot.sim(tmp.sim, xlim=c(-2, 3))
  
  name = paste(od, "sim4_", i, ".jpeg", sep="")
  jpeg(file = name, bg = "transparent", width = 2720, height = 1820, quality=100, family="CM Roman", res=300)
  par(mar=c(0,0,0,0))
  par(mfrow=c(1,1),oma=c(0,0,0,0))
  
  print(p + ggtitle(bquote(rho == .(cm$rho) ~~
                             delta == .(cm$delta) ~~
                             lambda == .(cm$lambda) ~~
                             gamma == .(cm$gamma) ~~
                             theta == .(cm$theta))))
  dev.off()
  par(mfrow=c(1,1))
}






#####################################
#### Create result matrices RMSE ####
#####################################


### Direct impacts ###

res4_rmse_dir<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res4_rmse_dir)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res4_rmse_dir<-cbind(res4_rmse_dir,tmp)

for(i in 1:nrow(res4_rmse_dir)){
  tmp.sim<-get(paste("sim4_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_d[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res4_rmse_dir[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res4_rmse_dir[i, 6:ncol(res4_rmse_dir)]<-t(m[,2])
}



### Separate into omv and non-omv

res4_rmse_dir1<-res4_rmse_dir[res4_rmse_dir$gamma==0,]
res4_rmse_dir2<-res4_rmse_dir[res4_rmse_dir$gamma>0,]

res4_rmse_dir1<-res4_rmse_dir1[,-which(names(res4_rmse_dir1)=="gamma")]
res4_rmse_dir2<-res4_rmse_dir2[,-which(names(res4_rmse_dir2)=="gamma")]


res4_rmse_dir1[,1:4]<-apply(res4_rmse_dir1[,1:4], 2, FUN=function(x) as.character(x))
res4_rmse_dir2[,1:4]<-apply(res4_rmse_dir2[,1:4], 2, FUN=function(x) as.character(x))





### Indirect impacts ###

res4_rmse_ind<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res4_rmse_ind)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res4_rmse_ind<-cbind(res4_rmse_ind,tmp)

for(i in 1:nrow(res4_rmse_ind)){
  tmp.sim<-get(paste("sim4_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sp[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res4_rmse_ind[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res4_rmse_ind[i, 6:ncol(res4_rmse_ind)]<-t(m[,2])
}


res4_rmse_ind<-res4_rmse_ind[,-which(names(res4_rmse_ind) %in% c("OLS_x1", "OLS_x2", "SEM_x1", "SEM_x2"))]

### Separate into omv and non-omv

res4_rmse_ind1<-res4_rmse_ind[res4_rmse_ind$gamma==0,]
res4_rmse_ind2<-res4_rmse_ind[res4_rmse_ind$gamma>0,]

res4_rmse_ind1<-res4_rmse_ind1[,-which(names(res4_rmse_ind1)=="gamma")]
res4_rmse_ind2<-res4_rmse_ind2[,-which(names(res4_rmse_ind2)=="gamma")]


res4_rmse_ind1[,1:4]<-apply(res4_rmse_ind1[,1:4], 2, FUN=function(x) as.character(x))
res4_rmse_ind2[,1:4]<-apply(res4_rmse_ind2[,1:4], 2, FUN=function(x) as.character(x))



#######################
### Wirte tex table ###
#######################


matout<-function(x, file="", digits=NULL, colav=T, threepart=T, stars=F,
                 code_before=NULL, code_after=NULL,
                 header1=NULL, header2=NULL, 
                 caption=NULL, label=NULL, note=NULL,...) {
  
  ### Options for xtable
  default_args = list(include.colnames=FALSE, only.contents=TRUE,
                      include.rownames=FALSE, hline.after=NULL, comment=FALSE,
                      print.results=FALSE, sanitize.text.function = function(z) z)
  passed_args = list(...)
  
  
  if(is.null(caption)){caption<-"Table"}
  
  c<-ncol(x)
  
  if(is.null(code_before)){
    code_before<-paste("\\begin{tabular}{", paste(rep("l", c), collapse = " "), "}")
  }
  if(is.null(code_after)){
    code_after<-paste("\\end{tabular}")
  }
  
  num<-which(sapply(x,FUN=function(x) is.numeric(x)))
  
  if(colav==T){
    av<-data.frame(t(colMeans(abs(x[,num[1]:ncol(x)]))))
    av2<-data.frame(t(colMeans(matrix(t(colMeans(abs(x[,num[1]:ncol(x)]))), 2))))
    
    if(stars){
      pa1<-colnames(av[seq(1, ncol(av), 2)])[which.min(abs(av[seq(1, ncol(av), 2)]))]
      pa2<-colnames(av[seq(2, ncol(av), 2)])[which.min(abs(av[seq(2, ncol(av), 2)]))]
      pa3<-colnames(av2)[which.min(av2)]
      
      #av[1,]<-round(av[1,], digits)
      av[1,]<-sprintf(paste("%.", digits, "f"), round(av[1,],digits))
      av2[1,]<-sprintf(paste("%.", digits, "f"), round(av2[1,],digits))
      
      
      av[1,colnames(av)==pa1]<-paste(av[1,colnames(av)==pa1], "\\ensuremath{^{*}}", sep="")
      av[1,colnames(av)==pa2]<-paste(av[1,colnames(av)==pa2], "\\ensuremath{^{*}}", sep="")
      av2[1,colnames(av2)==pa3]<-paste(av2[1,colnames(av2)==pa3], "\\ensuremath{^{*}}", sep="")
      
      
    }
    
    calling_args2 = c(list(x=xtable(av, digits=digits)),
                      c(passed_args,
                        default_args[setdiff(names(default_args), names(passed_args))]))
    
  }
  
  if(stars){
    p11<-apply(x[,seq(num[1], ncol(x), 2)], 1, FUN=function(z) which.min(abs(z)))
    p11<-colnames(x[,seq(num[1], ncol(x), 2)])[p11]
    
    p12<-apply(x[,seq(num[1]+1, ncol(x), 2)], 1, FUN=function(z) which.min(abs(z)))
    p12<-colnames(x[,seq(num[1]+1, ncol(x), 2)])[p12]
    
    #x[, num]<-round(x[, num], digits)
    #x[, num]<-apply(x[, num], 2, FUN=function(z) as.character(format(z, digits=digits, scientific=F)))
    x[, num]<-apply(x[, num], 2, FUN=function(z) sprintf(paste("%.", digits, "f"), round(z,digits)))
    
    
    for(i in 1:nrow(x)){
      x[i,colnames(x)==p11[i]]<-paste(x[i,colnames(x)==p11[i]], "\\ensuremath{^{*}}", sep="")
      x[i,colnames(x)==p12[i]]<-paste(x[i,colnames(x)==p12[i]], "\\ensuremath{^{*}}", sep="")
    }
    
    
  }
  
  calling_args = c(list(x=xtable(x, digits=digits)),
                   c(passed_args,
                     default_args[setdiff(names(default_args), names(passed_args))]))
  
  midrule<-paste(seq(num[1], max(num), 2), "-", seq(num[1]+1, max(num), 2), sep="")
  midrule<-c(paste("\\cmidrule(r){", midrule[1], "}", sep=""), 
             paste("\\cmidrule(lr){", midrule[-c(1, length(midrule))], "}", sep=""),
             paste("\\cmidrule(l){", midrule[length(midrule)], "}", sep=""))
  
  
  ### Output
  cat("\\begin{table}\n", file=file)
  cat("\\centering\n", file=file, append=T)
  
  if(threepart==T){
    cat("{\\begin{threeparttable}\n", file=file, append=T)
  }
  
  cat(paste0("\\caption{", caption, "}\n"), file=file, append=T)
  
  if(!is.null(label)){
    cat(paste0("\\label{", label, "}\n"), file=file, append=T)
  }
  
  cat(paste0(code_before, "\n"), file=file, append=T)
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  
  if(!is.null(header1)){
    cat(paste(header1), "\\\\ \n", file=file, append=T)
  }else{cat(paste(paste("M", seq(1, c,), collapse=" & "), "\\\\ \n"), file=file, append=T)}
  
  cat(midrule, file=file, append=T)
  
  if(!is.null(header2)){
    cat(paste(header2), "\\\\ \n", file=file, append=T)
  }
  
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  cat(do.call(print.xtable, calling_args), file=file, append=T)
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  
  if(colav==T){
    cat(paste("\\multicolumn{", (num[1]-1), "}{l}{Average} &", collapse=""), file=file, append=T)
    cat(do.call(print.xtable, calling_args2), file=file, append=T)
    cat(paste("\\multicolumn{", (num[1]-1), "}{l}{} &", collapse=""), file=file, append=T)
    cat(paste("\\multicolumn{2}{c}{", av2, "}", collapse=" & "), "\\\\ \n", file=file, append=T)
    cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  }
  
  cat(paste0(code_after, "\n"), file=file, append=T)
  
  if(!is.null(note)){
    cat("\\begin{tablenotes}\n", file=file, append=T)
    cat(paste("\\item \\scriptsize{", note, "}\n", collape=""), file=file, append=T)
    cat("\\end{tablenotes}\n", file=file, append=T)
  }
  
  if(threepart==T){
    cat("\\end{threeparttable}}\n", file=file, append=T)
  }
  
  cat("\\end{table}\n", file=file, append=T)
}



### Table 1

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 7), collapse=" & "), sep=" & " )

matout(res4_rmse_dir1, file=paste(od, "res4_rmse_dir1.tex", sep=""), digits=4, stars=T,
       caption="RMSE of direct impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res4_dir1",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 2

matout(res4_rmse_dir2, file=paste(od, "res4_rmse_dir2.tex", sep=""), digits=4, stars=T,
       caption="RMSE of direct impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res4_dir2",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")



### Table 3

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("SLX", "SAR", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 5), collapse=" & "), sep=" & " )

matout(res4_rmse_ind1, file=paste(od, "res4_rmse_ind1.tex", sep=""), digits=4, stars=T,
       caption="RMSE of indirect impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res4_ind1",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 4

matout(res4_rmse_ind2, file=paste(od, "res4_rmse_ind2.tex", sep=""), digits=4, stars=T,
       caption="RMSE of indirect impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res4_ind2",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")





#####################################
#### Create result matrices Bias ####
#####################################


### Direct impacts ###

res4_bias_dir<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res4_bias_dir)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res4_bias_dir<-cbind(res4_bias_dir,tmp)

for(i in 1:nrow(res4_bias_dir)){
  tmp.sim<-get(paste("sim4_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_d[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res4_bias_dir[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res4_bias_dir[i, 6:ncol(res4_bias_dir)]<-t(m[,1])
}



### Separate into omv and non-omv

res4_bias_dir1<-res4_bias_dir[res4_bias_dir$gamma==0,]
res4_bias_dir2<-res4_bias_dir[res4_bias_dir$gamma>0,]

res4_bias_dir1<-res4_bias_dir1[,-which(names(res4_bias_dir1)=="gamma")]
res4_bias_dir2<-res4_bias_dir2[,-which(names(res4_bias_dir2)=="gamma")]


res4_bias_dir1[,1:4]<-apply(res4_bias_dir1[,1:4], 2, FUN=function(x) as.character(x))
res4_bias_dir2[,1:4]<-apply(res4_bias_dir2[,1:4], 2, FUN=function(x) as.character(x))





### Indirect impacts ###

res4_bias_ind<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res4_bias_ind)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res4_bias_ind<-cbind(res4_bias_ind,tmp)

for(i in 1:nrow(res4_bias_ind)){
  tmp.sim<-get(paste("sim4_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sp[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res4_bias_ind[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res4_bias_ind[i, 6:ncol(res4_bias_ind)]<-t(m[,1])
}


res4_bias_ind<-res4_bias_ind[,-which(names(res4_bias_ind) %in% c("OLS_x1", "OLS_x2", "SEM_x1", "SEM_x2"))]

### Separate into omv and non-omv

res4_bias_ind1<-res4_bias_ind[res4_bias_ind$gamma==0,]
res4_bias_ind2<-res4_bias_ind[res4_bias_ind$gamma>0,]

res4_bias_ind1<-res4_bias_ind1[,-which(names(res4_bias_ind1)=="gamma")]
res4_bias_ind2<-res4_bias_ind2[,-which(names(res4_bias_ind2)=="gamma")]


res4_bias_ind1[,1:4]<-apply(res4_bias_ind1[,1:4], 2, FUN=function(x) as.character(x))
res4_bias_ind2[,1:4]<-apply(res4_bias_ind2[,1:4], 2, FUN=function(x) as.character(x))



#######################
### Wirte tex table ###
#######################



### Table 1

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 7), collapse=" & "), sep=" & " )

matout(res4_bias_dir1, file=paste(od, "res4_bias_dir1.tex", sep=""), digits=4, stars=T,
       caption="Bias of direct impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res4_dir1b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 2

matout(res4_bias_dir2, file=paste(od, "res4_bias_dir2.tex", sep=""), digits=4, stars=T,
       caption="Bias of direct impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res4_dir2b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")



### Table 3

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("SLX", "SAR", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 5), collapse=" & "), sep=" & " )

matout(res4_bias_ind1, file=paste(od, "res4_bias_ind1.tex", sep=""), digits=4, stars=T,
       caption="Bias of indirect impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res4_ind1b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 4

matout(res4_bias_ind2, file=paste(od, "res4_bias_ind2.tex", sep=""), digits=4, stars=T,
       caption="Bias of indirect impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res4_ind2b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")









##############################
#### Create combined plot ####
##############################


### SD Direct impacts ###

res4_sd_dir<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res4_sd_dir)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res4_sd_dir<-cbind(res4_sd_dir,tmp)

for(i in 1:nrow(res4_sd_dir)){
  tmp.sim<-get(paste("sim4_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sd[-seq(1, nrow(s$mat_d), by=3),1]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res4_sd_dir[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res4_sd_dir[i, 6:ncol(res4_sd_dir)]<-t(m)
}



### Separate into omv and non-omv

res4_sd_dir1<-res4_sd_dir[res4_sd_dir$gamma==0,]
res4_sd_dir2<-res4_sd_dir[res4_sd_dir$gamma>0,]

res4_sd_dir1<-res4_sd_dir1[,-which(names(res4_sd_dir1)=="gamma")]
res4_sd_dir2<-res4_sd_dir2[,-which(names(res4_sd_dir2)=="gamma")]


res4_sd_dir1[,1:4]<-apply(res4_sd_dir1[,1:4], 2, FUN=function(x) as.character(x))
res4_sd_dir2[,1:4]<-apply(res4_sd_dir2[,1:4], 2, FUN=function(x) as.character(x))





### SD Indirect impacts ###

res4_sd_ind<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res4_sd_ind)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res4_sd_ind<-cbind(res4_sd_ind,tmp)

for(i in 1:nrow(res4_sd_ind)){
  tmp.sim<-get(paste("sim4_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sd[-seq(1, nrow(s$mat_sd), by=3),2]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res4_sd_ind[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res4_sd_ind[i, 6:ncol(res4_sd_ind)]<-t(m)
}


res4_sd_ind<-res4_sd_ind[,-which(names(res4_sd_ind) %in% c("OLS_x1", "OLS_x2", "SEM_x1", "SEM_x2"))]

### Separate into omv and non-omv

res4_sd_ind1<-res4_sd_ind[res4_sd_ind$gamma==0,]
res4_sd_ind2<-res4_sd_ind[res4_sd_ind$gamma>0,]

res4_sd_ind1<-res4_sd_ind1[,-which(names(res4_sd_ind1)=="gamma")]
res4_sd_ind2<-res4_sd_ind2[,-which(names(res4_sd_ind2)=="gamma")]


res4_sd_ind1[,1:4]<-apply(res4_sd_ind1[,1:4], 2, FUN=function(x) as.character(x))
res4_sd_ind2[,1:4]<-apply(res4_sd_ind2[,1:4], 2, FUN=function(x) as.character(x))


### Expand strings
res4_bias_dir1$delta[which(res4_bias_dir1$delta=="0, 0")]<-c("0.0, 0.0")
res4_bias_ind1$delta[which(res4_bias_ind1$delta=="0, 0")]<-c("0.0, 0.0")
res4_bias_dir1$theta[which(res4_bias_dir1$theta=="0, 0")]<-c("0.0, 0.0")
res4_bias_ind1$theta[which(res4_bias_ind1$theta=="0, 0")]<-c("0.0, 0.0")

res4_bias_dir2$delta[which(res4_bias_dir2$delta=="0, 0")]<-c("0.0, 0.0")
res4_bias_ind2$delta[which(res4_bias_ind2$delta=="0, 0")]<-c("0.0, 0.0")
res4_bias_dir2$theta[which(res4_bias_dir2$theta=="0, 0")]<-c("0.0, 0.0")
res4_bias_ind2$theta[which(res4_bias_ind2$theta=="0, 0")]<-c("0.0, 0.0")


# res4_bias_dir1$delta<-sprintf("%-15s", res4_bias_dir1$delta)
# res4_bias_dir1$theta<-sprintf("%-15s", res4_bias_dir1$theta)
# res4_bias_ind1$delta<-sprintf("%-15s", res4_bias_ind1$delta)
# res4_bias_ind1$theta<-sprintf("%-15s", res4_bias_ind1$theta)




###########################
### Figure without omv  ###
###########################

tmp<-res4_bias_dir1[,grep("x1", colnames(res4_bias_dir1))]
tmp_sd<-res4_sd_dir1[,grep("x1", colnames(res4_sd_dir1))]


### Set up model Frame
mf_dir1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res4_bias_dir1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir1$Model<-rep(names,7)
mf_dir1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir1$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir1$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir1$Coefficient<-as.numeric(mf_dir1$Coefficient)
mf_dir1$SE<-as.numeric(mf_dir1$SE)


# Reorder levels
mf_dir1$Model <- factor(mf_dir1$Model, levels = rev(mf_dir1$Model[1:16]))
mf_dir1$modelName <- factor(mf_dir1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_all<-mf_dir1
mf_all$eff<-"Direct x1"


### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp1 <- ggplot(mf_dir1, aes(colour = modelName, shape = modelName))
zp1 <- zp1 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp1 <- zp1 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp1 <- zp1 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + ylim(xlim[1], xlim[2])
zp1 <- zp1 + theme(legend.title = element_blank())
zp1 <- zp1 + labs(y="Bias", x=element_blank()) 
zp1 <- zp1 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp1 <- zp1 + ggtitle(element_blank()) 
zp1 <- zp1 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp1)



### Figure without omv x2 ###

tmp<-res4_bias_dir1[,grep("x2", colnames(res4_bias_dir1))]
tmp_sd<-res4_sd_dir1[,grep("x2", colnames(res4_sd_dir1))]


### Set up model Frame
mf_dir1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res4_bias_dir1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir1$Model<-rep(names,7)
mf_dir1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir1$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir1$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir1$Coefficient<-as.numeric(mf_dir1$Coefficient)
mf_dir1$SE<-as.numeric(mf_dir1$SE)


# Reorder levels
mf_dir1$Model <- factor(mf_dir1$Model, levels = rev(mf_dir1$Model[1:16]))
mf_dir1$modelName <- factor(mf_dir1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_dir1$eff<-"Direct x2"
mf_all<-rbind(mf_all, mf_dir1)



### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp2 <- ggplot(mf_dir1, aes(colour = modelName, shape = modelName))
zp2 <- zp2 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp2 <- zp2 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp2 <- zp2 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp2 <- zp2 + coord_flip() + theme_bw()
zp2 <- zp2 + ylim(xlim[1], xlim[2])
zp2 <- zp2 + theme(legend.title = element_blank())
zp2 <- zp2 + labs(y="Bias", x=element_blank()) 
zp2 <- zp2 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp2 <- zp2 + ggtitle(element_blank()) 
zp2 <- zp2 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp2)




#### Indirect ####

### Figure without omv x1 ###

tmp<-res4_bias_ind1[,grep("x1", colnames(res4_bias_ind1))]
tmp_sd<-res4_sd_ind1[,grep("x1", colnames(res4_sd_ind1))]


### Set up model Frame
mf_ind1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res4_bias_ind1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind1$Model<-rep(names,7)
mf_ind1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind1$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind1$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind1$Coefficient<-as.numeric(mf_ind1$Coefficient)
mf_ind1$SE<-as.numeric(mf_ind1$SE)


# Reorder levels
mf_ind1$Model <- factor(mf_ind1$Model, levels = rev(mf_ind1$Model[1:16]))
mf_ind1$modelName <- factor(mf_ind1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind1$eff<-"Indirect x1"
mf_all<-rbind(mf_all, mf_ind1)


### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp3 <- ggplot(mf_ind1, aes(colour = modelName, shape = modelName))
zp3 <- zp3 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp3 <- zp3 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp3 <- zp3 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp3 <- zp3 + coord_flip() + theme_bw()
zp3 <- zp3 + ylim(xlim[1], xlim[2])
zp3 <- zp3 + theme(legend.title = element_blank())
zp3 <- zp3 + labs(y="Bias", x=element_blank()) 
zp3 <- zp3 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp3 <- zp3 + ggtitle(element_blank()) 
zp3 <- zp3 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp3)





### Figure without omv x2 ###

tmp<-res4_bias_ind1[,grep("x2", colnames(res4_bias_ind1))]
tmp_sd<-res4_sd_ind1[,grep("x2", colnames(res4_sd_ind1))]


### Set up model Frame
mf_ind1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res4_bias_ind1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind1$Model<-rep(names,7)
mf_ind1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind1$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind1$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind1$Coefficient<-as.numeric(mf_ind1$Coefficient)
mf_ind1$SE<-as.numeric(mf_ind1$SE)


# Reorder levels
mf_ind1$Model <- factor(mf_ind1$Model, levels = rev(mf_ind1$Model[1:16]))
mf_ind1$modelName <- factor(mf_ind1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind1$eff<-"Indirect x2"
mf_all<-rbind(mf_all, mf_ind1)

### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp4 <- ggplot(mf_ind1, aes(colour = modelName, shape = modelName))
zp4 <- zp4 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp4 <- zp4 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp4 <- zp4 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp4 <- zp4 + coord_flip() + theme_bw()
zp4 <- zp4 + ylim(xlim[1], xlim[2])
zp4 <- zp4 + theme(legend.title = element_blank())
zp4 <- zp4 + labs(y="Bias", x=element_blank()) 
zp4 <- zp4 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp4 <- zp4 + ggtitle(element_blank()) 
zp4 <- zp4 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp4)






#### Combine plots ####

mf_all$eff <- factor(mf_all$eff, levels = c("Direct x1", "Direct x2", "Indirect x1", "Indirect x2"))

mf_all$y_min<-NA
mf_all$y_min[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<--0.3
mf_all$y_min[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<--1.5

mf_all$y_max<-NA
mf_all$y_max[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<-0.7
mf_all$y_max[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<-2.2


### Plot

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

# Set limits (trim data manually)
xlim=c(-1.5, 2.2)

mf_all$lb<-mf_all$Coefficient - mf_all$SE*interval2
mf_all$ub<-mf_all$Coefficient + mf_all$SE*interval2

mf_all$ub[which(mf_all$ub>mf_all$y_max)]<-mf_all$y_max[which(mf_all$ub>mf_all$y_max)]
mf_all$lb[which(mf_all$lb<mf_all$y_min)]<-mf_all$y_min[which(mf_all$lb<mf_all$y_min)]



# Plot
zp_all <- ggplot(mf_all, aes(colour = modelName, shape = modelName))
zp_all <- zp_all + facet_grid(. ~ eff, scales="free_x")
zp_all <- zp_all + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp_all <- zp_all + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.7)
zp_all <- zp_all + geom_pointrange(aes(x = Model, y = Coefficient, ymin =lb,
                                       ymax = ub, 0.2),
                                   lwd = 0.7, position = position_dodge(width = 0.6),
                                   fill ="black")
zp_all <- zp_all + coord_flip() + theme_bw()
#zp_all <- zp_all + ylim(xlim[1], xlim[2])
zp_all <- zp_all + geom_blank(aes(y = y_min)) + geom_blank(aes(y = y_max))
zp_all <- zp_all + scale_y_continuous(expand = c(0,0) )
zp_all <- zp_all + theme(legend.title = element_blank())
zp_all <- zp_all + labs(y="Bias", x=element_blank()) 
zp_all <- zp_all + theme(legend.text = element_text(size = 20),
                         legend.position="bottom",
                         legend.key=element_blank(),
                         axis.text.x = element_text(size=16),
                         axis.text.y = element_text(size=16, colour="black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20),
                         plot.title = element_text(size=20, margin = margin(t = 10, b = -19 ), hjust=-0.235),
                         strip.background =element_blank(),
                         strip.text = element_text(size=20, colour="black"),
                         panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp_all <- zp_all + ggtitle(bquote(rho ~~~~~~ bold(delta) ~~~~~~~~~~~ lambda ~~~~~ bold(theta))) 
zp_all <- zp_all + guides(fill=guide_legend(ncol = 1), colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                          shape = guide_legend(reverse = T))
print(zp_all)




#### Combine plots ####


cairo_pdf(file=paste(od, "sim4_combined_1.pdf", sep=""), width=17.54, height=12.40, 
          bg = "white", family="CM Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
zp_all
dev.off()








########################
### Figure with omv  ###
########################

tmp<-res4_bias_dir2[,grep("x1", colnames(res4_bias_dir2))]
tmp_sd<-res4_sd_dir2[,grep("x1", colnames(res4_sd_dir2))]


### Set up model Frame
mf_dir2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res4_bias_dir2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir2$Model<-rep(names,7)
mf_dir2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir2$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir2$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir2$Coefficient<-as.numeric(mf_dir2$Coefficient)
mf_dir2$SE<-as.numeric(mf_dir2$SE)


# Reorder levels
mf_dir2$Model <- factor(mf_dir2$Model, levels = rev(mf_dir2$Model[1:16]))
mf_dir2$modelName <- factor(mf_dir2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_all<-mf_dir2
mf_all$eff<-"Direct x1"


### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp1 <- ggplot(mf_dir2, aes(colour = modelName, shape = modelName))
zp1 <- zp1 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp1 <- zp1 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp1 <- zp1 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + ylim(xlim[1], xlim[2])
zp1 <- zp1 + theme(legend.title = element_blank())
zp1 <- zp1 + labs(y="Bias", x=element_blank()) 
zp1 <- zp1 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp1 <- zp1 + ggtitle(element_blank()) 
zp1 <- zp1 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp1)





### Figure without omv x2 ###

tmp<-res4_bias_dir2[,grep("x2", colnames(res4_bias_dir2))]
tmp_sd<-res4_sd_dir2[,grep("x2", colnames(res4_sd_dir2))]


### Set up model Frame
mf_dir2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res4_bias_dir2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir2$Model<-rep(names,7)
mf_dir2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir2$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir2$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir2$Coefficient<-as.numeric(mf_dir2$Coefficient)
mf_dir2$SE<-as.numeric(mf_dir2$SE)


# Reorder levels
mf_dir2$Model <- factor(mf_dir2$Model, levels = rev(mf_dir2$Model[1:16]))
mf_dir2$modelName <- factor(mf_dir2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_dir2$eff<-"Direct x2"
mf_all<-rbind(mf_all, mf_dir2)



### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp2 <- ggplot(mf_dir2, aes(colour = modelName, shape = modelName))
zp2 <- zp2 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp2 <- zp2 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp2 <- zp2 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp2 <- zp2 + coord_flip() + theme_bw()
zp2 <- zp2 + ylim(xlim[1], xlim[2])
zp2 <- zp2 + theme(legend.title = element_blank())
zp2 <- zp2 + labs(y="Bias", x=element_blank()) 
zp2 <- zp2 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp2 <- zp2 + ggtitle(element_blank()) 
zp2 <- zp2 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp2)




#### Indirect ####

### Figure without omv x1 ###

tmp<-res4_bias_ind2[,grep("x1", colnames(res4_bias_ind2))]
tmp_sd<-res4_sd_ind2[,grep("x1", colnames(res4_sd_ind2))]


### Set up model Frame
mf_ind2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res4_bias_ind2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind2$Model<-rep(names,7)
mf_ind2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind2$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind2$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind2$Coefficient<-as.numeric(mf_ind2$Coefficient)
mf_ind2$SE<-as.numeric(mf_ind2$SE)


# Reorder levels
mf_ind2$Model <- factor(mf_ind2$Model, levels = rev(mf_ind2$Model[1:16]))
mf_ind2$modelName <- factor(mf_ind2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind2$eff<-"Indirect x1"
mf_all<-rbind(mf_all, mf_ind2)


### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp3 <- ggplot(mf_ind2, aes(colour = modelName, shape = modelName))
zp3 <- zp3 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp3 <- zp3 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp3 <- zp3 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp3 <- zp3 + coord_flip() + theme_bw()
zp3 <- zp3 + ylim(xlim[1], xlim[2])
zp3 <- zp3 + theme(legend.title = element_blank())
zp3 <- zp3 + labs(y="Bias", x=element_blank()) 
zp3 <- zp3 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp3 <- zp3 + ggtitle(element_blank()) 
zp3 <- zp3 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp3)





### Figure without omv x2 ###

tmp<-res4_bias_ind2[,grep("x2", colnames(res4_bias_ind2))]
tmp_sd<-res4_sd_ind2[,grep("x2", colnames(res4_sd_ind2))]


### Set up model Frame
mf_ind2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res4_bias_ind2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind2$Model<-rep(names,7)
mf_ind2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind2$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind2$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind2$Coefficient<-as.numeric(mf_ind2$Coefficient)
mf_ind2$SE<-as.numeric(mf_ind2$SE)


# Reorder levels
mf_ind2$Model <- factor(mf_ind2$Model, levels = rev(mf_ind2$Model[1:16]))
mf_ind2$modelName <- factor(mf_ind2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind2$eff<-"Indirect x2"
mf_all<-rbind(mf_all, mf_ind2)

### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp4 <- ggplot(mf_ind2, aes(colour = modelName, shape = modelName))
zp4 <- zp4 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp4 <- zp4 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp4 <- zp4 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp4 <- zp4 + coord_flip() + theme_bw()
zp4 <- zp4 + ylim(xlim[1], xlim[2])
zp4 <- zp4 + theme(legend.title = element_blank())
zp4 <- zp4 + labs(y="Bias", x=element_blank()) 
zp4 <- zp4 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp4 <- zp4 + ggtitle(element_blank()) 
zp4 <- zp4 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp4)






#### Combine plots ####

mf_all$eff <- factor(mf_all$eff, levels = c("Direct x1", "Direct x2", "Indirect x1", "Indirect x2"))

mf_all$y_min<-NA
mf_all$y_min[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<--0.3
mf_all$y_min[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<--1.5

mf_all$y_max<-NA
mf_all$y_max[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<-0.7
mf_all$y_max[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<-2.2


### Plot

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

# Set limits (trim data manually)
xlim=c(-1.5, 2.2)

mf_all$lb<-mf_all$Coefficient - mf_all$SE*interval2
mf_all$ub<-mf_all$Coefficient + mf_all$SE*interval2

mf_all$ub[which(mf_all$ub>mf_all$y_max)]<-mf_all$y_max[which(mf_all$ub>mf_all$y_max)]
mf_all$lb[which(mf_all$lb<mf_all$y_min)]<-mf_all$y_min[which(mf_all$lb<mf_all$y_min)]



# Plot
zp_all <- ggplot(mf_all, aes(colour = modelName, shape = modelName))
zp_all <- zp_all + facet_grid(. ~ eff, scales="free_x")
zp_all <- zp_all + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp_all <- zp_all + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.7)
zp_all <- zp_all + geom_pointrange(aes(x = Model, y = Coefficient, ymin =lb,
                                       ymax = ub, 0.2),
                                   lwd = 0.7, position = position_dodge(width = 0.6),
                                   fill ="black")
zp_all <- zp_all + coord_flip() + theme_bw()
#zp_all <- zp_all + ylim(xlim[1], xlim[2])
zp_all <- zp_all + geom_blank(aes(y = y_min)) + geom_blank(aes(y = y_max))
zp_all <- zp_all + scale_y_continuous(expand = c(0,0) )
zp_all <- zp_all + theme(legend.title = element_blank())
zp_all <- zp_all + labs(y="Bias", x=element_blank()) 
zp_all <- zp_all + theme(legend.text = element_text(size = 20),
                         legend.position="bottom",
                         legend.key=element_blank(),
                         axis.text.x = element_text(size=16),
                         axis.text.y = element_text(size=16, colour="black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20),
                         plot.title = element_text(size=20, margin = margin(t = 10, b = -19 ), hjust=-0.235),
                         strip.background =element_blank(),
                         strip.text = element_text(size=20, colour="black"),
                         panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp_all <- zp_all + ggtitle(bquote(rho ~~~~~~ bold(delta) ~~~~~~~~~~~ lambda ~~~~~ bold(theta))) 
zp_all <- zp_all + guides(fill=guide_legend(ncol = 1), colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                          shape = guide_legend(reverse = T))
print(zp_all)




#### Combine plots ####


cairo_pdf(file=paste(od, "sim4_combined_2.pdf", sep=""), width=17.54, height=12.40, 
          bg = "white", family="CM Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
zp_all
dev.off()






######################################
######################################
#### 2) Specification Test        ####
######################################
######################################

### Load program
source(paste(dd, "01g_Monte Carlo Simulation Spatial_Specification Test_Program.R", sep="/"))


######################################
#### Define possible combinations ####
######################################

beta<-c(0.2, 0.5)

rho<-c(0, 0.4, 0.8)
delta<-c("0, 0", "0.4, 0.7")
lambda<-c(0, 0.4, 0.8)
gamma<-c(0)
theta<-c("0, 0", "0.1, 0.8")


### Set up matrix of all combinations ###
# (Do once and save / load to use same matrix when restarting simulations)

nc<-2^5

l<-list(rho, delta, lambda, gamma, theta)

comb_mat_test1<-expand.grid(l)

colnames(comb_mat_test1)<-c("rho", "delta", "lambda", "gamma", "theta")
comb_mat_test1$delta<-as.character(comb_mat_test1$delta)
comb_mat_test1$theta<-as.character(comb_mat_test1$theta)

save(comb_mat_test1, file="Combinations_test1.RData")




###########################################
#### Start simulations of combinations ####
###########################################

# Set up manual loops (to spare capacity)


k<-c(1:36)

### Simumlation Run
for(i in k){

  simtest1<-NA # avoid using the last sim (in case of error)

  par<-comb_mat_test1[i,]
  deltai<-as.numeric(c(unlist(strsplit(par$delta, split=", "))))
  thetai<-as.numeric(c(unlist(strsplit(par$theta, split=", "))))

  simtest1<-spsim_lmtest(N=900, R=1000, seed=1573262181,
              ay=par$rho, ax1=deltai[1], ax2=deltai[2], au=par$lambda,
              bx1=beta[1], bwx1=thetai[1], bx2=beta[2], bwx2=thetai[2],
              bu=1, bwu=0, sdu=1,
              omv1=par$gamma, womv1=0)

  assign(paste0("simtest1_", i, sep=""), simtest1)



  save(list=c(paste0("simtest1_", i, sep="")), file=paste0("Simulation_test1_", i, ".RData", sep=""))

}


# In case you reload results:
# for(i in k){
#   load(paste("Simulation_test1_", i, ".RData", sep=""))
# }





##########################################
#### Create result matrices rejection ####
##########################################


### Rejection rates ###

res_test1<-comb_mat_test1

tmp<-data.frame(matrix(NA, ncol=5, nrow=36))
names(tmp)<-c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA")

res_test1<-cbind(res_test1,tmp)

for(i in 1:nrow(res_test1)){
  tmp.sim<-get(paste("simtest1_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary(tmp.sim)
  m<-s$mat[2,]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res_test1[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res_test1[i, 6:ncol(res_test1)]<-m
}



#######################
### Wirte tex table ###
#######################


matout2<-function(x, file="", digits=NULL, threepart=T, 
                  code_before=NULL, code_after=NULL,
                  header1=NULL, header2=NULL, 
                  caption=NULL, label=NULL, note=NULL,...) {
  
  ### Options for xtable
  default_args = list(include.colnames=FALSE, only.contents=TRUE,
                      include.rownames=FALSE, hline.after=NULL, comment=FALSE,
                      print.results=FALSE, sanitize.text.function = function(z) z)
  passed_args = list(...)
  
  
  if(is.null(caption)){caption<-"Table"}
  
  c<-ncol(x)
  
  if(is.null(code_before)){
    code_before<-paste("\\begin{tabular}{", paste(rep("l", c), collapse = " "), "}")
  }
  if(is.null(code_after)){
    code_after<-paste("\\end{tabular}")
  }
  
  num<-which(sapply(x,FUN=function(x) is.numeric(x)))
  
  
  calling_args = c(list(x=xtable(x, digits=digits)),
                   c(passed_args,
                     default_args[setdiff(names(default_args), names(passed_args))]))
  
  # midrule<-paste(seq(num[1], max(num), 2), "-", seq(num[1]+1, max(num), 2), sep="")
  # midrule<-c(paste("\\cmidrule(r){", midrule[1], "}", sep=""), 
  #            paste("\\cmidrule(lr){", midrule[-c(1, length(midrule))], "}", sep=""),
  #            paste("\\cmidrule(l){", midrule[length(midrule)], "}", sep=""))
  
  
  ### Output
  cat("\\begin{table}[ht!]\n", file=file)
  cat("\\centering\n", file=file, append=T)
  
  if(threepart==T){
    cat("{\\begin{threeparttable}\n", file=file, append=T)
  }
  
  cat(paste0("\\caption{", caption, "}\n"), file=file, append=T)
  
  if(!is.null(label)){
    cat(paste0("\\label{", label, "}\n"), file=file, append=T)
  }
  
  cat(paste0(code_before, "\n"), file=file, append=T)
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  
  if(!is.null(header1)){
    cat(paste(header1), "\\\\ \n", file=file, append=T)
  }else{cat(paste(paste("M", seq(1, c,), collapse=" & "), "\\\\ \n"), file=file, append=T)}
  
  # cat(midrule, file=file, append=T)
  
  if(!is.null(header2)){
    cat(paste(header2), "\\\\ \n", file=file, append=T)
  }
  
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  cat(do.call(print.xtable, calling_args), file=file, append=T)
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  
  
  cat(paste0(code_after, "\n"), file=file, append=T)
  
  if(!is.null(note)){
    cat("\\begin{tablenotes}\n", file=file, append=T)
    cat(paste("\\item \\scriptsize{", note, "}\n", collape=""), file=file, append=T)
    cat("\\end{tablenotes}\n", file=file, append=T)
  }
  
  if(threepart==T){
    cat("\\end{threeparttable}}\n", file=file, append=T)
  }
  
  cat("\\end{table}\n", file=file, append=T)
}



### Table 1


res_test1<-res_test1[,-which(colnames(res_test1)=="gamma")]

res_test1[,1:4]<-apply(res_test1[,1:4], 2, FUN=function(x) as.character(x))


header1<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "),
               paste(paste("\\tc{", c("$LM_\\lambda$","$LM_\\rho$","$LM^*_\\lambda$", "$LM^*_\\rho$", "$LM_{\\lambda\\rho}$"), "}", sep =""), collapse=" & "), sep=" & ")

matout2(res_test1, file=paste(od, "res_test1.tex", sep=""), digits=4,
        caption="Rejection rates of $H_0$ (Lagrange multiplier test). ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
        label="table:res_test1",
        code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
        header1=header1, 
        note="Number of observations=900, repetitions=1000. $LM$= Lagrange multiplier test, $LM^*$= Robust Lagrange multiplier test, each for $H_0$: $\\lambda=0$, $H_0$: $\\rho=0$, $H_0$: $\\lambda,\\rho=0$.")




######################################
######################################
#### 3) Dist NB Weights           ####
######################################
######################################

source(paste(dd, "01d_Monte Carlo Simulation Spatial_Program_distnb.R", sep="/"))


######################################
#### Define possible combinations ####
######################################

beta<-c(0.2, 0.5)

rho<-c(0, 0.5)
delta<-c("0, 0", "0.4, 0.7")
lambda<-c(0, 0.5)
gamma<-c(0, 0.3)
theta<-c("0, 0", "0.1, 0.8")


### Set up matrix of all combinations ###

nc<-2^5

l<-list(rho, delta, lambda, gamma, theta)

comb_mat4<-expand.grid(l)

colnames(comb_mat4)<-c("rho", "delta", "lambda", "gamma", "theta")
comb_mat4$delta<-as.character(comb_mat4$delta)
comb_mat4$theta<-as.character(comb_mat4$theta)



###########################################
#### Start simulations of combinations ####
###########################################

# Set up manual loops (to spare capacity)


k<-c(1:32)

### Simumlation Run
for(i in k){

  sim15<-NA # avoid using the last sim (in case of error)

  par<-comb_mat4[i,]
  deltai<-as.numeric(c(unlist(strsplit(par$delta, split=", "))))
  thetai<-as.numeric(c(unlist(strsplit(par$theta, split=", "))))
  sim15<-spsim(N=900, R=1000, seed=1573262181,
              ay=par$rho, ax1=deltai[1], ax2=deltai[2], au=par$lambda,
              bx1=beta[1], bwx1=thetai[1], bx2=beta[2], bwx2=thetai[2],
              bu=1, bwu=0, sdu=1,
              omv1=par$gamma, womv1=0)
  assign(paste0("sim15_", i, sep=""), sim15)



  save(list=c(paste0("sim15_", i, sep="")), file=paste0("Simulation15_", i, ".RData", sep=""))

}




for(i in k){
  load(paste("Simulation15_", i, ".RData", sep=""))
}





#######################################
#### Make plots of each Simulation ####
#######################################


### Jpeg ###

for(i in k){
  cm<-comb_mat4[i,]
  
  tmp.sim<-get(paste("sim15_", i, sep=""))
  
  p<-plot.sim(tmp.sim, xlim=c(-2, 3))
  
  name = paste(od, "sim15_", i, ".jpeg", sep="")
  jpeg(file = name, bg = "transparent", width = 2720, height = 1820, quality=100, family="CM Roman", res=300)
  par(mar=c(0,0,0,0))
  par(mfrow=c(1,1),oma=c(0,0,0,0))
  
  print(p + ggtitle(bquote(rho == .(cm$rho) ~~
                             delta == .(cm$delta) ~~
                             lambda == .(cm$lambda) ~~
                             gamma == .(cm$gamma) ~~
                             theta == .(cm$theta))))
  dev.off()
  par(mfrow=c(1,1))
}





#####################################
#### Create result matrices RMSE ####
#####################################


### Direct impacts ###

res15_rmse_dir<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res15_rmse_dir)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res15_rmse_dir<-cbind(res15_rmse_dir,tmp)

for(i in 1:nrow(res15_rmse_dir)){
  tmp.sim<-get(paste("sim15_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_d[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res15_rmse_dir[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res15_rmse_dir[i, 6:ncol(res15_rmse_dir)]<-t(m[,2])
}



### Separate into omv and non-omv

res15_rmse_dir1<-res15_rmse_dir[res15_rmse_dir$gamma==0,]
res15_rmse_dir2<-res15_rmse_dir[res15_rmse_dir$gamma>0,]

res15_rmse_dir1<-res15_rmse_dir1[,-which(names(res15_rmse_dir1)=="gamma")]
res15_rmse_dir2<-res15_rmse_dir2[,-which(names(res15_rmse_dir2)=="gamma")]


res15_rmse_dir1[,1:4]<-apply(res15_rmse_dir1[,1:4], 2, FUN=function(x) as.character(x))
res15_rmse_dir2[,1:4]<-apply(res15_rmse_dir2[,1:4], 2, FUN=function(x) as.character(x))





### Indirect impacts ###

res15_rmse_ind<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res15_rmse_ind)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res15_rmse_ind<-cbind(res15_rmse_ind,tmp)

for(i in 1:nrow(res15_rmse_ind)){
  tmp.sim<-get(paste("sim15_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sp[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res15_rmse_ind[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res15_rmse_ind[i, 6:ncol(res15_rmse_ind)]<-t(m[,2])
}


res15_rmse_ind<-res15_rmse_ind[,-which(names(res15_rmse_ind) %in% c("OLS_x1", "OLS_x2", "SEM_x1", "SEM_x2"))]

### Separate into omv and non-omv

res15_rmse_ind1<-res15_rmse_ind[res15_rmse_ind$gamma==0,]
res15_rmse_ind2<-res15_rmse_ind[res15_rmse_ind$gamma>0,]

res15_rmse_ind1<-res15_rmse_ind1[,-which(names(res15_rmse_ind1)=="gamma")]
res15_rmse_ind2<-res15_rmse_ind2[,-which(names(res15_rmse_ind2)=="gamma")]


res15_rmse_ind1[,1:4]<-apply(res15_rmse_ind1[,1:4], 2, FUN=function(x) as.character(x))
res15_rmse_ind2[,1:4]<-apply(res15_rmse_ind2[,1:4], 2, FUN=function(x) as.character(x))



#######################
### Wirte tex table ###
#######################


matout<-function(x, file="", digits=NULL, colav=T, threepart=T, stars=F,
                 code_before=NULL, code_after=NULL,
                 header1=NULL, header2=NULL, 
                 caption=NULL, label=NULL, note=NULL,...) {
  
  ### Options for xtable
  default_args = list(include.colnames=FALSE, only.contents=TRUE,
                      include.rownames=FALSE, hline.after=NULL, comment=FALSE,
                      print.results=FALSE, sanitize.text.function = function(z) z)
  passed_args = list(...)
  
  
  if(is.null(caption)){caption<-"Table"}
  
  c<-ncol(x)
  
  if(is.null(code_before)){
    code_before<-paste("\\begin{tabular}{", paste(rep("l", c), collapse = " "), "}")
  }
  if(is.null(code_after)){
    code_after<-paste("\\end{tabular}")
  }
  
  num<-which(sapply(x,FUN=function(x) is.numeric(x)))
  
  if(colav==T){
    av<-data.frame(t(colMeans(abs(x[,num[1]:ncol(x)]))))
    av2<-data.frame(t(colMeans(matrix(t(colMeans(abs(x[,num[1]:ncol(x)]))), 2))))
    
    if(stars){
      pa1<-colnames(av[seq(1, ncol(av), 2)])[which.min(abs(av[seq(1, ncol(av), 2)]))]
      pa2<-colnames(av[seq(2, ncol(av), 2)])[which.min(abs(av[seq(2, ncol(av), 2)]))]
      pa3<-colnames(av2)[which.min(av2)]
      
      #av[1,]<-round(av[1,], digits)
      av[1,]<-sprintf(paste("%.", digits, "f"), round(av[1,],digits))
      av2[1,]<-sprintf(paste("%.", digits, "f"), round(av2[1,],digits))
      
      
      av[1,colnames(av)==pa1]<-paste(av[1,colnames(av)==pa1], "\\ensuremath{^{*}}", sep="")
      av[1,colnames(av)==pa2]<-paste(av[1,colnames(av)==pa2], "\\ensuremath{^{*}}", sep="")
      av2[1,colnames(av2)==pa3]<-paste(av2[1,colnames(av2)==pa3], "\\ensuremath{^{*}}", sep="")
      
      
    }
    
    calling_args2 = c(list(x=xtable(av, digits=digits)),
                      c(passed_args,
                        default_args[setdiff(names(default_args), names(passed_args))]))
    
  }
  
  if(stars){
    p11<-apply(x[,seq(num[1], ncol(x), 2)], 1, FUN=function(z) which.min(abs(z)))
    p11<-colnames(x[,seq(num[1], ncol(x), 2)])[p11]
    
    p12<-apply(x[,seq(num[1]+1, ncol(x), 2)], 1, FUN=function(z) which.min(abs(z)))
    p12<-colnames(x[,seq(num[1]+1, ncol(x), 2)])[p12]
    
    #x[, num]<-round(x[, num], digits)
    #x[, num]<-apply(x[, num], 2, FUN=function(z) as.character(format(z, digits=digits, scientific=F)))
    x[, num]<-apply(x[, num], 2, FUN=function(z) sprintf(paste("%.", digits, "f"), round(z,digits)))
    
    
    for(i in 1:nrow(x)){
      x[i,colnames(x)==p11[i]]<-paste(x[i,colnames(x)==p11[i]], "\\ensuremath{^{*}}", sep="")
      x[i,colnames(x)==p12[i]]<-paste(x[i,colnames(x)==p12[i]], "\\ensuremath{^{*}}", sep="")
    }
    
    
  }
  
  calling_args = c(list(x=xtable(x, digits=digits)),
                   c(passed_args,
                     default_args[setdiff(names(default_args), names(passed_args))]))
  
  midrule<-paste(seq(num[1], max(num), 2), "-", seq(num[1]+1, max(num), 2), sep="")
  midrule<-c(paste("\\cmidrule(r){", midrule[1], "}", sep=""), 
             paste("\\cmidrule(lr){", midrule[-c(1, length(midrule))], "}", sep=""),
             paste("\\cmidrule(l){", midrule[length(midrule)], "}", sep=""))
  
  
  ### Output
  cat("\\begin{table}\n", file=file)
  cat("\\centering\n", file=file, append=T)
  
  if(threepart==T){
    cat("{\\begin{threeparttable}\n", file=file, append=T)
  }
  
  cat(paste0("\\caption{", caption, "}\n"), file=file, append=T)
  
  if(!is.null(label)){
    cat(paste0("\\label{", label, "}\n"), file=file, append=T)
  }
  
  cat(paste0(code_before, "\n"), file=file, append=T)
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  
  if(!is.null(header1)){
    cat(paste(header1), "\\\\ \n", file=file, append=T)
  }else{cat(paste(paste("M", seq(1, c,), collapse=" & "), "\\\\ \n"), file=file, append=T)}
  
  cat(midrule, file=file, append=T)
  
  if(!is.null(header2)){
    cat(paste(header2), "\\\\ \n", file=file, append=T)
  }
  
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  cat(do.call(print.xtable, calling_args), file=file, append=T)
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  
  if(colav==T){
    cat(paste("\\multicolumn{", (num[1]-1), "}{l}{Average} &", collapse=""), file=file, append=T)
    cat(do.call(print.xtable, calling_args2), file=file, append=T)
    cat(paste("\\multicolumn{", (num[1]-1), "}{l}{} &", collapse=""), file=file, append=T)
    cat(paste("\\multicolumn{2}{c}{", av2, "}", collapse=" & "), "\\\\ \n", file=file, append=T)
    cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  }
  
  cat(paste0(code_after, "\n"), file=file, append=T)
  
  if(!is.null(note)){
    cat("\\begin{tablenotes}\n", file=file, append=T)
    cat(paste("\\item \\scriptsize{", note, "}\n", collape=""), file=file, append=T)
    cat("\\end{tablenotes}\n", file=file, append=T)
  }
  
  if(threepart==T){
    cat("\\end{threeparttable}}\n", file=file, append=T)
  }
  
  cat("\\end{table}\n", file=file, append=T)
}



### Table 1

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 7), collapse=" & "), sep=" & " )

matout(res15_rmse_dir1, file=paste(od, "res15_rmse_dir1.tex", sep=""), digits=4, stars=T,
       caption="RMSE of direct impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res15_dir1",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 2

matout(res15_rmse_dir2, file=paste(od, "res15_rmse_dir2.tex", sep=""), digits=4, stars=T,
       caption="RMSE of direct impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res15_dir2",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")



### Table 3

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("SLX", "SAR", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 5), collapse=" & "), sep=" & " )

matout(res15_rmse_ind1, file=paste(od, "res15_rmse_ind1.tex", sep=""), digits=4, stars=T,
       caption="RMSE of indirect impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res15_ind1",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 4

matout(res15_rmse_ind2, file=paste(od, "res15_rmse_ind2.tex", sep=""), digits=4, stars=T,
       caption="RMSE of indirect impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res15_ind2",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")




#####################################
#### Create result matrices Bias ####
#####################################


### Direct impacts ###

res15_bias_dir<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res15_bias_dir)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res15_bias_dir<-cbind(res15_bias_dir,tmp)

for(i in 1:nrow(res15_bias_dir)){
  tmp.sim<-get(paste("sim15_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_d[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res15_bias_dir[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res15_bias_dir[i, 6:ncol(res15_bias_dir)]<-t(m[,1])
}



### Separate into omv and non-omv

res15_bias_dir1<-res15_bias_dir[res15_bias_dir$gamma==0,]
res15_bias_dir2<-res15_bias_dir[res15_bias_dir$gamma>0,]

res15_bias_dir1<-res15_bias_dir1[,-which(names(res15_bias_dir1)=="gamma")]
res15_bias_dir2<-res15_bias_dir2[,-which(names(res15_bias_dir2)=="gamma")]


res15_bias_dir1[,1:4]<-apply(res15_bias_dir1[,1:4], 2, FUN=function(x) as.character(x))
res15_bias_dir2[,1:4]<-apply(res15_bias_dir2[,1:4], 2, FUN=function(x) as.character(x))





### Indirect impacts ###

res15_bias_ind<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res15_bias_ind)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res15_bias_ind<-cbind(res15_bias_ind,tmp)

for(i in 1:nrow(res15_bias_ind)){
  tmp.sim<-get(paste("sim15_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sp[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res15_bias_ind[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res15_bias_ind[i, 6:ncol(res15_bias_ind)]<-t(m[,1])
}


res15_bias_ind<-res15_bias_ind[,-which(names(res15_bias_ind) %in% c("OLS_x1", "OLS_x2", "SEM_x1", "SEM_x2"))]

### Separate into omv and non-omv

res15_bias_ind1<-res15_bias_ind[res15_bias_ind$gamma==0,]
res15_bias_ind2<-res15_bias_ind[res15_bias_ind$gamma>0,]

res15_bias_ind1<-res15_bias_ind1[,-which(names(res15_bias_ind1)=="gamma")]
res15_bias_ind2<-res15_bias_ind2[,-which(names(res15_bias_ind2)=="gamma")]


res15_bias_ind1[,1:4]<-apply(res15_bias_ind1[,1:4], 2, FUN=function(x) as.character(x))
res15_bias_ind2[,1:4]<-apply(res15_bias_ind2[,1:4], 2, FUN=function(x) as.character(x))



#######################
### Wirte tex table ###
#######################



### Table 1

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 7), collapse=" & "), sep=" & " )

matout(res15_bias_dir1, file=paste(od, "res15_bias_dir1.tex", sep=""), digits=4, stars=T,
       caption="Bias of direct impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res15_dir1b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 2

matout(res15_bias_dir2, file=paste(od, "res15_bias_dir2.tex", sep=""), digits=4, stars=T,
       caption="Bias of direct impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res15_dir2b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")



### Table 3

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("SLX", "SAR", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 5), collapse=" & "), sep=" & " )

matout(res15_bias_ind1, file=paste(od, "res15_bias_ind1.tex", sep=""), digits=4, stars=T,
       caption="Bias of indirect impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res15_ind1b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 4

matout(res15_bias_ind2, file=paste(od, "res15_bias_ind2.tex", sep=""), digits=4, stars=T,
       caption="Bias of indirect impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res15_ind2b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")




##############################
#### Create combined plot ####
##############################


### SD Direct impacts ###

res15_sd_dir<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res15_sd_dir)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res15_sd_dir<-cbind(res15_sd_dir,tmp)

for(i in 1:nrow(res15_sd_dir)){
  tmp.sim<-get(paste("sim15_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sd[-seq(1, nrow(s$mat_d), by=3),1]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res15_sd_dir[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res15_sd_dir[i, 6:ncol(res15_sd_dir)]<-t(m)
}



### Separate into omv and non-omv

res15_sd_dir1<-res15_sd_dir[res15_sd_dir$gamma==0,]
res15_sd_dir2<-res15_sd_dir[res15_sd_dir$gamma>0,]

res15_sd_dir1<-res15_sd_dir1[,-which(names(res15_sd_dir1)=="gamma")]
res15_sd_dir2<-res15_sd_dir2[,-which(names(res15_sd_dir2)=="gamma")]


res15_sd_dir1[,1:4]<-apply(res15_sd_dir1[,1:4], 2, FUN=function(x) as.character(x))
res15_sd_dir2[,1:4]<-apply(res15_sd_dir2[,1:4], 2, FUN=function(x) as.character(x))





### SD Indirect impacts ###

res15_sd_ind<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res15_sd_ind)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res15_sd_ind<-cbind(res15_sd_ind,tmp)

for(i in 1:nrow(res15_sd_ind)){
  tmp.sim<-get(paste("sim15_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sd[-seq(1, nrow(s$mat_sd), by=3),2]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res15_sd_ind[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res15_sd_ind[i, 6:ncol(res15_sd_ind)]<-t(m)
}


res15_sd_ind<-res15_sd_ind[,-which(names(res15_sd_ind) %in% c("OLS_x1", "OLS_x2", "SEM_x1", "SEM_x2"))]

### Separate into omv and non-omv

res15_sd_ind1<-res15_sd_ind[res15_sd_ind$gamma==0,]
res15_sd_ind2<-res15_sd_ind[res15_sd_ind$gamma>0,]

res15_sd_ind1<-res15_sd_ind1[,-which(names(res15_sd_ind1)=="gamma")]
res15_sd_ind2<-res15_sd_ind2[,-which(names(res15_sd_ind2)=="gamma")]


res15_sd_ind1[,1:4]<-apply(res15_sd_ind1[,1:4], 2, FUN=function(x) as.character(x))
res15_sd_ind2[,1:4]<-apply(res15_sd_ind2[,1:4], 2, FUN=function(x) as.character(x))


### Expand strings
res15_bias_dir1$delta[which(res15_bias_dir1$delta=="0, 0")]<-c("0.0, 0.0")
res15_bias_ind1$delta[which(res15_bias_ind1$delta=="0, 0")]<-c("0.0, 0.0")
res15_bias_dir1$theta[which(res15_bias_dir1$theta=="0, 0")]<-c("0.0, 0.0")
res15_bias_ind1$theta[which(res15_bias_ind1$theta=="0, 0")]<-c("0.0, 0.0")

res15_bias_dir2$delta[which(res15_bias_dir2$delta=="0, 0")]<-c("0.0, 0.0")
res15_bias_ind2$delta[which(res15_bias_ind2$delta=="0, 0")]<-c("0.0, 0.0")
res15_bias_dir2$theta[which(res15_bias_dir2$theta=="0, 0")]<-c("0.0, 0.0")
res15_bias_ind2$theta[which(res15_bias_ind2$theta=="0, 0")]<-c("0.0, 0.0")


# res15_bias_dir1$delta<-sprintf("%-15s", res15_bias_dir1$delta)
# res15_bias_dir1$theta<-sprintf("%-15s", res15_bias_dir1$theta)
# res15_bias_ind1$delta<-sprintf("%-15s", res15_bias_ind1$delta)
# res15_bias_ind1$theta<-sprintf("%-15s", res15_bias_ind1$theta)




###########################
### Figure without omv  ###
###########################

tmp<-res15_bias_dir1[,grep("x1", colnames(res15_bias_dir1))]
tmp_sd<-res15_sd_dir1[,grep("x1", colnames(res15_sd_dir1))]


### Set up model Frame
mf_dir1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res15_bias_dir1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir1$Model<-rep(names,7)
mf_dir1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir1$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir1$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir1$Coefficient<-as.numeric(mf_dir1$Coefficient)
mf_dir1$SE<-as.numeric(mf_dir1$SE)


# Reorder levels
mf_dir1$Model <- factor(mf_dir1$Model, levels = rev(mf_dir1$Model[1:16]))
mf_dir1$modelName <- factor(mf_dir1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_all<-mf_dir1
mf_all$eff<-"Direct x1"


### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp1 <- ggplot(mf_dir1, aes(colour = modelName, shape = modelName))
zp1 <- zp1 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp1 <- zp1 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp1 <- zp1 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + ylim(xlim[1], xlim[2])
zp1 <- zp1 + theme(legend.title = element_blank())
zp1 <- zp1 + labs(y="Bias", x=element_blank()) 
zp1 <- zp1 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp1 <- zp1 + ggtitle(element_blank()) 
zp1 <- zp1 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp1)



### Figure without omv x2 ###

tmp<-res15_bias_dir1[,grep("x2", colnames(res15_bias_dir1))]
tmp_sd<-res15_sd_dir1[,grep("x2", colnames(res15_sd_dir1))]


### Set up model Frame
mf_dir1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res15_bias_dir1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir1$Model<-rep(names,7)
mf_dir1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir1$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir1$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir1$Coefficient<-as.numeric(mf_dir1$Coefficient)
mf_dir1$SE<-as.numeric(mf_dir1$SE)


# Reorder levels
mf_dir1$Model <- factor(mf_dir1$Model, levels = rev(mf_dir1$Model[1:16]))
mf_dir1$modelName <- factor(mf_dir1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_dir1$eff<-"Direct x2"
mf_all<-rbind(mf_all, mf_dir1)



### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp2 <- ggplot(mf_dir1, aes(colour = modelName, shape = modelName))
zp2 <- zp2 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp2 <- zp2 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp2 <- zp2 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp2 <- zp2 + coord_flip() + theme_bw()
zp2 <- zp2 + ylim(xlim[1], xlim[2])
zp2 <- zp2 + theme(legend.title = element_blank())
zp2 <- zp2 + labs(y="Bias", x=element_blank()) 
zp2 <- zp2 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp2 <- zp2 + ggtitle(element_blank()) 
zp2 <- zp2 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp2)




#### Indirect ####

### Figure without omv x1 ###

tmp<-res15_bias_ind1[,grep("x1", colnames(res15_bias_ind1))]
tmp_sd<-res15_sd_ind1[,grep("x1", colnames(res15_sd_ind1))]


### Set up model Frame
mf_ind1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res15_bias_ind1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind1$Model<-rep(names,7)
mf_ind1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind1$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind1$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind1$Coefficient<-as.numeric(mf_ind1$Coefficient)
mf_ind1$SE<-as.numeric(mf_ind1$SE)


# Reorder levels
mf_ind1$Model <- factor(mf_ind1$Model, levels = rev(mf_ind1$Model[1:16]))
mf_ind1$modelName <- factor(mf_ind1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind1$eff<-"Indirect x1"
mf_all<-rbind(mf_all, mf_ind1)


### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp3 <- ggplot(mf_ind1, aes(colour = modelName, shape = modelName))
zp3 <- zp3 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp3 <- zp3 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp3 <- zp3 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp3 <- zp3 + coord_flip() + theme_bw()
zp3 <- zp3 + ylim(xlim[1], xlim[2])
zp3 <- zp3 + theme(legend.title = element_blank())
zp3 <- zp3 + labs(y="Bias", x=element_blank()) 
zp3 <- zp3 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp3 <- zp3 + ggtitle(element_blank()) 
zp3 <- zp3 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp3)





### Figure without omv x2 ###

tmp<-res15_bias_ind1[,grep("x2", colnames(res15_bias_ind1))]
tmp_sd<-res15_sd_ind1[,grep("x2", colnames(res15_sd_ind1))]


### Set up model Frame
mf_ind1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res15_bias_ind1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind1$Model<-rep(names,7)
mf_ind1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind1$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind1$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind1$Coefficient<-as.numeric(mf_ind1$Coefficient)
mf_ind1$SE<-as.numeric(mf_ind1$SE)


# Reorder levels
mf_ind1$Model <- factor(mf_ind1$Model, levels = rev(mf_ind1$Model[1:16]))
mf_ind1$modelName <- factor(mf_ind1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind1$eff<-"Indirect x2"
mf_all<-rbind(mf_all, mf_ind1)

### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp4 <- ggplot(mf_ind1, aes(colour = modelName, shape = modelName))
zp4 <- zp4 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp4 <- zp4 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp4 <- zp4 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp4 <- zp4 + coord_flip() + theme_bw()
zp4 <- zp4 + ylim(xlim[1], xlim[2])
zp4 <- zp4 + theme(legend.title = element_blank())
zp4 <- zp4 + labs(y="Bias", x=element_blank()) 
zp4 <- zp4 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp4 <- zp4 + ggtitle(element_blank()) 
zp4 <- zp4 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp4)






#### Combine plots ####

mf_all$eff <- factor(mf_all$eff, levels = c("Direct x1", "Direct x2", "Indirect x1", "Indirect x2"))

mf_all$y_min<-NA
mf_all$y_min[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<--0.3
mf_all$y_min[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<--1.5

mf_all$y_max<-NA
mf_all$y_max[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<-0.7
mf_all$y_max[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<-2.2


### Plot

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

# Set limits (trim data manually)
xlim=c(-1.5, 2.2)

mf_all$lb<-mf_all$Coefficient - mf_all$SE*interval2
mf_all$ub<-mf_all$Coefficient + mf_all$SE*interval2

mf_all$ub[which(mf_all$ub>mf_all$y_max)]<-mf_all$y_max[which(mf_all$ub>mf_all$y_max)]
mf_all$lb[which(mf_all$lb<mf_all$y_min)]<-mf_all$y_min[which(mf_all$lb<mf_all$y_min)]



# Plot
zp_all <- ggplot(mf_all, aes(colour = modelName, shape = modelName))
zp_all <- zp_all + facet_grid(. ~ eff, scales="free_x")
zp_all <- zp_all + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp_all <- zp_all + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.7)
zp_all <- zp_all + geom_pointrange(aes(x = Model, y = Coefficient, ymin =lb,
                                       ymax = ub, 0.2),
                                   lwd = 0.7, position = position_dodge(width = 0.6),
                                   fill ="black")
zp_all <- zp_all + coord_flip() + theme_bw()
#zp_all <- zp_all + ylim(xlim[1], xlim[2])
zp_all <- zp_all + geom_blank(aes(y = y_min)) + geom_blank(aes(y = y_max))
zp_all <- zp_all + scale_y_continuous(expand = c(0,0) )
zp_all <- zp_all + theme(legend.title = element_blank())
zp_all <- zp_all + labs(y="Bias", x=element_blank()) 
zp_all <- zp_all + theme(legend.text = element_text(size = 20),
                         legend.position="bottom",
                         legend.key=element_blank(),
                         axis.text.x = element_text(size=16),
                         axis.text.y = element_text(size=16, colour="black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20),
                         plot.title = element_text(size=20, margin = margin(t = 10, b = -19 ), hjust=-0.235),
                         strip.background =element_blank(),
                         strip.text = element_text(size=20, colour="black"),
                         panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp_all <- zp_all + ggtitle(bquote(rho ~~~~~~ bold(delta) ~~~~~~~~~~~ lambda ~~~~~ bold(theta))) 
zp_all <- zp_all + guides(fill=guide_legend(ncol = 1), colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                          shape = guide_legend(reverse = T))
print(zp_all)




#### Combine plots ####


cairo_pdf(file=paste(od, "sim15_combined_1.pdf", sep=""), width=17.54, height=12.40, 
          bg = "white", family="CM Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
zp_all
dev.off()




########################
### Figure with omv  ###
########################

tmp<-res15_bias_dir2[,grep("x1", colnames(res15_bias_dir2))]
tmp_sd<-res15_sd_dir2[,grep("x1", colnames(res15_sd_dir2))]


### Set up model Frame
mf_dir2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res15_bias_dir2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir2$Model<-rep(names,7)
mf_dir2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir2$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir2$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir2$Coefficient<-as.numeric(mf_dir2$Coefficient)
mf_dir2$SE<-as.numeric(mf_dir2$SE)


# Reorder levels
mf_dir2$Model <- factor(mf_dir2$Model, levels = rev(mf_dir2$Model[1:16]))
mf_dir2$modelName <- factor(mf_dir2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_all<-mf_dir2
mf_all$eff<-"Direct x1"


### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp1 <- ggplot(mf_dir2, aes(colour = modelName, shape = modelName))
zp1 <- zp1 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp1 <- zp1 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp1 <- zp1 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + ylim(xlim[1], xlim[2])
zp1 <- zp1 + theme(legend.title = element_blank())
zp1 <- zp1 + labs(y="Bias", x=element_blank()) 
zp1 <- zp1 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp1 <- zp1 + ggtitle(element_blank()) 
zp1 <- zp1 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp1)




### Figure without omv x2 ###

tmp<-res15_bias_dir2[,grep("x2", colnames(res15_bias_dir2))]
tmp_sd<-res15_sd_dir2[,grep("x2", colnames(res15_sd_dir2))]


### Set up model Frame
mf_dir2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res15_bias_dir2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir2$Model<-rep(names,7)
mf_dir2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir2$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir2$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir2$Coefficient<-as.numeric(mf_dir2$Coefficient)
mf_dir2$SE<-as.numeric(mf_dir2$SE)


# Reorder levels
mf_dir2$Model <- factor(mf_dir2$Model, levels = rev(mf_dir2$Model[1:16]))
mf_dir2$modelName <- factor(mf_dir2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_dir2$eff<-"Direct x2"
mf_all<-rbind(mf_all, mf_dir2)



### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp2 <- ggplot(mf_dir2, aes(colour = modelName, shape = modelName))
zp2 <- zp2 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp2 <- zp2 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp2 <- zp2 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp2 <- zp2 + coord_flip() + theme_bw()
zp2 <- zp2 + ylim(xlim[1], xlim[2])
zp2 <- zp2 + theme(legend.title = element_blank())
zp2 <- zp2 + labs(y="Bias", x=element_blank()) 
zp2 <- zp2 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp2 <- zp2 + ggtitle(element_blank()) 
zp2 <- zp2 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp2)




#### Indirect ####

### Figure without omv x1 ###

tmp<-res15_bias_ind2[,grep("x1", colnames(res15_bias_ind2))]
tmp_sd<-res15_sd_ind2[,grep("x1", colnames(res15_sd_ind2))]


### Set up model Frame
mf_ind2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res15_bias_ind2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind2$Model<-rep(names,7)
mf_ind2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind2$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind2$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind2$Coefficient<-as.numeric(mf_ind2$Coefficient)
mf_ind2$SE<-as.numeric(mf_ind2$SE)


# Reorder levels
mf_ind2$Model <- factor(mf_ind2$Model, levels = rev(mf_ind2$Model[1:16]))
mf_ind2$modelName <- factor(mf_ind2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind2$eff<-"Indirect x1"
mf_all<-rbind(mf_all, mf_ind2)


### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp3 <- ggplot(mf_ind2, aes(colour = modelName, shape = modelName))
zp3 <- zp3 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp3 <- zp3 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp3 <- zp3 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp3 <- zp3 + coord_flip() + theme_bw()
zp3 <- zp3 + ylim(xlim[1], xlim[2])
zp3 <- zp3 + theme(legend.title = element_blank())
zp3 <- zp3 + labs(y="Bias", x=element_blank()) 
zp3 <- zp3 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp3 <- zp3 + ggtitle(element_blank()) 
zp3 <- zp3 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp3)





### Figure without omv x2 ###

tmp<-res15_bias_ind2[,grep("x2", colnames(res15_bias_ind2))]
tmp_sd<-res15_sd_ind2[,grep("x2", colnames(res15_sd_ind2))]


### Set up model Frame
mf_ind2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res15_bias_ind2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind2$Model<-rep(names,7)
mf_ind2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind2$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind2$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind2$Coefficient<-as.numeric(mf_ind2$Coefficient)
mf_ind2$SE<-as.numeric(mf_ind2$SE)


# Reorder levels
mf_ind2$Model <- factor(mf_ind2$Model, levels = rev(mf_ind2$Model[1:16]))
mf_ind2$modelName <- factor(mf_ind2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind2$eff<-"Indirect x2"
mf_all<-rbind(mf_all, mf_ind2)

### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp4 <- ggplot(mf_ind2, aes(colour = modelName, shape = modelName))
zp4 <- zp4 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp4 <- zp4 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp4 <- zp4 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp4 <- zp4 + coord_flip() + theme_bw()
zp4 <- zp4 + ylim(xlim[1], xlim[2])
zp4 <- zp4 + theme(legend.title = element_blank())
zp4 <- zp4 + labs(y="Bias", x=element_blank()) 
zp4 <- zp4 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp4 <- zp4 + ggtitle(element_blank()) 
zp4 <- zp4 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp4)






#### Combine plots ####

mf_all$eff <- factor(mf_all$eff, levels = c("Direct x1", "Direct x2", "Indirect x1", "Indirect x2"))

mf_all$y_min<-NA
mf_all$y_min[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<--0.3
mf_all$y_min[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<--1.5

mf_all$y_max<-NA
mf_all$y_max[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<-0.7
mf_all$y_max[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<-2.2


### Plot

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

# Set limits (trim data manually)
xlim=c(-1.5, 2.2)

mf_all$lb<-mf_all$Coefficient - mf_all$SE*interval2
mf_all$ub<-mf_all$Coefficient + mf_all$SE*interval2

mf_all$ub[which(mf_all$ub>mf_all$y_max)]<-mf_all$y_max[which(mf_all$ub>mf_all$y_max)]
mf_all$lb[which(mf_all$lb<mf_all$y_min)]<-mf_all$y_min[which(mf_all$lb<mf_all$y_min)]



# Plot
zp_all <- ggplot(mf_all, aes(colour = modelName, shape = modelName))
zp_all <- zp_all + facet_grid(. ~ eff, scales="free_x")
zp_all <- zp_all + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp_all <- zp_all + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.7)
zp_all <- zp_all + geom_pointrange(aes(x = Model, y = Coefficient, ymin =lb,
                                       ymax = ub, 0.2),
                                   lwd = 0.7, position = position_dodge(width = 0.6),
                                   fill ="black")
zp_all <- zp_all + coord_flip() + theme_bw()
#zp_all <- zp_all + ylim(xlim[1], xlim[2])
zp_all <- zp_all + geom_blank(aes(y = y_min)) + geom_blank(aes(y = y_max))
zp_all <- zp_all + scale_y_continuous(expand = c(0,0) )
zp_all <- zp_all + theme(legend.title = element_blank())
zp_all <- zp_all + labs(y="Bias", x=element_blank()) 
zp_all <- zp_all + theme(legend.text = element_text(size = 20),
                         legend.position="bottom",
                         legend.key=element_blank(),
                         axis.text.x = element_text(size=16),
                         axis.text.y = element_text(size=16, colour="black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20),
                         plot.title = element_text(size=20, margin = margin(t = 10, b = -19 ), hjust=-0.235),
                         strip.background =element_blank(),
                         strip.text = element_text(size=20, colour="black"),
                         panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp_all <- zp_all + ggtitle(bquote(rho ~~~~~~ bold(delta) ~~~~~~~~~~~ lambda ~~~~~ bold(theta))) 
zp_all <- zp_all + guides(fill=guide_legend(ncol = 1), colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                          shape = guide_legend(reverse = T))
print(zp_all)




#### Combine plots ####


cairo_pdf(file=paste(od, "sim15_combined_2.pdf", sep=""), width=17.54, height=12.40, 
          bg = "white", family="CM Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
zp_all
dev.off()





######################################
######################################
#### 4) K NB Weights              ####
######################################
######################################

source(paste(dd, "01e_Monte Carlo Simulation Spatial_Program_knbW.R", sep="/"))


######################################
#### Define possible combinations ####
######################################

beta<-c(0.2, 0.5)

rho<-c(0, 0.5)
delta<-c("0, 0", "0.4, 0.7")
lambda<-c(0, 0.5)
gamma<-c(0, 0.3)
theta<-c("0, 0", "0.1, 0.8")


### Set up matrix of all combinations ###

nc<-2^5

l<-list(rho, delta, lambda, gamma, theta)

comb_mat4<-expand.grid(l)

colnames(comb_mat4)<-c("rho", "delta", "lambda", "gamma", "theta")
comb_mat4$delta<-as.character(comb_mat4$delta)
comb_mat4$theta<-as.character(comb_mat4$theta)



###########################################
#### Start simulations of combinations ####
###########################################

# Set up manual loops (to spare capacity)


k<-c(1:32)

### Simumlation Run
for(i in k){

  sim16<-NA # avoid using the last sim (in case of error)

  par<-comb_mat4[i,]
  deltai<-as.numeric(c(unlist(strsplit(par$delta, split=", "))))
  thetai<-as.numeric(c(unlist(strsplit(par$theta, split=", "))))
  sim16<-spsim(N=900, R=1000, seed=1573262181,
              ay=par$rho, ax1=deltai[1], ax2=deltai[2], au=par$lambda,
              bx1=beta[1], bwx1=thetai[1], bx2=beta[2], bwx2=thetai[2],
              bu=1, bwu=0, sdu=1,
              omv1=par$gamma, womv1=0)
  assign(paste0("sim16_", i, sep=""), sim16)



  save(list=c(paste0("sim16_", i, sep="")), file=paste0("Simulation16_", i, ".RData", sep=""))

}


for(i in k){
  load(paste("Simulation16_", i, ".RData", sep=""))
}





#######################################
#### Make plots of each Simulation ####
#######################################


### Jpeg ###

for(i in k){
  cm<-comb_mat4[i,]
  
  tmp.sim<-get(paste("sim16_", i, sep=""))
  
  p<-plot.sim(tmp.sim, xlim=c(-2, 3))
  
  name = paste(od, "sim16_", i, ".jpeg", sep="")
  jpeg(file = name, bg = "transparent", width = 2720, height = 1820, quality=100, family="CM Roman", res=300)
  par(mar=c(0,0,0,0))
  par(mfrow=c(1,1),oma=c(0,0,0,0))
  
  print(p + ggtitle(bquote(rho == .(cm$rho) ~~
                             delta == .(cm$delta) ~~
                             lambda == .(cm$lambda) ~~
                             gamma == .(cm$gamma) ~~
                             theta == .(cm$theta))))
  dev.off()
  par(mfrow=c(1,1))
}





#####################################
#### Create result matrices RMSE ####
#####################################


### Direct impacts ###

res16_rmse_dir<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res16_rmse_dir)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res16_rmse_dir<-cbind(res16_rmse_dir,tmp)

for(i in 1:nrow(res16_rmse_dir)){
  tmp.sim<-get(paste("sim16_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_d[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res16_rmse_dir[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res16_rmse_dir[i, 6:ncol(res16_rmse_dir)]<-t(m[,2])
}



### Separate into omv and non-omv

res16_rmse_dir1<-res16_rmse_dir[res16_rmse_dir$gamma==0,]
res16_rmse_dir2<-res16_rmse_dir[res16_rmse_dir$gamma>0,]

res16_rmse_dir1<-res16_rmse_dir1[,-which(names(res16_rmse_dir1)=="gamma")]
res16_rmse_dir2<-res16_rmse_dir2[,-which(names(res16_rmse_dir2)=="gamma")]


res16_rmse_dir1[,1:4]<-apply(res16_rmse_dir1[,1:4], 2, FUN=function(x) as.character(x))
res16_rmse_dir2[,1:4]<-apply(res16_rmse_dir2[,1:4], 2, FUN=function(x) as.character(x))





### Indirect impacts ###

res16_rmse_ind<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res16_rmse_ind)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res16_rmse_ind<-cbind(res16_rmse_ind,tmp)

for(i in 1:nrow(res16_rmse_ind)){
  tmp.sim<-get(paste("sim16_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sp[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res16_rmse_ind[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res16_rmse_ind[i, 6:ncol(res16_rmse_ind)]<-t(m[,2])
}


res16_rmse_ind<-res16_rmse_ind[,-which(names(res16_rmse_ind) %in% c("OLS_x1", "OLS_x2", "SEM_x1", "SEM_x2"))]

### Separate into omv and non-omv

res16_rmse_ind1<-res16_rmse_ind[res16_rmse_ind$gamma==0,]
res16_rmse_ind2<-res16_rmse_ind[res16_rmse_ind$gamma>0,]

res16_rmse_ind1<-res16_rmse_ind1[,-which(names(res16_rmse_ind1)=="gamma")]
res16_rmse_ind2<-res16_rmse_ind2[,-which(names(res16_rmse_ind2)=="gamma")]


res16_rmse_ind1[,1:4]<-apply(res16_rmse_ind1[,1:4], 2, FUN=function(x) as.character(x))
res16_rmse_ind2[,1:4]<-apply(res16_rmse_ind2[,1:4], 2, FUN=function(x) as.character(x))



#######################
### Wirte tex table ###
#######################


matout<-function(x, file="", digits=NULL, colav=T, threepart=T, stars=F,
                 code_before=NULL, code_after=NULL,
                 header1=NULL, header2=NULL, 
                 caption=NULL, label=NULL, note=NULL,...) {
  
  ### Options for xtable
  default_args = list(include.colnames=FALSE, only.contents=TRUE,
                      include.rownames=FALSE, hline.after=NULL, comment=FALSE,
                      print.results=FALSE, sanitize.text.function = function(z) z)
  passed_args = list(...)
  
  
  if(is.null(caption)){caption<-"Table"}
  
  c<-ncol(x)
  
  if(is.null(code_before)){
    code_before<-paste("\\begin{tabular}{", paste(rep("l", c), collapse = " "), "}")
  }
  if(is.null(code_after)){
    code_after<-paste("\\end{tabular}")
  }
  
  num<-which(sapply(x,FUN=function(x) is.numeric(x)))
  
  if(colav==T){
    av<-data.frame(t(colMeans(abs(x[,num[1]:ncol(x)]))))
    av2<-data.frame(t(colMeans(matrix(t(colMeans(abs(x[,num[1]:ncol(x)]))), 2))))
    
    if(stars){
      pa1<-colnames(av[seq(1, ncol(av), 2)])[which.min(abs(av[seq(1, ncol(av), 2)]))]
      pa2<-colnames(av[seq(2, ncol(av), 2)])[which.min(abs(av[seq(2, ncol(av), 2)]))]
      pa3<-colnames(av2)[which.min(av2)]
      
      #av[1,]<-round(av[1,], digits)
      av[1,]<-sprintf(paste("%.", digits, "f"), round(av[1,],digits))
      av2[1,]<-sprintf(paste("%.", digits, "f"), round(av2[1,],digits))
      
      
      av[1,colnames(av)==pa1]<-paste(av[1,colnames(av)==pa1], "\\ensuremath{^{*}}", sep="")
      av[1,colnames(av)==pa2]<-paste(av[1,colnames(av)==pa2], "\\ensuremath{^{*}}", sep="")
      av2[1,colnames(av2)==pa3]<-paste(av2[1,colnames(av2)==pa3], "\\ensuremath{^{*}}", sep="")
      
      
    }
    
    calling_args2 = c(list(x=xtable(av, digits=digits)),
                      c(passed_args,
                        default_args[setdiff(names(default_args), names(passed_args))]))
    
  }
  
  if(stars){
    p11<-apply(x[,seq(num[1], ncol(x), 2)], 1, FUN=function(z) which.min(abs(z)))
    p11<-colnames(x[,seq(num[1], ncol(x), 2)])[p11]
    
    p12<-apply(x[,seq(num[1]+1, ncol(x), 2)], 1, FUN=function(z) which.min(abs(z)))
    p12<-colnames(x[,seq(num[1]+1, ncol(x), 2)])[p12]
    
    #x[, num]<-round(x[, num], digits)
    #x[, num]<-apply(x[, num], 2, FUN=function(z) as.character(format(z, digits=digits, scientific=F)))
    x[, num]<-apply(x[, num], 2, FUN=function(z) sprintf(paste("%.", digits, "f"), round(z,digits)))
    
    
    for(i in 1:nrow(x)){
      x[i,colnames(x)==p11[i]]<-paste(x[i,colnames(x)==p11[i]], "\\ensuremath{^{*}}", sep="")
      x[i,colnames(x)==p12[i]]<-paste(x[i,colnames(x)==p12[i]], "\\ensuremath{^{*}}", sep="")
    }
    
    
  }
  
  calling_args = c(list(x=xtable(x, digits=digits)),
                   c(passed_args,
                     default_args[setdiff(names(default_args), names(passed_args))]))
  
  midrule<-paste(seq(num[1], max(num), 2), "-", seq(num[1]+1, max(num), 2), sep="")
  midrule<-c(paste("\\cmidrule(r){", midrule[1], "}", sep=""), 
             paste("\\cmidrule(lr){", midrule[-c(1, length(midrule))], "}", sep=""),
             paste("\\cmidrule(l){", midrule[length(midrule)], "}", sep=""))
  
  
  ### Output
  cat("\\begin{table}\n", file=file)
  cat("\\centering\n", file=file, append=T)
  
  if(threepart==T){
    cat("{\\begin{threeparttable}\n", file=file, append=T)
  }
  
  cat(paste0("\\caption{", caption, "}\n"), file=file, append=T)
  
  if(!is.null(label)){
    cat(paste0("\\label{", label, "}\n"), file=file, append=T)
  }
  
  cat(paste0(code_before, "\n"), file=file, append=T)
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  
  if(!is.null(header1)){
    cat(paste(header1), "\\\\ \n", file=file, append=T)
  }else{cat(paste(paste("M", seq(1, c,), collapse=" & "), "\\\\ \n"), file=file, append=T)}
  
  cat(midrule, file=file, append=T)
  
  if(!is.null(header2)){
    cat(paste(header2), "\\\\ \n", file=file, append=T)
  }
  
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  cat(do.call(print.xtable, calling_args), file=file, append=T)
  cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  
  if(colav==T){
    cat(paste("\\multicolumn{", (num[1]-1), "}{l}{Average} &", collapse=""), file=file, append=T)
    cat(do.call(print.xtable, calling_args2), file=file, append=T)
    cat(paste("\\multicolumn{", (num[1]-1), "}{l}{} &", collapse=""), file=file, append=T)
    cat(paste("\\multicolumn{2}{c}{", av2, "}", collapse=" & "), "\\\\ \n", file=file, append=T)
    cat("\\hline \\\\[-1.8ex]\n", file=file, append=T)
  }
  
  cat(paste0(code_after, "\n"), file=file, append=T)
  
  if(!is.null(note)){
    cat("\\begin{tablenotes}\n", file=file, append=T)
    cat(paste("\\item \\scriptsize{", note, "}\n", collape=""), file=file, append=T)
    cat("\\end{tablenotes}\n", file=file, append=T)
  }
  
  if(threepart==T){
    cat("\\end{threeparttable}}\n", file=file, append=T)
  }
  
  cat("\\end{table}\n", file=file, append=T)
}



### Table 1

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 7), collapse=" & "), sep=" & " )

matout(res16_rmse_dir1, file=paste(od, "res16_rmse_dir1.tex", sep=""), digits=4, stars=T,
       caption="RMSE of direct impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res16_dir1",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 2

matout(res16_rmse_dir2, file=paste(od, "res16_rmse_dir2.tex", sep=""), digits=4, stars=T,
       caption="RMSE of direct impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res16_dir2",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")



### Table 3

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("SLX", "SAR", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 5), collapse=" & "), sep=" & " )

matout(res16_rmse_ind1, file=paste(od, "res16_rmse_ind1.tex", sep=""), digits=4, stars=T,
       caption="RMSE of indirect impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res16_ind1",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 4

matout(res16_rmse_ind2, file=paste(od, "res16_rmse_ind2.tex", sep=""), digits=4, stars=T,
       caption="RMSE of indirect impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res16_ind2",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest RMSE for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")




#####################################
#### Create result matrices Bias ####
#####################################


### Direct impacts ###

res16_bias_dir<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res16_bias_dir)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res16_bias_dir<-cbind(res16_bias_dir,tmp)

for(i in 1:nrow(res16_bias_dir)){
  tmp.sim<-get(paste("sim16_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_d[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res16_bias_dir[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res16_bias_dir[i, 6:ncol(res16_bias_dir)]<-t(m[,1])
}



### Separate into omv and non-omv

res16_bias_dir1<-res16_bias_dir[res16_bias_dir$gamma==0,]
res16_bias_dir2<-res16_bias_dir[res16_bias_dir$gamma>0,]

res16_bias_dir1<-res16_bias_dir1[,-which(names(res16_bias_dir1)=="gamma")]
res16_bias_dir2<-res16_bias_dir2[,-which(names(res16_bias_dir2)=="gamma")]


res16_bias_dir1[,1:4]<-apply(res16_bias_dir1[,1:4], 2, FUN=function(x) as.character(x))
res16_bias_dir2[,1:4]<-apply(res16_bias_dir2[,1:4], 2, FUN=function(x) as.character(x))





### Indirect impacts ###

res16_bias_ind<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res16_bias_ind)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res16_bias_ind<-cbind(res16_bias_ind,tmp)

for(i in 1:nrow(res16_bias_ind)){
  tmp.sim<-get(paste("sim16_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sp[-seq(1, nrow(s$mat_d), by=3),]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res16_bias_ind[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res16_bias_ind[i, 6:ncol(res16_bias_ind)]<-t(m[,1])
}


res16_bias_ind<-res16_bias_ind[,-which(names(res16_bias_ind) %in% c("OLS_x1", "OLS_x2", "SEM_x1", "SEM_x2"))]

### Separate into omv and non-omv

res16_bias_ind1<-res16_bias_ind[res16_bias_ind$gamma==0,]
res16_bias_ind2<-res16_bias_ind[res16_bias_ind$gamma>0,]

res16_bias_ind1<-res16_bias_ind1[,-which(names(res16_bias_ind1)=="gamma")]
res16_bias_ind2<-res16_bias_ind2[,-which(names(res16_bias_ind2)=="gamma")]


res16_bias_ind1[,1:4]<-apply(res16_bias_ind1[,1:4], 2, FUN=function(x) as.character(x))
res16_bias_ind2[,1:4]<-apply(res16_bias_ind2[,1:4], 2, FUN=function(x) as.character(x))



#######################
### Wirte tex table ###
#######################



### Table 1

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 7), collapse=" & "), sep=" & " )

matout(res16_bias_dir1, file=paste(od, "res16_bias_dir1.tex", sep=""), digits=4, stars=T,
       caption="Bias of direct impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res16_dir1b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 2

matout(res16_bias_dir2, file=paste(od, "res16_bias_dir2.tex", sep=""), digits=4, stars=T,
       caption="Bias of direct impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res16_dir2b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")



### Table 3

header1<-paste(" & & & & ",
               paste("\\multicolumn{2}{c}{", c("SLX", "SAR", "SAC", "SDM", "SDEM"), "}", collapse=" & "), sep="")

header2<-paste(paste("$\\rho$", "${\\bm\\delta}$", "$\\lambda$", "${\\bm\\theta}$", sep=" & "), paste(rep("\\tc{$\\bm x_1$}  & \\tc{$\\bm x_2$}", 5), collapse=" & "), sep=" & " )

matout(res16_bias_ind1, file=paste(od, "res16_bias_ind1.tex", sep=""), digits=4, stars=T,
       caption="Bias of indirect impacts without omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0, 0)^\\intercal$.", 
       label="table:res16_ind1b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")


### Table 4

matout(res16_bias_ind2, file=paste(od, "res16_bias_ind2.tex", sep=""), digits=4, stars=T,
       caption="Bias of indirect impacts with omv. ${\\bm \\beta}=(0.2, 0.5)^\\intercal$, ${\\bm \\gamma}=(0.3, 0)^\\intercal$.", 
       label="table:res16_ind2b",
       code_before=c(paste("\\begin{tabular}{llll", paste(rep("D{.}{.}{1.5}", 14), collapse=" "), "}", collapse=" ")),
       header1=header1, header2=header2,
       note="\\ensuremath{^{*}} Lowest bias for ${\\bm x}_k$ within the parameter combination. Number of observations=900, repetitions=1000.")






##############################
#### Create combined plot ####
##############################


### SD Direct impacts ###

res16_sd_dir<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res16_sd_dir)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res16_sd_dir<-cbind(res16_sd_dir,tmp)

for(i in 1:nrow(res16_sd_dir)){
  tmp.sim<-get(paste("sim16_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sd[-seq(1, nrow(s$mat_d), by=3),1]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res16_sd_dir[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res16_sd_dir[i, 6:ncol(res16_sd_dir)]<-t(m)
}



### Separate into omv and non-omv

res16_sd_dir1<-res16_sd_dir[res16_sd_dir$gamma==0,]
res16_sd_dir2<-res16_sd_dir[res16_sd_dir$gamma>0,]

res16_sd_dir1<-res16_sd_dir1[,-which(names(res16_sd_dir1)=="gamma")]
res16_sd_dir2<-res16_sd_dir2[,-which(names(res16_sd_dir2)=="gamma")]


res16_sd_dir1[,1:4]<-apply(res16_sd_dir1[,1:4], 2, FUN=function(x) as.character(x))
res16_sd_dir2[,1:4]<-apply(res16_sd_dir2[,1:4], 2, FUN=function(x) as.character(x))





### SD Indirect impacts ###

res16_sd_ind<-comb_mat4

tmp<-data.frame(matrix(NA, ncol=14, nrow=nrow(res16_sd_ind)))
names(tmp)<-c("OLS_x1", "OLS_x2", "SLX_x1", "SLX_x2", "SAR_x1", "SAR_x2", "SEM_x1", "SEM_x2", "SAC_x1",
              "SAC_x2", "SDM_x1", "SDM_x2", "SDEM_x1", "SDEM_x2")

res16_sd_ind<-cbind(res16_sd_ind,tmp)

for(i in 1:nrow(res16_sd_ind)){
  tmp.sim<-get(paste("sim16_", i, sep=""))
  par<-tmp.sim$parameters
  s<-summary.sim(tmp.sim)
  m<-s$mat_sd[-seq(1, nrow(s$mat_sd), by=3),2]
  
  # Test for correct row
  par2<-as.character(c(par[1],paste(par[2],par[3], sep=", "), par[4], par[11], paste(par[6],par[8], sep=", ")))
  if(!identical(as.character(res16_sd_ind[i, 1:5]), par2)){
    cat("\n", "Parameters not identical at ", i, "\n")
    break
  }
  
  res16_sd_ind[i, 6:ncol(res16_sd_ind)]<-t(m)
}


res16_sd_ind<-res16_sd_ind[,-which(names(res16_sd_ind) %in% c("OLS_x1", "OLS_x2", "SEM_x1", "SEM_x2"))]

### Separate into omv and non-omv

res16_sd_ind1<-res16_sd_ind[res16_sd_ind$gamma==0,]
res16_sd_ind2<-res16_sd_ind[res16_sd_ind$gamma>0,]

res16_sd_ind1<-res16_sd_ind1[,-which(names(res16_sd_ind1)=="gamma")]
res16_sd_ind2<-res16_sd_ind2[,-which(names(res16_sd_ind2)=="gamma")]


res16_sd_ind1[,1:4]<-apply(res16_sd_ind1[,1:4], 2, FUN=function(x) as.character(x))
res16_sd_ind2[,1:4]<-apply(res16_sd_ind2[,1:4], 2, FUN=function(x) as.character(x))


### Expand strings
res16_bias_dir1$delta[which(res16_bias_dir1$delta=="0, 0")]<-c("0.0, 0.0")
res16_bias_ind1$delta[which(res16_bias_ind1$delta=="0, 0")]<-c("0.0, 0.0")
res16_bias_dir1$theta[which(res16_bias_dir1$theta=="0, 0")]<-c("0.0, 0.0")
res16_bias_ind1$theta[which(res16_bias_ind1$theta=="0, 0")]<-c("0.0, 0.0")

res16_bias_dir2$delta[which(res16_bias_dir2$delta=="0, 0")]<-c("0.0, 0.0")
res16_bias_ind2$delta[which(res16_bias_ind2$delta=="0, 0")]<-c("0.0, 0.0")
res16_bias_dir2$theta[which(res16_bias_dir2$theta=="0, 0")]<-c("0.0, 0.0")
res16_bias_ind2$theta[which(res16_bias_ind2$theta=="0, 0")]<-c("0.0, 0.0")


# res16_bias_dir1$delta<-sprintf("%-15s", res16_bias_dir1$delta)
# res16_bias_dir1$theta<-sprintf("%-15s", res16_bias_dir1$theta)
# res16_bias_ind1$delta<-sprintf("%-15s", res16_bias_ind1$delta)
# res16_bias_ind1$theta<-sprintf("%-15s", res16_bias_ind1$theta)




###########################
### Figure without omv  ###
###########################

tmp<-res16_bias_dir1[,grep("x1", colnames(res16_bias_dir1))]
tmp_sd<-res16_sd_dir1[,grep("x1", colnames(res16_sd_dir1))]


### Set up model Frame
mf_dir1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res16_bias_dir1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir1$Model<-rep(names,7)
mf_dir1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir1$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir1$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir1$Coefficient<-as.numeric(mf_dir1$Coefficient)
mf_dir1$SE<-as.numeric(mf_dir1$SE)


# Reorder levels
mf_dir1$Model <- factor(mf_dir1$Model, levels = rev(mf_dir1$Model[1:16]))
mf_dir1$modelName <- factor(mf_dir1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_all<-mf_dir1
mf_all$eff<-"Direct x1"


### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp1 <- ggplot(mf_dir1, aes(colour = modelName, shape = modelName))
zp1 <- zp1 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp1 <- zp1 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp1 <- zp1 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + ylim(xlim[1], xlim[2])
zp1 <- zp1 + theme(legend.title = element_blank())
zp1 <- zp1 + labs(y="Bias", x=element_blank()) 
zp1 <- zp1 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp1 <- zp1 + ggtitle(element_blank()) 
zp1 <- zp1 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp1)




### Figure without omv x2 ###

tmp<-res16_bias_dir1[,grep("x2", colnames(res16_bias_dir1))]
tmp_sd<-res16_sd_dir1[,grep("x2", colnames(res16_sd_dir1))]


### Set up model Frame
mf_dir1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res16_bias_dir1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir1$Model<-rep(names,7)
mf_dir1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir1$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir1$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir1$Coefficient<-as.numeric(mf_dir1$Coefficient)
mf_dir1$SE<-as.numeric(mf_dir1$SE)


# Reorder levels
mf_dir1$Model <- factor(mf_dir1$Model, levels = rev(mf_dir1$Model[1:16]))
mf_dir1$modelName <- factor(mf_dir1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_dir1$eff<-"Direct x2"
mf_all<-rbind(mf_all, mf_dir1)



### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp2 <- ggplot(mf_dir1, aes(colour = modelName, shape = modelName))
zp2 <- zp2 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp2 <- zp2 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp2 <- zp2 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp2 <- zp2 + coord_flip() + theme_bw()
zp2 <- zp2 + ylim(xlim[1], xlim[2])
zp2 <- zp2 + theme(legend.title = element_blank())
zp2 <- zp2 + labs(y="Bias", x=element_blank()) 
zp2 <- zp2 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp2 <- zp2 + ggtitle(element_blank()) 
zp2 <- zp2 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp2)




#### Indirect ####

### Figure without omv x1 ###

tmp<-res16_bias_ind1[,grep("x1", colnames(res16_bias_ind1))]
tmp_sd<-res16_sd_ind1[,grep("x1", colnames(res16_sd_ind1))]


### Set up model Frame
mf_ind1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res16_bias_ind1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind1$Model<-rep(names,7)
mf_ind1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind1$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind1$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind1$Coefficient<-as.numeric(mf_ind1$Coefficient)
mf_ind1$SE<-as.numeric(mf_ind1$SE)


# Reorder levels
mf_ind1$Model <- factor(mf_ind1$Model, levels = rev(mf_ind1$Model[1:16]))
mf_ind1$modelName <- factor(mf_ind1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind1$eff<-"Indirect x1"
mf_all<-rbind(mf_all, mf_ind1)


### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp3 <- ggplot(mf_ind1, aes(colour = modelName, shape = modelName))
zp3 <- zp3 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp3 <- zp3 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp3 <- zp3 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp3 <- zp3 + coord_flip() + theme_bw()
zp3 <- zp3 + ylim(xlim[1], xlim[2])
zp3 <- zp3 + theme(legend.title = element_blank())
zp3 <- zp3 + labs(y="Bias", x=element_blank()) 
zp3 <- zp3 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp3 <- zp3 + ggtitle(element_blank()) 
zp3 <- zp3 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp3)





### Figure without omv x2 ###

tmp<-res16_bias_ind1[,grep("x2", colnames(res16_bias_ind1))]
tmp_sd<-res16_sd_ind1[,grep("x2", colnames(res16_sd_ind1))]


### Set up model Frame
mf_ind1<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind1)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res16_bias_ind1[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind1$Model<-rep(names,7)
mf_ind1$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind1$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind1$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind1$Coefficient<-as.numeric(mf_ind1$Coefficient)
mf_ind1$SE<-as.numeric(mf_ind1$SE)


# Reorder levels
mf_ind1$Model <- factor(mf_ind1$Model, levels = rev(mf_ind1$Model[1:16]))
mf_ind1$modelName <- factor(mf_ind1$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind1$eff<-"Indirect x2"
mf_all<-rbind(mf_all, mf_ind1)

### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp4 <- ggplot(mf_ind1, aes(colour = modelName, shape = modelName))
zp4 <- zp4 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp4 <- zp4 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp4 <- zp4 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp4 <- zp4 + coord_flip() + theme_bw()
zp4 <- zp4 + ylim(xlim[1], xlim[2])
zp4 <- zp4 + theme(legend.title = element_blank())
zp4 <- zp4 + labs(y="Bias", x=element_blank()) 
zp4 <- zp4 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp4 <- zp4 + ggtitle(element_blank()) 
zp4 <- zp4 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp4)






#### Combine plots ####

mf_all$eff <- factor(mf_all$eff, levels = c("Direct x1", "Direct x2", "Indirect x1", "Indirect x2"))

mf_all$y_min<-NA
mf_all$y_min[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<--0.3
mf_all$y_min[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<--1.5

mf_all$y_max<-NA
mf_all$y_max[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<-0.7
mf_all$y_max[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<-2.2


### Plot

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

# Set limits (trim data manually)
xlim=c(-1.5, 2.2)

mf_all$lb<-mf_all$Coefficient - mf_all$SE*interval2
mf_all$ub<-mf_all$Coefficient + mf_all$SE*interval2

mf_all$ub[which(mf_all$ub>mf_all$y_max)]<-mf_all$y_max[which(mf_all$ub>mf_all$y_max)]
mf_all$lb[which(mf_all$lb<mf_all$y_min)]<-mf_all$y_min[which(mf_all$lb<mf_all$y_min)]



# Plot
zp_all <- ggplot(mf_all, aes(colour = modelName, shape = modelName))
zp_all <- zp_all + facet_grid(. ~ eff, scales="free_x")
zp_all <- zp_all + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp_all <- zp_all + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.7)
zp_all <- zp_all + geom_pointrange(aes(x = Model, y = Coefficient, ymin =lb,
                                       ymax = ub, 0.2),
                                   lwd = 0.7, position = position_dodge(width = 0.6),
                                   fill ="black")
zp_all <- zp_all + coord_flip() + theme_bw()
#zp_all <- zp_all + ylim(xlim[1], xlim[2])
zp_all <- zp_all + geom_blank(aes(y = y_min)) + geom_blank(aes(y = y_max))
zp_all <- zp_all + scale_y_continuous(expand = c(0,0) )
zp_all <- zp_all + theme(legend.title = element_blank())
zp_all <- zp_all + labs(y="Bias", x=element_blank()) 
zp_all <- zp_all + theme(legend.text = element_text(size = 20),
                         legend.position="bottom",
                         legend.key=element_blank(),
                         axis.text.x = element_text(size=16),
                         axis.text.y = element_text(size=16, colour="black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20),
                         plot.title = element_text(size=20, margin = margin(t = 10, b = -19 ), hjust=-0.235),
                         strip.background =element_blank(),
                         strip.text = element_text(size=20, colour="black"),
                         panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp_all <- zp_all + ggtitle(bquote(rho ~~~~~~ bold(delta) ~~~~~~~~~~~ lambda ~~~~~ bold(theta))) 
zp_all <- zp_all + guides(fill=guide_legend(ncol = 1), colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                          shape = guide_legend(reverse = T))
print(zp_all)




#### Combine plots ####


cairo_pdf(file=paste(od, "sim16_combined_1.pdf", sep=""), width=17.54, height=12.40, 
          bg = "white", family="CM Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
zp_all
dev.off()








########################
### Figure with omv  ###
########################

tmp<-res16_bias_dir2[,grep("x1", colnames(res16_bias_dir2))]
tmp_sd<-res16_sd_dir2[,grep("x1", colnames(res16_sd_dir2))]


### Set up model Frame
mf_dir2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res16_bias_dir2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir2$Model<-rep(names,7)
mf_dir2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir2$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir2$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir2$Coefficient<-as.numeric(mf_dir2$Coefficient)
mf_dir2$SE<-as.numeric(mf_dir2$SE)


# Reorder levels
mf_dir2$Model <- factor(mf_dir2$Model, levels = rev(mf_dir2$Model[1:16]))
mf_dir2$modelName <- factor(mf_dir2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_all<-mf_dir2
mf_all$eff<-"Direct x1"


### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp1 <- ggplot(mf_dir2, aes(colour = modelName, shape = modelName))
zp1 <- zp1 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp1 <- zp1 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp1 <- zp1 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + ylim(xlim[1], xlim[2])
zp1 <- zp1 + theme(legend.title = element_blank())
zp1 <- zp1 + labs(y="Bias", x=element_blank()) 
zp1 <- zp1 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp1 <- zp1 + ggtitle(element_blank()) 
zp1 <- zp1 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp1)




### Figure without omv x2 ###

tmp<-res16_bias_dir2[,grep("x2", colnames(res16_bias_dir2))]
tmp_sd<-res16_sd_dir2[,grep("x2", colnames(res16_sd_dir2))]


### Set up model Frame
mf_dir2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_dir2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res16_bias_dir2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_dir2$Model<-rep(names,7)
mf_dir2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_dir2$Coefficient<-c(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], tmp[,6], tmp[,7])
mf_dir2$SE<-c(tmp_sd[,1], tmp_sd[,2], tmp_sd[,3], tmp_sd[,4], tmp_sd[,5], tmp_sd[,6], tmp_sd[,7])

mf_dir2$Coefficient<-as.numeric(mf_dir2$Coefficient)
mf_dir2$SE<-as.numeric(mf_dir2$SE)


# Reorder levels
mf_dir2$Model <- factor(mf_dir2$Model, levels = rev(mf_dir2$Model[1:16]))
mf_dir2$modelName <- factor(mf_dir2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_dir2$eff<-"Direct x2"
mf_all<-rbind(mf_all, mf_dir2)



### Plot
xlim=c(-0.5, 0.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp2 <- ggplot(mf_dir2, aes(colour = modelName, shape = modelName))
zp2 <- zp2 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp2 <- zp2 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp2 <- zp2 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp2 <- zp2 + coord_flip() + theme_bw()
zp2 <- zp2 + ylim(xlim[1], xlim[2])
zp2 <- zp2 + theme(legend.title = element_blank())
zp2 <- zp2 + labs(y="Bias", x=element_blank()) 
zp2 <- zp2 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp2 <- zp2 + ggtitle(element_blank()) 
zp2 <- zp2 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp2)




#### Indirect ####

### Figure without omv x1 ###

tmp<-res16_bias_ind2[,grep("x1", colnames(res16_bias_ind2))]
tmp_sd<-res16_sd_ind2[,grep("x1", colnames(res16_sd_ind2))]


### Set up model Frame
mf_ind2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res16_bias_ind2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind2$Model<-rep(names,7)
mf_ind2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind2$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind2$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind2$Coefficient<-as.numeric(mf_ind2$Coefficient)
mf_ind2$SE<-as.numeric(mf_ind2$SE)


# Reorder levels
mf_ind2$Model <- factor(mf_ind2$Model, levels = rev(mf_ind2$Model[1:16]))
mf_ind2$modelName <- factor(mf_ind2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind2$eff<-"Indirect x1"
mf_all<-rbind(mf_all, mf_ind2)


### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp3 <- ggplot(mf_ind2, aes(colour = modelName, shape = modelName))
zp3 <- zp3 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp3 <- zp3 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp3 <- zp3 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp3 <- zp3 + coord_flip() + theme_bw()
zp3 <- zp3 + ylim(xlim[1], xlim[2])
zp3 <- zp3 + theme(legend.title = element_blank())
zp3 <- zp3 + labs(y="Bias", x=element_blank()) 
zp3 <- zp3 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp3 <- zp3 + ggtitle(element_blank()) 
zp3 <- zp3 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp3)





### Figure without omv x2 ###

tmp<-res16_bias_ind2[,grep("x2", colnames(res16_bias_ind2))]
tmp_sd<-res16_sd_ind2[,grep("x2", colnames(res16_sd_ind2))]


### Set up model Frame
mf_ind2<-data.frame(matrix(NA, ncol=4, nrow=16*7))
colnames(mf_ind2)<-c("Model", "Coefficient", "SE", "modelName")

names<-apply(res16_bias_ind2[,1:4], 1, function(x) paste(x, collapse="    "))
names<-paste(paste(c(1:16), ")", sep=""), names, sep="    ")
mf_ind2$Model<-rep(names,7)
mf_ind2$modelName<-rep(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM"), each=16)

mf_ind2$Coefficient<-c(rep(NA, 16), tmp[,1], tmp[,2], rep(NA, 16), tmp[,3], tmp[,4], tmp[,5])
mf_ind2$SE<-c(rep(NA, 16), tmp_sd[,1], tmp_sd[,2], rep(NA, 16), tmp_sd[,3], tmp_sd[,4], tmp_sd[,5])

mf_ind2$Coefficient<-as.numeric(mf_ind2$Coefficient)
mf_ind2$SE<-as.numeric(mf_ind2$SE)


# Reorder levels
mf_ind2$Model <- factor(mf_ind2$Model, levels = rev(mf_ind2$Model[1:16]))
mf_ind2$modelName <- factor(mf_ind2$modelName, levels = rev(c("OLS", "SLX", "SAR", "SEM", "SAC", "SDM", "SDEM")))

mf_ind2$eff<-"Indirect x2"
mf_all<-rbind(mf_all, mf_ind2)

### Plot
xlim=c(-1.5, 2.8)

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier


# Plot
zp4 <- ggplot(mf_ind2, aes(colour = modelName, shape = modelName))
zp4 <- zp4 + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp4 <- zp4 + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.8)
zp4 <- zp4 + geom_pointrange(aes(x = Model, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2, 0.2),
                             lwd = 0.8, position = position_dodge(width = 1/2),
                             fill ="black")
zp4 <- zp4 + coord_flip() + theme_bw()
zp4 <- zp4 + ylim(xlim[1], xlim[2])
zp4 <- zp4 + theme(legend.title = element_blank())
zp4 <- zp4 + labs(y="Bias", x=element_blank()) 
zp4 <- zp4 + theme(legend.text = element_text(size = 20),
                   legend.position="bottom",
                   legend.key=element_blank(),
                   axis.text.x = element_text(size=20),
                   axis.text.y = element_text(size=20, colour="black"),
                   axis.title.x = element_text(size=20),
                   axis.title.y = element_text(size=20),
                   plot.title = element_text(size=20),
                   panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp4 <- zp4 + ggtitle(element_blank()) 
zp4 <- zp4 + guides(colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                    shape = guide_legend(reverse = T))
print(zp4)






#### Combine plots ####

mf_all$eff <- factor(mf_all$eff, levels = c("Direct x1", "Direct x2", "Indirect x1", "Indirect x2"))

mf_all$y_min<-NA
mf_all$y_min[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<--0.3
mf_all$y_min[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<--1.5

mf_all$y_max<-NA
mf_all$y_max[which(mf_all$eff %in% c("Direct x1", "Direct x2"))]<-0.7
mf_all$y_max[which(mf_all$eff %in% c("Indirect x1", "Indirect x2"))]<-2.2


### Plot

# Confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

# Set limits (trim data manually)
xlim=c(-1.5, 2.2)

mf_all$lb<-mf_all$Coefficient - mf_all$SE*interval2
mf_all$ub<-mf_all$Coefficient + mf_all$SE*interval2

mf_all$ub[which(mf_all$ub>mf_all$y_max)]<-mf_all$y_max[which(mf_all$ub>mf_all$y_max)]
mf_all$lb[which(mf_all$lb<mf_all$y_min)]<-mf_all$y_min[which(mf_all$lb<mf_all$y_min)]



# Plot
zp_all <- ggplot(mf_all, aes(colour = modelName, shape = modelName))
zp_all <- zp_all + facet_grid(. ~ eff, scales="free_x")
zp_all <- zp_all + geom_hline(yintercept = 0, colour = "black", lty = 2)
zp_all <- zp_all + scale_shape_manual(values=c(15,16,17,18,3,4,8)) + scale_size_manual(values=0.7)
zp_all <- zp_all + geom_pointrange(aes(x = Model, y = Coefficient, ymin =lb,
                                       ymax = ub, 0.2),
                                   lwd = 0.7, position = position_dodge(width = 0.6),
                                   fill ="black")
zp_all <- zp_all + coord_flip() + theme_bw()
#zp_all <- zp_all + ylim(xlim[1], xlim[2])
zp_all <- zp_all + geom_blank(aes(y = y_min)) + geom_blank(aes(y = y_max))
zp_all <- zp_all + scale_y_continuous(expand = c(0,0) )
zp_all <- zp_all + theme(legend.title = element_blank())
zp_all <- zp_all + labs(y="Bias", x=element_blank()) 
zp_all <- zp_all + theme(legend.text = element_text(size = 20),
                         legend.position="bottom",
                         legend.key=element_blank(),
                         axis.text.x = element_text(size=16),
                         axis.text.y = element_text(size=16, colour="black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20),
                         plot.title = element_text(size=20, margin = margin(t = 10, b = -19 ), hjust=-0.235),
                         strip.background =element_blank(),
                         strip.text = element_text(size=20, colour="black"),
                         panel.grid.major = element_line(size = (0.5)),panel.grid.minor = element_blank())
zp_all <- zp_all + ggtitle(bquote(rho ~~~~~~ bold(delta) ~~~~~~~~~~~ lambda ~~~~~ bold(theta))) 
zp_all <- zp_all + guides(fill=guide_legend(ncol = 1), colour = guide_legend(override.aes = list(linetype=0), reverse=T),
                          shape = guide_legend(reverse = T))
print(zp_all)




#### Combine plots ####


cairo_pdf(file=paste(od, "sim16_combined_2.pdf", sep=""), width=17.54, height=12.40, 
          bg = "white", family="CM Roman")
par(mar=c(0,0,0,0))
par(mfrow=c(1,1),oma=c(0,0,0,0))
zp_all
dev.off()





