## **************************************************************************
##
##    (c) 2020 Guillaume Guénard
##        Université de Montréal, Montreal, Quebec, Canada
##
##    **Main analysis R script**
##
##    This file is part of HydroMFHM
##
##    HydroMFHM is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    HydroMFHM is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with HydroMFHM.  If not, see <https://www.gnu.org/licenses/>.
##
## **************************************************************************
##
## rm(list=ls())
library(DBI)
library(magrittr)
library(yaml)
library(ape)
library(MPSEM)
library(codep)
library(glmnet)
library(parallel)
source("./Aux/HydroMFHM_Aux.R")
##
### Loading configuration
psqlCon <-
  read_yaml("../config/credentials.yml") %>%
  connectToDB
## DBI::dbListTables(conn=psqlCon)
## disconnectFromDB(psqlCon)
## psqlCon <- connectToDB(cfg,TRUE)
## DBI::dbListTables(conn=psqlCon)
## disconnectFromDB(psqlCon,TRUE)
##
load("../Data/Info.rda")
load("../Data/Flow_metrics.rda")
load("../Data/HNetTrees.rda")
load("../Data/Density.rda")
##
## dgr1 <- TreeHNet[["expanded"]] %>% Phylo2DirectedGraph
PEM1 <-
  Density$fused$density %>%
  log1p %>%
  t %>%
  PEM.forcedSimple(
    x = Density$fused$size %>% as.matrix,
    w = TreeHNet[["expanded"]] %>% Phylo2DirectedGraph
  )
## PEM1$a == 0 : The changes in trait values appear to be neutrally driven
##
## all(rownames(Info[["river"]]) == rownames(Density[["fused"]]$density))
## gcd1 <- Info[["river"]][,c("lat","lon")] %>% gcd.hf2(from=.,to=.)
emap1 <-
  Info[["river"]][,c("lat","lon")] %>%
  gcd.hf2(from=.,to=.) %>%
  eigenmap(opt.coord=Info[["river"]][,c("lat","lon")])
## summary(emap1)
## plot(log10(emap1$lambda),type="l")
##
### Bilinear model data
dat <- list()
##
### Model matrix: among-species descriptors
dat[["Z"]] <- list()
dat[["Z"]]$data <-
  Density$fused$size %>%
  cbind(intercept=1, length=.) %>%
  cbind(PEM1$u) %>%
  as.matrix
## sqrt(colSums(dat[["Z"]]$data^2))
dat[["Z"]]$term <- c("Constant","Trait",rep("Phylogeny",ncol(PEM1$u)))
##
### Model matrix: among-site descriptors
env <- c("log_depth","log1p_velocity","log1p_d50","degree-day","log_tot_phos",
          "macro_covr")
flm <- c("MA3","FH1","TA2","TH2","RA7","nRA1","RL2","MA60")
##
dat[["X"]] <- list()
##
dat[["X"]]$data <-
  Info[["river"]][,env] %>%
  cbind(intercept=1,.) %>%
  cbind(Flow_metrics$values[,flm]) %>%
  cbind(emap1$U) %>%
  as.matrix
## sqrt(colSums(dat[["X"]]$data^2))
dat[["X"]]$term <- c("Constant", rep("Environment",length(env)),
                     rep("Flow",length(flm)), rep("Space",ncol(emap1$U)))
##
### Model matrix: crossed descriptors
## Rows    : nrow(dat[["Z"]]$data)*nrow(dat[["X"]]$data)
## Columns : ncol(dat[["Z"]]$data)*ncol(dat[["X"]]$data)
dat[["blm"]] <- list()
##
dat[["blm"]]$ZkronX <- dat[["Z"]]$data %x% dat[["X"]]$data
## dim(dat[["blm"]]$ZkronX)
dat[["blm"]]$y <- Density$fused$density
dim(dat[["blm"]]$y) <-
  dat[["blm"]]$y %>%
  {c(nrow(.)*ncol(.),1L)}
rownames(dat[["blm"]]$y) <-
  Density$fused$density %>%
  {paste(
    rep(rownames(.),ncol(.)),
    rep(colnames(.),each=nrow(Info$river)),sep=" -> "
  )}
##
dimnames(dat[["blm"]]$ZkronX) <-
  list(
    dat[["X"]]$data %>%
      rownames %>%
      rep(dat[["Z"]]$data %>% nrow) %>%
      paste(
        dat[["Z"]]$data %>%
          rownames %>%
          rep(each=dat[["X"]]$data %>% nrow),
        sep=" -> "
      ),
    dat[["X"]]$data %>%
      colnames %>%
      rep(dat[["Z"]]$data %>% ncol) %>%
      paste(
        dat[["Z"]]$data %>%
          colnames %>%
          rep(each=dat[["X"]]$data %>% ncol),
        sep=" -> "
      )
  )
## dat[["blm"]]$ZkronX %>% colnames %>% head(100L)
##
dat[["blm"]][["mfcol"]] <-
  data.frame(
    X=dat[["X"]]$term %>%
      rep(ncol(dat[["Z"]]$data)) %>%
      as.factor %>%
      .[-1L],
    Z=dat[["Z"]]$term %>%
      rep(each=ncol(dat[["X"]]$data)) %>%
      as.factor %>%
      .[-1L]
  )
##
dat[["blm"]][["mm"]] <-
  cbind(
    contr.treatment(levels(dat$blm$mfcol[,1L]))[dat$blm$mfcol[,1L],],
    contr.treatment(levels(dat$blm$mfcol[,2L]))[dat$blm$mfcol[,2L],2L:1L]
  ) %>%
  cbind(
    .[,1L]*.[,4L], .[,1L]*.[,5L],
    .[,2L]*.[,4L], .[,2L]*.[,5L],
    .[,3L]*.[,4L], .[,3L]*.[,5L]
  )
colnames(dat$blm$mm)[6L:11L] <-
  c("Environment:Trait","Environment:Phylogeny",
    "Flow:Trait","Flow:Phylogeny",
    "Space:Trait","Space:Phylogeny")
rownames(dat$blm$mm) <- NULL
##
## dat$blm$mm %>% head
## dat$blm$mm %>% dim
## dat$blm$mm %>% {11/(1+exp(.%*%cbind(c(-2,1,-5,3,1,2,-6,-2,1,-3,0))))}
### Here cam the idea to use invlogit as a model the dispatch the
### regularization factors (epsilon).
##
### Display the spatial eigenvectors
if(FALSE) {
  rng <- emap1$U %>% range
  cols <- rainbow(1200)[1L:1000L]
  for(i in 1L:ncol(emap1$U)) {
    plot(
      x=Info$river$lon, y=Info$river$lat, pch=21L, main=i, asp=1,
      bg=cols[round(999*(emap1$U[,i] - rng[1L])/(rng[2L] - rng[1L])) + 1]
    )
    if(is.null(locator(1))) break
  }
  rm(rng, cols, i)
}
##
### Display the eigenvalues
if(FALSE) {
  1L:length(emap1$lambda) %>%
    plot(y=emap1$lambda,log="y",type="p",ylab="Eigenvalue")
  l_min <- 400
  abline(v=which(emap1$lambda<l_min)[1L])
  abline(h=emap1$lambda[which(emap1$lambda<l_min)[1L]])
  rm(l_min)
}
##
### Calculation of the Elastic Net regressions models
par <- dat$blm$mm %>% ncol %>% rep(0, .)
dat$blm[["glmnet"]] <-
  dat$blm$ZkronX[,-1L] %>%
  glmnet(
    y = dat$blm$y,
    alpha=0.9,
    penalty.factor=c(0.5,1/(1+exp(-dat$blm$mm%*%par))),
    family="poisson"
  )
## dat$blm$glmnet %>% plot
dat$blm[["glm0"]] <- glm(dat$blm$y~1,family=poisson())
## coef(dat$blm$glmnet,s=0.05) %>% {.!=0} %>% mean %>% {.*100}
##
### Fitted values with a penalty of 0.05
dat$blm[["fitval"]] <-
  coef(dat$blm$glmnet, s=0.05) %>%
  as.matrix %>%
  {dat$blm$ZkronX%*%.} %>%
  exp
dim(dat$blm$fitval) <- dim(Density$fused$density)
dimnames(dat$blm$fitval) <- dimnames(Density$fused$density)
##
{
  rng <- range(Density$fused$density,dat$blm$fitval) %>% log1p
  sp_brace <- 0.5
  for(i in unique(Density$fused$species))
    sp_brace <- c(sp_brace,max(which(Density$fused$species==i))+0.5)
  spnms <- Density$fused$species %>% unique %>% strsplit(" ") %>%
    lapply(
      function(x) {
        if(length(x)==2L) {
          if(x[2L]=="ssp.") x <- paste(x[1L],x[2L],sep=". ")
          else x <- paste(substr(x[1L],1L,1L),x[2L],sep=". ")
        } else x <- paste(substr(x[1L],1L,1L),". ",x[2L]," - ",x[5L],sep="")
        x
      }
    ) %>% unlist
  rnms <- Info$river %>% rownames %>% strsplit('',fixed=TRUE) %>%
    lapply(
      function(x) {
        if(any(x=="(")) x <- x[-(which(x=="("):length(x))]
        for(j in 2L:length(x)) {
          # j <- 4L
          if(any(x[j]==LETTERS)) {
            x <- c(x[1L:j-1L]," ",x[j:length(x)])
            break
          }
        }
        paste(x,collapse="")
      }
    ) %>% unlist
  ##
  cols <- rainbow(1200L)[1L:1000L]
  X11(height=7.25,width=10.0)
  par(fig=c(0.15,0.5,0.15,1),mar=c(2,2,2,1))
  image(z=Density$fused$density %>% log1p,
        x=1L:nrow(Density$fused$density), y=1L:ncol(Density$fused$density),
        zlim=rng, col=cols, axes=FALSE, xlab="", ylab="",
        main="Observed density")
  abline(h=sp_brace)
  axis(1,at=1L:nrow(Density$fused$density),las=2,tick=FALSE,labels=rnms,
       cex.axis=0.6)
  axis(2,at=sp_brace%>%{(1:length(spnms))*(.[1L]+.[length(.)-1])/length(spnms)},
       labels=spnms,las=1,tick=FALSE,cex.axis=0.65,font=3) ; box()
  sp_brace %>% {
    arrows(x0=-0.7,x1=0.25,y0=(1:length(spnms))*(.[1L]+.[length(.)-1L])/length(spnms),
           y1=0.5*(.[1L:(length(.)-1L)]+.[2L:length(.)]),xpd=TRUE,length=0.035)
  }
  par(fig=c(0.5,0.85,0.15,1),mar=c(2,1,2,2),new=TRUE)
  image(z=dat$blm$fitval %>% log1p,
        x=1L:nrow(Density$fused$density), y=1L:ncol(Density$fused$density),
        zlim=rng, col=cols, axes=FALSE, xlab="", ylab="",
        main="Fitted density")
  abline(h=sp_brace)
  axis(1,at=1L:nrow(Density$fused$density),las=2,tick=FALSE,labels=rnms,
       cex.axis=0.6)
  par(fig=c(0.85,1,0.15,1),mar=c(2,0,2,4.5),new=TRUE)
  image(z=t(matrix(seq(rng[1L],rng[2L],length.out=1000L),1000L,10L)),
        x=1L:10L,y=seq(rng[1L],rng[2L],length.out=1000L),
        zlim=rng,col=cols,axes=FALSE,xlab="",ylab="")
  axis(side=4L,at=log1p(c(0,1,2,5,10,18)),label=c(c(0,1,2,5,10,18)),las=1)
  box()
  mtext(side=4L,outer=TRUE,text=expression(fish%.%(100~m^{2})^{-1}),line=-2)
  rm(rng,sp_brace,spnms,rnms,i)
}
dev.copy2eps(file="../Image/Preliminary model fit.eps")
dev.off()
##
### A temptative coefficient of determination based on likelihood
{
  dat$blm[["metrics"]] <- c()
  ##
  dat$blm$metrics["L_perfect"] <-
    Density$fused$density %>%
    as.numeric %>%
    PoissonPDF_log(.,lambda=.) %>%
    sum
  ##
  dat$blm$metrics["L_model"] <-
    Density$fused$density %>%
    as.numeric %>%
    PoissonPDF_log(lambda=dat$blm$fitval%>%as.numeric) %>%
    sum
  ##
  dat$blm$metrics["L_null"] <-
    Density$fused$density %>%
    as.numeric %>%
    PoissonPDF_log(lambda=dat$blm$glm0$coef%>%exp) %>%
    sum
  ##
  dat$blm$metrics["ML_Rsqr"] <- dat$blm$metrics %>%
    {1-(.["L_perfect"] - .["L_model"])/(.["L_perfect"]-.["L_null"])}
}
##
### Explanatory figure for the likelihood-based determination coefficient
{
  X11(width=5.0,height=3.5)
  par(mar=c(4.5,4.5,3.5,8))
  dat$blm$metrics %>%
    {plot(y=1-(.["L_perfect"]-seq(.["L_perfect"],1.25*.["L_null"],length.out=2))/(.["L_perfect"]-.["L_null"]),
          x=seq(.["L_perfect"],1.25*.["L_null"],length.out=2), type="l",
          ylim=c(-0.25,1.1), ylab=expression(italic(R)^2),
          xlab="Log Likelihood",xlim=c(-2500,-200),las=1)}
  
  abline(h=1,col=gray(0.7))
  dat$blm$metrics %>% {abline(v=.["L_perfect"],col=gray(0.7))}
  dat$blm$metrics %>% {points(y=1, x=.["L_perfect"], pch=21L, bg="black")}
  abline(h=0,col=gray(0.7))
  dat$blm$metrics %>% {abline(v=.["L_null"],col=gray(0.7))}
  dat$blm$metrics %>% {points(y=0, x=.["L_null"], pch=21L, bg="black")}
  dat$blm$metrics %>%
    {abline(h=1-(.["L_perfect"]-.["L_model"])/(.["L_perfect"]-.["L_null"]),
            col=gray(0.5))}
  dat$blm$metrics %>% {abline(v=.["L_model"], col=gray(0.5))}
  dat$blm$metrics %>%
    {points(y=1-(.["L_perfect"]-.["L_model"])/(.["L_perfect"]-.["L_null"]),
            x=.["L_model"], pch=21L, bg="black")}
  dat$blm$metrics %>%
    {axis(3L,at=c(.["L_null"],.["L_model"],.["L_perfect"]),
          labels=paste(c("null","Model","«perfect»"),"\nlikelihood"))}
  dat$blm$metrics %>% {
    axis(
      4L,
      at=c(1,1-(.["L_perfect"]-.["L_model"])/(.["L_perfect"]-.["L_null"]),0),
      label=c("Perfect model","Useful model","Threshold value"),
      las=1
    )}
}
##
dev.copy2eps(file="../Image/Likelihood-based Determination coefficient.eps")
dev.off()
##
### Splitting the cross-validation data into groups of sites and species
cvdat <- list()
##
### Species groups
cvdat[["sp"]] <- list(n=4L)
cvdat$sp[["grp"]] <-
  cvdat$sp$n %>% {1L:.} %>%
  rep(length.out=Density$fused$species %>% unique %>% length) %>%
  LETTERS[.]
names(cvdat$sp$grp) <-
  PEM1$u[,1L] %>%
  tapply(Density$fused$species,head,n=1) %>%
  sort %>%
  names
cvdat$sp[["ss"]] <-
  Density$fused$species %>%
  match(cvdat$sp$grp %>% names) %>%
  cvdat$sp$grp[.]
##
### Site groups
cvdat[["st"]] <- list(n=4L)
cvdat$st[["grp"]] <-
  cvdat$st$n %>% {1L:.} %>%
  rep(length.out=Info$river %>% nrow) %>%
  LETTERS[.]
names(cvdat$st$grp) <-
  Info$river %>% rownames
##
## Info$river$flow_regime %>% as.character %>% tapply(cvdat$st$grp,unique)
##
### Cross-validation groups
cvdat[["groups"]] <- list()
for(i in 1L:cvdat$sp$n) for(j in 1L:cvdat$st$n) {
  ## i <- j <- 1L
  lab <- paste("Sp_",LETTERS[i],"-St_",LETTERS[j],sep="")
  ii <- cvdat$sp$ss!=LETTERS[i]
  jj <- cvdat$st$grp!=LETTERS[j]
  cvdat$groups[[lab]] <- list()
  cvdat$groups[[lab]]$y <-
    Density$fused$density[jj,ii]
  dim(cvdat$groups[[lab]]$y) <- c(cvdat$groups[[lab]]$y %>% length, 1L)
  emap_tmp <-
    Info$river[jj,c("lat","lon")] %>%
    gcd.hf2(from=.,to=.) %>%
    eigenmap(opt.coord=Info$river[jj,c("lat","lon")])
  X_tmp <-
    dat$X %>%
    {.$data[jj,.$term!="Space"]} %>%
    cbind(emap_tmp$U) %>%
    list(data=.)
  X_tmp$term <-
    dat$X$term %>%
    {.[.!="Space"]} %>%
    c(rep("Space",ncol(emap_tmp$U)))
  Xp_tmp <-
    dat$X %>%
    {.$data[!jj,.$term!="Space"]} %>%
    cbind(
      emap_tmp %>%
        eigenmap.score(
          Info$river[,c("lat","lon")] %>%
            {gcd.hf2(from=.[jj,],to=.[!jj,])}
        )
    )
  grloc_tmp <-
    TreeHNet$expanded %>%
    getGraphLocations(
      target=.$tip.label[!ii]
    )
  PEM_tmp <-
    PEM.forcedSimple(
      y=Density$fused$density[jj,ii] %>% log1p %>% t,
      x=Density$fused$size[ii] %>%
        as.matrix,
      w=grloc_tmp$x
    )
  Z_tmp <-
    dat$Z$data[ii,1L:2L] %>%
    cbind(PEM_tmp$u) %>%
    list(data=.)
  Z_tmp$term <-
    dat$Z$term %>%
    {.[.!="Phylogeny"]} %>%
    c(rep("Phylogeny",ncol(PEM_tmp$u)))
  Zp_tmp <-
    dat$Z %>%
    {.$data[!ii,.$term!="Phylogeny"]} %>%
    cbind(
      PEM_tmp %>%
        Locations2PEMscores(grloc_tmp) %>%
        {.$scores}
    )
  cvdat$groups[[lab]]$descriptors <- Z_tmp$data%x%X_tmp$data
  cvdat$groups[[lab]]$predictors <- Zp_tmp%x%Xp_tmp
  mfcol_tmp <-
    data.frame(
      X=X_tmp$term %>%
        rep(Z_tmp$data %>% ncol),
      Z=Z_tmp$term %>%
        rep(each=X_tmp$data %>% ncol)
    )[-1L,]
  mm_tmp <-
    cbind(
      mfcol_tmp[,1L] %>%
        levels %>%
        contr.treatment %>%
        .[mfcol_tmp[,1L],],
      mfcol_tmp[,2L] %>%
        levels %>%
        contr.treatment %>%
        .[mfcol_tmp[,2L],2L:1L]
    )
  mm_tmp %<>%
    cbind(
      mm_tmp[,1L]*mm_tmp[,4L],
      mm_tmp[,1L]*mm_tmp[,5L],
      mm_tmp[,2L]*mm_tmp[,4L],
      mm_tmp[,2L]*mm_tmp[,5L],
      mm_tmp[,3L]*mm_tmp[,4L],
      mm_tmp[,3L]*mm_tmp[,5L]
    )
  colnames(mm_tmp)[6L:11] <-
    c("Environment:Trait","Environment:Phylogeny",
      "Flow:Trait","Flow:Phylogeny",
      "Space:Trait","Space:Phylogeny")
  rownames(mm_tmp) <- NULL
  cvdat$groups[[lab]]$mfcol <- mfcol_tmp
  cvdat$groups[[lab]]$mm <- mm_tmp
  cvdat$groups[[lab]]$obs_y <-
    Density$fused$density[!jj,!ii]
  dim(cvdat$groups[[lab]]$obs_y) <-
    c(cvdat$groups[[lab]]$obs_y %>% length, 1L)
}
rm(i,j,lab,ii,jj,emap_tmp,X_tmp,Xp_tmp,grloc_tmp,PEM_tmp,Z_tmp,Zp_tmp,
   mfcol_tmp,mm_tmp)
##
cvdat[["obs"]] <-
  cvdat$groups %>%
  lapply(function(x) x$obs_y) %>%
  unlist
##
cvdat[["mean"]] <-
  cvdat$groups %>%
  lapply(
    function(x) {
      tmp <- glm(x$y~1, family=poisson())
      (x$predictors[,1L,drop=FALSE]%*%tmp$coef) %>% exp
    }
  ) %>% unlist
##
cvdat[["metrics"]] <- c()
##
### The likelihood for the perfect predictions has been calculated previously
cvdat$metrics["L_perfect"] <-
  dat$blm$metrics["L_perfect"]
##
### The likelihood of mean-only models calculated for individual CV groups
cvdat$metrics["L_null"] <-
  cvdat[["obs"]] %>%
  PoissonPDF_log(lambda=cvdat[["mean"]]) %>%
  sum
##
save(cvdat, file = "../Data/cvdat.rda")
##
### Values for the different cv groups must be calculated in parallel
if(FALSE) {
  ## cl <- detectCores() %>% makeForkCluster
  cl <- makeForkCluster(4L)
  glm <- cvdat$groups %>% realpha(cl, alpha=0.9, penalty.coef=rep(0,11L))
  val <-
    objf_lambda(
      par=-1,cl=cl,grp=cvdat$groups,glm=glm,obs=cvdat$obs,
      L_perfect=cvdat$metrics["L_perfect"],L_null_CV=cvdat$metrics["L_null"]
    )
  ##
  rm(glm,val)
  ##
  res <-
    objf_alpha(
      par=c(0,rep(0,11L)),cl=cl,grp=cvdat$groups,obs=cvdat$obs,
      L_perfect=cvdat$metrics["L_perfect"],L_null_CV=cvdat$metrics["L_null"]
    )
  rm(res)
  ## cl %>% stopCluster
  ## rm(cl)
}
##
if(FALSE) {
  cl <- makeForkCluster(4L)
  opt_alpha <-
    optim(
      par=rep(0,12L),f=objf_alpha,method="BFGS",cl=cl,grp=cvdat$groups,
      obs=cvdat$obs,L_perfect=cvdat$metrics["L_perfect"],
      L_null_CV=cvdat$metrics["L_null"],control=list(maxit=1000)
    )
  glm <- cvdat$groups %>%
    realpha(
      cl=cl,
      alpha=(1+exp(-opt_alpha$par[1L]))^-1,
      penalty.coef=opt_alpha$par[-1L]
    )
  opt_lambda <-
    optim(
      par=-5,f=objf_lambda,method="Brent",cl=cl,grp=cvdat$groups,glm=glm,
      obs=cvdat$obs,L_perfect=cvdat$metrics["L_perfect"],
      L_null_CV=cvdat$metrics["L_null"],control=list(maxit=1000),
      lower=-10,upper=0
    )
  save(opt_alpha,glm,opt_lambda,file="../Data/Regularization parameters.rda")
  cl %>% stopCluster
  rm(cl)
} else load(file="../Data/Regularization parameters.rda")
##
if(FALSE) {
  cl <- makeForkCluster(4L)
  opt_alpha_sann <-
    optim(
      par=rep(0,12L),f=objf_alpha,method="SANN",cl=cl,grp=cvdat$groups,
      obs=cvdat$obs,L_perfect=cvdat$metrics["L_perfect"],
      L_null_CV=cvdat$metrics["L_null"],control=list(maxit=1000)
    )
  glm_sann <- cvdat$groups %>%
    realpha(
      cl=cl,
      alpha=(1+exp(-opt_alpha_sann$par[1L]))^-1,
      penalty.coef=opt_alpha_sann$par[-1L]
    )
  opt_lambda_sann <-
    optim(
      par=-5,f=objf_lambda,method="Brent",cl=cl,grp=cvdat$groups,glm=glm_sann,
      obs=cvdat$obs,L_perfect=cvdat$metrics["L_perfect"],
      L_null_CV=cvdat$metrics["L_null"],control=list(maxit=1000),
      lower=-10,upper=0
    )
  save(opt_alpha_sann,glm_sann,opt_lambda_sann,
       file="../Data/Regularization parameters (sann).rda")
  cl %>% stopCluster
  rm(cl)
} else load(file="../Data/Regularization parameters (sann).rda")
### Global optimization yielded a worst solution than gradient descent.
##
cvdat[["ptype"]] <-
  apply(dat$blm$mm,1L,paste,collapse=" ") %>%
  unique %>%
  sort(decreasing=FALSE) %>%
  strsplit(" ") %>%
  sapply(as.numeric) %>% t
colnames(cvdat[["ptype"]]) <-
  colnames(dat$blm$mm)
##
cvdat[["cont"]] <-
  dat$blm$mfcol[,1L]:dat$blm$mfcol[,2L]
##
swap <-
  c(`Constant:Constant`="(Intercept)",`Constant:Phylogeny`="Phylogeny",
    `Constant:Trait`="Trait",`Environment:Constant`="Environment",
    `Flow:Constant`="Flow",`Space:Constant`="Space")
lvl <-
  cvdat$cont %>%
  levels %>%
  sapply(function(x,swap) if(x %in% names(swap)) swap[x] else x, swap=swap)
##
cvdat$cont %<>%
  as.numeric %>%
  {lvl[.]} %>%
  as.factor %<>%
  tapply(.,.,length)
##
rownames(cvdat$ptype) <-
  cvdat$ptype %>%
  apply(
    1L,
    function(x,cn) {
      cn[which(!!x)[length(which(!!x))]]
    },
    cn=colnames(.)
  )
##
rm(swap,lvl)
##
"%s: %f\n" %>%
  sprintf(
    rownames(cvdat$ptype),
    exp(opt_lambda$par)*(1+exp(-cvdat$ptype%*%opt_alpha$par[-1L]))^-1
  ) %>%
  cat
"%s: %f\n" %>%
  sprintf(
    rownames(cvdat$ptype),
    as.numeric(exp(opt_lambda$par)*(1+exp(-cvdat$ptype%*%opt_alpha$par[-1L]))^-1)/
      cvdat$cont[rownames(cvdat$ptype)]
  ) %>%
  cat
##
cvdat[["regularization"]] <- list()
cvdat$regularization[["alpha"]] <-
  (1+exp(-opt_alpha$par[1L]))^-1
cvdat$regularization[["lambda"]] <-
  exp(opt_lambda$par)
cvdat$regularization[["penalty.coef"]] <-
  opt_alpha$par[-1L]
names(cvdat$regularization[["penalty.coef"]]) <-
  colnames(dat$blm$mm)
cvdat$regularization[["penalty.factor"]] <-
  rbind(`(Intercept)`=0.5,(1+exp(-cvdat$ptype%*%opt_alpha$par[-1L]))^-1)
cvdat$regularization[["total.penalty"]] <-
  cvdat$regularization[["penalty.factor"]] %>%
  {opt_lambda$par %>% exp*.}
##
dat[["reg"]] <- list()
dat$reg[["glmnet"]] <-
  dat$blm$ZkronX[,-1L] %>%
  glmnet(
    y=dat$blm$y,
    alpha=cvdat$regularization$alpha,
    penalty.factor=c(0.5,(1+exp(-dat$blm$mm%*%opt_alpha$par[-1L]))^-1),
    family="poisson"
  )
## dat$reg[["glmnet"]] %>% plot
dat$reg[["coefficients"]] <-
  dat$reg[["glmnet"]] %>%
  coef(
    s=exp(opt_lambda$par),
    exact=TRUE,
    x=dat$blm$ZkronX[,-1L],
    y=dat$blm$y,
    penalty.factor=c(0.5,(1+exp(-dat$blm$mm%*%opt_alpha$par[-1L]))^-1)
  ) %>% as.matrix
## dat$reg$coefficients %>% head(100L)
##
dat$reg[["fitval"]] <-
  dat$blm$ZkronX%*%dat$reg$coefficients %>%
  exp
dim(dat$reg[["fitval"]]) <-
  dim(Density$fused$density)
dimnames(dat$reg[["fitval"]]) <-
  dimnames(Density$fused$density)
##
{
  rng <- range(Density$fused$density,dat$reg$fitval) %>% log1p
  sp_brace <- 0.5
  for(i in unique(Density$fused$species))
    sp_brace <- c(sp_brace,max(which(Density$fused$species==i))+0.5)
  spnms <- Density$fused$species %>% unique %>% strsplit(" ") %>%
    lapply(
      function(x) {
        if(length(x)==2L) {
          if(x[2L]=="ssp.") x <- paste(x[1L],x[2L],sep=". ")
          else x <- paste(substr(x[1L],1L,1L),x[2L],sep=". ")
        } else x <- paste(substr(x[1L],1L,1L),". ",x[2L]," - ",x[5L],sep="")
        x
      }
    ) %>% unlist
  rnms <- Info$river %>% rownames %>% strsplit('',fixed=TRUE) %>%
    lapply(
      function(x) {
        if(any(x=="(")) x <- x[-(which(x=="("):length(x))]
        for(j in 2L:length(x)) {
          # j <- 4L
          if(any(x[j]==LETTERS)) {
            x <- c(x[1L:j-1L]," ",x[j:length(x)])
            break
          }
        }
        paste(x,collapse="")
      }
    ) %>% unlist
  ##
  cols <- rainbow(1200L)[1L:1000L]
  X11(height=7.25,width=10.0)
  par(fig=c(0.15,0.5,0.15,1),mar=c(2,2,2,1))
  image(z=Density$fused$density %>% log1p,
        x=1L:nrow(Density$fused$density), y=1L:ncol(Density$fused$density),
        zlim=rng, col=cols, axes=FALSE, xlab="", ylab="",
        main="Observed density")
  abline(h=sp_brace)
  axis(1,at=1L:nrow(Density$fused$density),las=2,tick=FALSE,labels=rnms,
       cex.axis=0.6)
  axis(2,at=sp_brace%>%{(1:length(spnms))*(.[1L]+.[length(.)-1])/length(spnms)},
       labels=spnms,las=1,tick=FALSE,cex.axis=0.65,font=3) ; box()
  sp_brace %>% {
    arrows(x0=-0.7,x1=0.25,y0=(1:length(spnms))*(.[1L]+.[length(.)-1L])/length(spnms),
           y1=0.5*(.[1L:(length(.)-1L)]+.[2L:length(.)]),xpd=TRUE,length=0.035)
  }
  par(fig=c(0.5,0.85,0.15,1),mar=c(2,1,2,2),new=TRUE)
  image(z=dat$reg$fitval %>% log1p,
        x=1L:nrow(Density$fused$density), y=1L:ncol(Density$fused$density),
        zlim=rng, col=cols, axes=FALSE, xlab="", ylab="",
        main="Fitted density")
  abline(h=sp_brace)
  axis(1,at=1L:nrow(Density$fused$density),las=2,tick=FALSE,labels=rnms,
       cex.axis=0.6)
  par(fig=c(0.85,1,0.15,1),mar=c(2,0,2,4.5),new=TRUE)
  image(z=t(matrix(seq(rng[1L],rng[2L],length.out=1000L),1000L,10L)),
        x=1L:10L,y=seq(rng[1L],rng[2L],length.out=1000L),
        zlim=rng,col=cols,axes=FALSE,xlab="",ylab="")
  axis(side=4L,at=log1p(c(0,1,2,5,10,18)),label=c(c(0,1,2,5,10,18)),las=1)
  box()
  mtext(side=4L,outer=TRUE,text=expression(fish%.%(100~m^{2})^{-1}),line=-2)
  rm(rng,sp_brace,spnms,rnms,i)
}
dev.copy2eps(file="../Image/Regularized model fit.eps")
dev.off()
##
{
  dat$reg[["metrics"]] <- c()
  ##
  dat$reg$metrics["L_perfect"] <-
    dat$blm$metrics["L_perfect"]
  ##
  dat$reg$metrics["L_model"] <-
    Density$fused$density %>%
    as.numeric %>%
    PoissonPDF_log(lambda=dat$reg$fitval%>%as.numeric) %>%
    sum
  ##
  dat$reg$metrics["L_null"] <-
    dat$blm$metrics["L_null"]
  ##
  dat$reg$metrics["ML_Rsqr"] <- dat$reg$metrics %>%
    {1-(.["L_perfect"] - .["L_model"])/(.["L_perfect"]-.["L_null"])}
}
## dat$blm$metrics
## dat$reg$metrics
##

### Explanatory figure for the likelihood-based determination coefficient
{
  X11(width=5.0,height=3.5)
  par(mar=c(4.5,4.5,3.5,8))
  dat$reg$metrics %>%
    {plot(y=1-(.["L_perfect"]-seq(.["L_perfect"],1.25*.["L_null"],length.out=2))/(.["L_perfect"]-.["L_null"]),
          x=seq(.["L_perfect"],1.25*.["L_null"],length.out=2), type="l",
          ylim=c(-0.25,1.1), ylab=expression(italic(R)^2),
          xlab="Log Likelihood",xlim=c(-2500,-200),las=1)}
  
  abline(h=1,col=gray(0.7))
  dat$reg$metrics %>% {abline(v=.["L_perfect"],col=gray(0.7))}
  dat$reg$metrics %>% {points(y=1, x=.["L_perfect"], pch=21L, bg="black")}
  abline(h=0,col=gray(0.7))
  dat$reg$metrics %>% {abline(v=.["L_null"],col=gray(0.7))}
  dat$reg$metrics %>% {points(y=0, x=.["L_null"], pch=21L, bg="black")}
  dat$reg$metrics %>%
    {abline(h=1-(.["L_perfect"]-.["L_model"])/(.["L_perfect"]-.["L_null"]),
            col=gray(0.5))}
  dat$reg$metrics %>% {abline(v=.["L_model"], col=gray(0.5))}
  dat$reg$metrics %>%
    {points(y=1-(.["L_perfect"]-.["L_model"])/(.["L_perfect"]-.["L_null"]),
            x=.["L_model"], pch=21L, bg="black")}
  dat$reg$metrics %>%
    {axis(3L,at=c(.["L_null"],.["L_model"],.["L_perfect"]),
          labels=paste(c("null","Model","«perfect»"),"\nlikelihood"))}
  dat$reg$metrics %>% {
    axis(
      4L,
      at=c(1,1-(.["L_perfect"]-.["L_model"])/(.["L_perfect"]-.["L_null"]),0),
      label=c("Perfect model","Useful model","Threshold value"),
      las=1
    )}
}
##
dev.copy2eps(file="../Image/Likelihood-based Determination coefficient.eps")
dev.off()
##
dat$reg[["B"]] <-
  dat$reg$coefficients
dim(dat$reg$B) <- c(ncol(dat$X$data),ncol(dat$Z$data))
dimnames(dat$reg$B) <-
  list(colnames(dat$X$data),colnames(dat$Z$data))
## 100*mean(dat$reg$B!=0)
##
dat$reg[["sel_cont"]] <-
  matrix(
    NA,
    dat$X$term%>%unique%>%length,
    dat$Z$term%>%unique%>%length,
    dimnames=list(dat$X$term%>%unique,dat$Z$term%>%unique)
  )
dat$reg[["tot_cont"]] <-
  dat$reg[["sel_cont"]]
for(i in dat$X$term%>%unique)
  for(j in dat$Z$term%>%unique) {
    dat$reg$sel_cont[i,j] <- (!!dat$reg$B[dat$X$term==i,dat$Z$term==j])%>%sum
    dat$reg$tot_cont[i,j] <- sum(dat$X$term==i)*sum(dat$Z$term==j)
  }
rm(i,j)
## dat$reg$sel_cont
## dat$reg$tot_cont
## dat$reg %>% {100*.$sel_cont/.$tot_cont}
dat %>% Apply_nonzero_B("Environment","Phylogeny",colSums)
dat %>% Apply_nonzero_B("Environment","Phylogeny",rowSums)
dat %>% Apply_nonzero_B("Flow","Phylogeny",colSums)
dat %>% Apply_nonzero_B("Flow","Phylogeny",rowSums)
dat %>% Apply_nonzero_B("Space","Phylogeny",which)
##
### The intercept of the model
dat$reg$B[1L,1L,drop=FALSE]
##
### Marginal effects of the environmental variables
dat %>% list_which_nonzero("Environment","Constant")
## Degree-day         (+ effect)
## Total phosphorus   (+ effect)
##
### Marginal effect of the flow descriptors
dat %>% list_which_nonzero("Flow","Constant")
## MA3: variability in daily flow (+ effect)
##
### Marginal effect of space
dat %>% list_which_nonzero("Space","Constant")
## No effect was found
##
### Marginal effect of fish length
dat %>% list_which_nonzero("Constant","Trait",TRUE)
## No effect found
##
### Marginal effect of the phylogeny
dat %>% list_which_nonzero("Constant","Phylogeny",TRUE)
## V1 acting on its own
##
### Environment-length interactions
dat %>% list_which_nonzero("Environment","Trait")
## Effect of median substrate size on density descreases with fish size
## Effect of degree-day on density descreases with fish size
## Effect of total phosphorus on density decrease with fish size
##
### Flow-length interactions
dat %>% list_which_nonzero("Flow","Trait")
## Effect of the variability in daily flow on density decreases with fish size
##
par(mar=c(5,7,2,2))
dat %>%
  interaction_display(
    mod="reg", X=c("log1p_d50","degree-day","log_tot_phos","MA3"),
    labels=c("log(D50+1)","Degree-day","log(Total P)","MA3"),
    Z="length", xlab="Effect of fish length on density",
    signif=c(2L,0L,3L,1L)
  )
##
dev.copy2eps(file="../Image/Environment on the Length-Density relationship.eps")
dev.off()
##
### Space-length interactions
dat %>% list_which_nonzero("Space","Trait")
## A total of 11 terms.
## The absence of a margin effect for length combined with the presence of so
## many interaction terms entails that effect of length on density is markedly
## spatially heterogeneous.
##
### Displaying these marginal effects and interactions for description purposes
##
### Data to help display the spatial eigenfunction for each eigenvector
if(FALSE) {
  space_grid <- list()
  space_grid[["gcd"]] <-
    Info$river[,c("lat","lon")] %>%
    gcd.hf2(.,.)
  for(i in c("lat","lon"))
    space_grid[[i]] <-
    Info$river[,i] %>%
    {seq(floor(min(.)*10)/10-1.5,
         ceiling(max(.)*10)/10+1.5,0.1)}
  space_grid[["crd"]] <-
    space_grid %>% {
      cbind(
        Lat=rep(.$lat,length(.$lon)),
        Lon=rep(.$lon,each=length(.$lat))
      )}
  space_grid[["gcdscan"]] <-
    space_grid %>% {
      gcd.hf2(
        from=Info$river[,c("lat","lon")],
        to=.$crd
      )}
  space_grid[["emapscan"]] <-
    emap1 %>% eigenmap.score(space_grid[["gcdscan"]])
  save(space_grid, file="../Data/space_grid.rda")
}  else load(file="../Data/space_grid.rda")
##
##
X11(height=3.4, width=12.9)
## dev.size("in")
par(mar=c(5,5,2,2))
Space_display(
  dat, space_grid, col=rainbow(1200)[1L:1000L],
  mod="reg", Z="length", digits=2L, title="Total length",
  xlab="Longitude", ylab="Latitude", xleg=-110, yleg=51.5
)
##
dev.copy2eps(file="../Image/Length-Density Spatial variation.eps")
dev.off()
##
### Environment-phylogeny interactions
dat %>% list_which_nonzero("Environment","Phylogeny")
dat %>% list_which_nonzero("Environment","Phylogeny",TRUE)
## Some environmental variables, namely current velocity, median substrate grain
## size, degree-day, total phosphorus, and macrophyte cover, affect fish density
## differently depending on their phylogenetic patterns of parenthood. 
##
## In addition to the (+) marginal effect of the variability in daily flow on
## fish density, that same variable affect species densities differently
## depending on their phylogenetic making.
##
### Flow-phylogeny interactions
dat %>% list_which_nonzero("Flow","Phylogeny")
##
if(FALSE) {
  phylo_grid <- list()
  ##
  phylo_grid[["tree"]] <- TreeHNet$raw
  ## rownames(PEM1$u)==colnames(Density$fused$density)
  phylo_grid[["pem"]] <-
    Density$fused$species %>%
    match(unique(.),.) %>%
    PEM1$u[.,]
  ## rownames(phylo_grid$pem) %>% cbind(phylo_grid$tree$tip.label)
  rownames(phylo_grid$pem) <-
    phylo_grid$tree$tip.label
  save(phylo_grid,file="../Data/phylo_grid.rda")
} else load(file="../Data/phylo_grid.rda")
##
X11(width=5.0,height=6.5)
par(mar=c(4,1,2,1), fig = c(0,1,0,1))
dat %>%
  phylo_display(
    phylo_grid=phylo_grid, mod="reg", X="intercept",
    labels="", phylo.width=0.7,
    cex.phylo=0.85, bg="grey", std=FALSE
  )
dev.copy2eps(file="../Image/Density phylogenetic variation (marginal).eps")
dev.off()
##
X11(width=9.0,height=7.25)
par(mar=c(4,1,3,1), fig = c(0,1,0,1))
##
dat %>%
  phylo_display(
    phylo_grid=phylo_grid, mod="reg",
    X=c("log1p_velocity","log1p_d50","degree-day"),
    labels=c("log(Velocity+1)","log(D50+1)","Degree-day"),
    cex.phylo=0.85, bg="grey"
  )
dev.copy2eps(file="../Image/Environment-Density phylogenetic variation 1.eps")
##
dat %>%
  phylo_display(
    phylo_grid=phylo_grid, mod="reg",
    X=c("log_tot_phos","macro_covr","MA3"),
    labels=c("log(Total P)","Macrophyte","Flow variability\n(MA3)"),
    cex.phylo=0.85, bg="grey"
  )
dev.copy2eps(file="../Image/Environment-Density phylogenetic variation 2.eps")
dev.off()
##
## Space-phylogeny interactions
dat %>% list_which_nonzero("Space","Phylogeny")
## Very few to be told here: a single term was highlighted as useful at making
## predictions.
##
X11(height=3.4, width=12.9)
## dev.size("in")
par(mar=c(5,5,2,2))
Space_display(
  dat, space_grid, col=rainbow(1200)[1L:1000L],
  mod="reg", Z="V_12", digits=2L, title="PEM: V_12",
  xlab="Longitude", ylab="Latitude", xleg=-110, yleg=51.5
)
##
X11(width=9.0,height=7.25)
par(mar=c(4,1,3,1), fig = c(0,1,0,1))
phylo_display(
  dat, phylo_grid=phylo_grid, mod="reg",X="dbMEM21",
  cex.phylo=0.85, bg="grey"
)
##
### Discussion
##






##
disconnectFromDB(psqlCon)
rm(psqlCon)
##
