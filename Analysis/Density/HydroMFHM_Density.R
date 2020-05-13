## **************************************************************************
##
##    (c) 2020 Guillaume Guénard
##        Université de Montréal, Montreal, Quebec, Canada
##
##    **River densities R script**
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
source("../Aux/HydroMFHM_Aux.R")
##
### Loading configuration
psqlCon <-
  read_yaml("../../config/credentials.yml") %>%
  connectToDB
## DBI::dbListTables(conn=psqlCon)
## disconnectFromDB(psqlCon)
## psqlCon <- connectToDB(cfg,TRUE)
## DBI::dbListTables(conn=psqlCon)
## disconnectFromDB(psqlCon,TRUE)
##
load("../../Data/Info.rda")
load(file="../../Data/AllTrees.rda")
Taxonomy <-
  "../../Data/BCF-Taxonomy.csv" %>%
  read.csv(stringsAsFactor=FALSE)
SpeciesL <-
  "../../Data/Species List 2010_2013.csv" %>%
  read.csv(as.is=TRUE) %>%
  {.[!is.na(.[,"Tree.name"]),]}
##
Density <- list()
##
Density[["site"]] <-
  "../../Data/densites_SPCL_sites.csv" %>%
  read.csv(as.is=TRUE, sep=";", dec=",")
rownames(Density[["site"]]) <- Density[["site"]][,"Site"]
Density[["site"]] %<>% {.[rownames(Info[["site"]]),-(1L:9L)]}
## nrow(Density) ; nrow(Info[["site"]])
## 100*sum(colSums(Density[["site"]])<=sqrt(.Machine$double.eps))/ncol(Density[["site"]])
## 23.98754% of species/size class have never been observed.
##
### Grouping species X size classes by rivers
##
### First: getting rid of the species X size classes with no contingence
Density[["site"]] %<>% {.[,colSums(.)>sqrt(.Machine$double.eps)]}
## Density[["site"]] %>% apply(2L,range) %>% t
##
### Second: grouping fish densities by rivers
Density[["river"]] <-
  matrix(NA,nrow(Info[["river"]]),ncol(Density[["site"]]),
         dimnames=list(rownames(Info[["river"]]),
                       colnames(Density[["site"]])))
for(i in 1L:ncol(Density[["site"]])) {
  ## i <- 1L
  for(j in rownames(Info[["river"]])) {
    ## j <- rownames(Info[["river"]])[1L]
    dd <- Density[["site"]][Info[["site"]][,"river"]==j,i]
    if(sum(dd)<sqrt(.Machine$double.eps)) {
      Density[["river"]][j,i] <- 0
    } else {
      ee <- data.frame(
        log_depth=log(Info[["site"]][Info[["site"]][,"river"]==j,"depth"]),
        log1p_velocity=log1p(Info[["site"]][Info[["site"]][,"river"]==j,"velocity"]),
        log1p_d50=log1p(Info[["site"]][Info[["site"]][,"river"]==j,"d50"]))
      Density[["river"]][j,i] <-
        glm(dd~.,data=ee,family=poisson) %>%
        predict(newdata=Info[["river"]][j,c("log_depth","log1p_velocity","log1p_d50")]) %>%
        exp
    }
  }
}
rm(i,j,dd,ee)
## Density[["river"]] %>% apply(2L,range) %>% t
##
### Third: merging species and size class information
tmp <- Density[["river"]] %>% colnames %>% strsplit("_")
Density[["spszcl"]] <-
  tmp %>% length %>%
  {list(species=character(.), size=numeric(.))}
for(i in 1L:length(tmp)) {
  Density[["spszcl"]]$species[i] <- tmp[[i]][1L]
  Density[["spszcl"]]$size[i] <- as.numeric(tmp[[i]][2L])
}
rm(tmp,i)
##
Density[["sizeclasses"]] <-
  Density[["spszcl"]]$size %>%
  {data.frame(from=c(3,5*(1:(max(.)-1L))),to=5*(1:max(.)))}
Density[["sizeclasses"]][,"mean"] <-
  Density[["sizeclasses"]][,1L:2L] %>% rowMeans
## Density[["spszcl"]]$size %>% unique
## Density[["spszcl"]]$species %>% unique
Density[["spszcl"]]$species %<>%
  {SpeciesL[match(.,SpeciesL[,"Concatenated.name"]),"Tree.name"]}
colnames(Density[["river"]]) <-
  Density[["spszcl"]] %>%
  {paste(.$species, "-", .$size)}
Density[["river"]] %<>%
  as.matrix %>% {100*.}   ## Density in fish /100 m²
##
## Density[["river"]] %>% summary
## Density[["river"]] %>% as.numeric %>% density %>% plot
## Density[["river"]] %>% as.numeric %>% density(from=0,to=0.6) %>% plot(xlim=c(0,0.6))
## Density[["river"]] %>% as.numeric %>% sort %>% tail(n=100L)
## Density[["river"]] %>% {100*sum(as.numeric(.)==0)/(nrow(.)*ncol(.))}
## 83.5041% of observation have 0 density (once grouped using a Poisson GLM)
## Used to be 83.41628% while using a Gaussian GLM (Ecosphere paper)
##
### Fourth: to stand a chance to model anything at all (since it is utterly
###         zero-inflated), we will fuse the size classes that too little
###         contingence, separately for each species.
mincont <- 10   ## Set minimum contingency == 10%
Density[["fused"]] <- list()
##
for (sp in unique(Density[["spszcl"]]$species)) {
  ## sp <- unique(Density[["spszcl"]]$species)[1L]
  ## sp <- "Petromyzon marinus"
  sc <- Density[["spszcl"]] %>%
    {Density[["sizeclasses"]][.$size[.$species==sp],3L]}
  tmp <- Density %>%
    {.[["river"]][,.[["spszcl"]]$species==sp,drop=FALSE]}
  cont <- tmp %>% {100*colSums(.!=0)/nrow(.)} ## Percent contingency
  ##
  while(length(cont)>1L) {
    i <- which.min(cont)
    if(cont[i] < mincont) {
      ii <- c(i-1L,i+1L) ; ii <- ii[ii>=1L] ; ii <- ii[ii<=length(cont)]
      if(length(ii)!=1L) ii <- ii[which.min(cont[ii])]
      sc[ii] <- exp((log(sc[i])*cont[i] + log(sc[ii])*cont[ii])/(cont[i] + cont[ii]))
      sc <- sc[-i]
      tmp[,ii] <- tmp[,i] + tmp[,ii] ; tmp <- tmp[,-i,drop=FALSE]
      cont[ii] <- cont[i] + cont[ii] ; cont <- cont[-i]
    } else break
  }
  ##
  colnames(tmp) <- paste(sp," <",1L:ncol(tmp),">",sep="")
  ##
  Density[["fused"]][["density"]] %<>% cbind(tmp)
  Density[["fused"]][["size"]] %<>% c(sc)
}
rm(sp,sc,tmp,cont,i,ii)
names(Density[["fused"]][["size"]]) <-
  Density[["fused"]][["density"]] %>%
  colnames
##
## Density[["fused"]][["density"]] %>% {100*sum(as.numeric(.)==0)/(nrow(.)*ncol(.))}
cont <-
  Density[["fused"]][["density"]] %>%
  {100*colSums(.>0)/nrow(.)}
mincont <- 7.5
## 100*sum(cont>mincont)/length(cont)
Density[["fused"]][["size"]] %<>%
  {.[cont>mincont]}
Density[["fused"]][["density"]] %<>%
  {.[,cont>mincont]}
## Density[["fused"]][["density"]] %>% {100*sum(as.numeric(.)==0)/(nrow(.)*ncol(.))}
rm(cont,mincont)
##
## Density[["fused"]][["density"]] %>% {100*colSums(.>0)/nrow(.)}
tmp <-
  Density[["fused"]][["density"]] %>%
  colnames %>%
  strsplit(" <")
Density[["fused"]][["species"]] <-
  tmp %>% length %>% character
for(i in 1L:length(tmp))
  Density[["fused"]][["species"]][i] <- tmp[[i]][1L]
rm(tmp,i)
##
## Density[["fused"]][["species"]] %>% unique %>% length ## 48 species
## Density[["fused"]][["species"]] %>% length            ## 142 species/size classes
## Density[["fused"]][["species"]] %>% {tapply(.,.,length)} %>% as.factor %>% summary
##
## Density %>% str
##
### Fifth: Arrange trees and order data with respect to tree order
TreeHNet <- list()
##
TreeHNet[["raw"]] <-
  SupplSpeciesTree %>%
  drop.tip(
    tip=.$tip.label %>%
      match(Density[["fused"]][["species"]] %>% unique) %>%
      is.na %>%
      SupplSpeciesTree$tip.label[.]
  )
##
for(i in 1L:TreeHNet[["raw"]]$Nnode) {
  ## i=26L
  TreeHNet[["raw"]] %>%
    root(n=length(.$tip.label)+i) %>%
    plot
  if(is.null(locator(1))) {
    TreeHNet[["raw"]] %<>%
      root(n=length(.$tip.label)+i)
    break
  }
}
rm(i)
## TreeHNet[["raw"]] %>% plot
##
TreeHNet[["expanded"]] <- TreeHNet[["raw"]]
for (sp in unique(Density[["fused"]][["species"]])) {
  ## sp <- unique(Density[["fused"]][["species"]])[1L]
  n <- sum(Density[["fused"]][["species"]]==sp)
  if(n>1L) {
    str <-
      LETTERS[1L:n] %>%
      paste(":0", sep="") %>%
      paste(collapse=",") %>%
      paste("(", ., ");", sep="") %>%
      read.tree(text = .)
    str$tip.label <-
      Density[["fused"]] %>%
      {colnames(.$density)[.$species==sp]}
    TreeHNet[["expanded"]] %<>%
      bind.tree(str,where=which(.$tip.label==sp))
  } else {
    TreeHNet[["expanded"]]$tip.label[TreeHNet[["expanded"]]$tip.label==sp] <-
      Density[["fused"]] %>%
      {colnames(.$density)[.$species==sp]}
  }
}
rm(sp,n,str)
TreeHNet[["expanded"]]$node.label <- NULL  ## Drop the (useless) node labels.
## 
## TreeHNet[["expanded"]] %>% plot
TreeHNet[["expanded"]] %<>% write.tree %>% read.tree(text=.)
TreeHNet[["expanded"]]$tip.label %<>%
  strsplit("_") %>%
  lapply(function(x) paste(x, collapse=" ")) %>%
  unlist
## TreeHNet[["expanded"]] %>% plot
##
### Using Lethenteron appendix as an outgroup for rooting
## TreeHNet[["expanded"]]$edge[TreeHNet[["expanded"]]$edge[,2L]==140,2L]
##
save(TreeHNet,file="../../Data/HNetTrees.rda")
##
X11(width=5.0,height=7.0)
par(mar=c(2,2,2,2))
TreeHNet[["raw"]] %>% plot(cex=0.8)
dev.copy2eps(file="../../Image/Phylogenetic tree.eps")
dev.off()
##
ord <-
  TreeHNet[["expanded"]]$tip.label %>%
  match(colnames(Density[["fused"]]$density))
Density[["fused"]]$size %<>% {.[ord]}
Density[["fused"]]$species %<>% {.[ord]}
Density[["fused"]]$density %<>% {.[,ord]}
rm(ord)
##
## all(TreeHNet[["expanded"]]$tip.label==colnames(Density[["fused"]]$density))
## Density[["fused"]]$species %>% {all(substr(TreeHNet[["expanded"]]$tip.label,1,nchar(.))==.)}
## all(TreeHNet[["expanded"]]$tip.label==labels(Density[["fused"]]$size))
##
## Density[["fused"]]$density %>% as.numeric %>% {mean(.==0)}
## Density[["fused"]]$density %>% as.numeric %>% density(from=0,to=0.05,bw=0.005) %>% plot(xlim=c(0,0.05))
##
save(Density, file="../../Data/Density.rda")
##
disconnectFromDB(psqlCon)
rm(psqlCon)
##
