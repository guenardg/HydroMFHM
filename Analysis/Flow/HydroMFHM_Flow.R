## **************************************************************************
##
##    (c) 2020 Guillaume Guénard
##        Université de Montréal, Montreal, Quebec, Canada
##
##    **Flow metrics R script**
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
##
library(DBI)
library(magrittr)
library(yaml)
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
load("../../Data/Key Flow R objects 2016-02-17.RData")
##
### Raw data are not needed here.
rm(flow.D.9709,flow.D.9713,flow.D.bku,flow.H.9709,flow.H.9713,flow.H.bku)
##
Flow_metrics <- list()
##
Flow_metrics[["description"]] <-
  "../../Data/Flow_metrics.csv" %>%
  read.csv(as.is=TRUE,sep=";")
rownames(Flow_metrics[["description"]]) <-
  Flow_metrics[["description"]][["Code"]]
Flow_metrics[["description"]] %<>% {.[-1L]}
## Flow_metrics[["description"]] %>% head
##
## summary(flow.indices.HN)
## summary(flow.indices.HNn)
## rownames(flow.indices.HN)
## rownames(Info[["river"]])
nms <- Info[["river"]] %>% nrow %>% character
names(nms) <- rownames(Info[["river"]])
nms[] <- c("ab.kananaskis.pocaterra","ab.elbow","ab.castle","ab.waterton.reg",
           "on.magpie","on.batchawana","on.goulais","on.aubinadong",
           "on.mississagi","qc.picanoc","qc.kiamika","qc.noire","qc.nicolet",
           "qc.ste_anne","qc.coaticook","qc.eaton","qc.st_francois",
           "qc.becancour.inverness","qc.aux_saumons","qc.etchemin","qc.du_sud",
           "qc.st_jean","qc.petit_saguenay","qc.ouelle","qc.du_loup_sjdk",
           "nb.gulquac","nb.dee","nb.serpentine")
## flow.indices.HN[nms,]
## flow.indices.HNn[nms,]
## ncol(flow.indices.HN)
## ncol(flow.indices.HNn)
Flow_metrics[["values"]] <-
  matrix(NA,nrow(flow.indices.HN),0L,
         dimnames=list(rownames(flow.indices.HN),NULL)) %>%
  data.frame
##
cn <- list(colnames(flow.indices.HN),colnames(flow.indices.HNn))
for(i in 1L:ncol(flow.indices.HN)) {
  ## i <- 1L
  nan <- !is.nan(flow.indices.HN[,i])
  na <- !is.na(flow.indices.HN[,i])
  Flow_metrics[["values"]][[cn[[1L]][i]]] <- flow.indices.HN[,i]
  if(any(flow.indices.HN[nan&na,i] != flow.indices.HNn[nan&na,i]))
    Flow_metrics[["values"]][[cn[[2L]][i]]] <- flow.indices.HNn[,i]
}
rm(cn,i,nan,na)
##
macnaughton_et_al <- c("MA3","nML6","FH1","DL12","DH6","TA2","TH2","RA7",
                       "nRA1","RL2","MA60")
Flow_metrics[["values"]] %<>% {.[nms,macnaughton_et_al]}
rownames(Flow_metrics[["values"]]) <- labels(nms)
tmp <- macnaughton_et_al %>%
  sapply(function(x) if(substr(x,1L,1L)=="n") substr(x,2L,nchar(x)) else x)
Flow_metrics[["description"]] %<>% {.[tmp,]}
rownames(Flow_metrics[["description"]]) <- macnaughton_et_al
rm(macnaughton_et_al,tmp,nms)
##
save(Flow_metrics, file="../../Data/Flow_metrics.rda")
##
disconnectFromDB(psqlCon)
rm(psqlCon)
##
