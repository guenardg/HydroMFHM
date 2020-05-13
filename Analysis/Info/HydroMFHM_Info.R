## **************************************************************************
##
##    (c) 2020 Guillaume Guénard
##        Université de Montréal, Montreal, Quebec, Canada
##
##    **Info R script**
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
Info <- list()
##
### Get information about the rivers
Info[["river"]] <-
  paste("SELECT river, province, dam, nseg",
        "FROM number_of_segments_electric",
        "WHERE river != 'Beaurivage' ORDER BY river, province;") %>%
  dbGetQuery(psqlCon,.)
rownames(Info[["river"]]) <- Info[["river"]][,"river"]
Info[["river"]] <- Info[["river"]][,-1L]
##
### Get information about the sites
Info[["site"]] <-
  paste("SELECT *",
        "FROM all_sites_visual",
        "WHERE river != 'Beaurivage'") %>%
  dbGetQuery(psqlCon,.)
rownames(Info[["site"]]) <- Info[["site"]][,"site"]
Info[["site"]] <- Info[["site"]][,-c(3L,4L,6L,8L)]
##
### Get information about the species that have been observed
Info[["species"]] <-
  paste("SELECT river, value AS species, COUNT(*) AS n",
        "FROM species_occurrence_visual",
        "GROUP BY river, value ORDER BY river, value") %>%
  dbGetQuery(psqlCon,.)
##
### Get the total density file (to extract site information)
Info[["density"]] <-
  "../../Data/Dtb_Guillaume_Guenard.csv" %>%
  read.csv(stringsAsFactors=FALSE)
##
### Getting the site variables right
x <- c("Natural","RunOfTheRiver","Peaking","Storage")
y <- match(rownames(Info[["site"]]),Info[["density"]][,"Site"])
Info[["site"]][,"flow_regime"] <-
  as.factor(x[Info[["density"]][y,"FlowRegimeM"]])
## Following a long mail thread, that river's Flow regime had been changed
Info[["site"]][Info[["site"]][,"river"]=="Waterton","flow_regime"] <- "Storage"
Info[["site"]][,"depth"] <- Info[["density"]][y,"Depth"]
Info[["site"]][,"velocity"] <- Info[["density"]][y,"Velocity"]
Info[["site"]][,"d50"] <- Info[["density"]][y,"D50"]
Info[["site"]][,"CovOK"] <-
  !(is.na(Info[["site"]][,"depth"]) |
    is.na(Info[["site"]][,"velocity"]) |
    is.na(Info[["site"]][,"d50"]))
Info[["density"]] <- NULL
rm(x,y)
##
### River's mid sampling point to be used as marker for large-scale maps
tmp <- 0.5*(Info[["site"]][,"lat_from"] + Info[["site"]][,"lat_to"])
tmp %<>% tapply(Info[["site"]][,"river"],mean)
Info[["river"]][names(tmp),"lat"] <- tmp
tmp <- 0.5*(Info[["site"]][,"lon_from"] + Info[["site"]][,"lon_to"])
tmp %<>% tapply(Info[["site"]][,"river"],mean)
Info[["river"]][names(tmp),"lon"] <- tmp
rm(tmp)
##
### Set Excluded data aside
Info[["excluded"]] <- Info[["site"]][!Info[["site"]][,"CovOK"],]
Info[["site"]] <- Info[["site"]][Info[["site"]][,"CovOK"],]
Info[["site"]][,"flow_regime"] %<>% as.character %>% as.factor
tmp <- Info[["site"]][,"flow_regime"] %>%
  as.character %>%
  tapply(Info[["site"]][,"river"],head,n=1)
Info[["river"]][names(tmp),"flow_regime"] <- tmp %>% as.factor
Info[["site"]] %<>% {.[,colnames(.)!="CovOK"]}
Info[["excluded"]] %<>% {.[,colnames(.)!="CovOK"]}
rm(tmp)
##
### Organize sites by rivers west -> east
Info[["site"]] <-
  Info[["river"]][,"lon"] %>%
  sort(index.return=TRUE) %>%
  {rownames(Info[["river"]])[.$ix]} %>%
  sapply(function(x,y) which(y==x), y=Info[["site"]][,"river"]) %>%
  unlist %>%
  Info[["site"]][.,]
##
### Calculate the number of sites per rivers
nsite <-
  Info[["site"]][,"river"] %>%
  tapply(.,.,length)
## nsite[which.min(nsite)]
## nsite[nsite==max(nsite)]
Info[["river"]][,"N_OK"] <- NA
Info[["river"]][labels(nsite)[[1L]],"N_OK"] <- nsite
rm(nsite)
##
### Organize rivers west -> east
Info[["river"]] %<>%
  {.[sort(.[,"lon"],index.return=TRUE)$ix,]}
##
EnvMod <- list()
##
### Here, depth, velocity, and substrate were aggregared by river because the
### analysis is to be performed on the basis 
##
EnvMod[["depth"]] <- lm(log(depth)~river-1, data=Info[["site"]])
## EnvMod[["depth"]] %>% residuals %>% shapiro.test
## EnvMod[["depth"]] %>% residuals %>% density %>% plot
## EnvMod[["depth"]] %>% summary
Info[["river"]]["log_depth"] <-
  EnvMod[["depth"]] %>%
  coefficients %>%
  {.[paste("river", Info[["river"]] %>% rownames, sep="")]} %>%
  as.numeric
##
EnvMod[["velocity"]] <- lm(log1p(velocity)~river-1, data=Info[["site"]])
## EnvMod[["velocity"]] %>% residuals %>% shapiro.test
## EnvMod[["velocity"]] %>% residuals %>% density %>% plot
## EnvMod[["velocity"]] %>% summary
Info[["river"]]["log1p_velocity"] <-
  EnvMod[["velocity"]] %>%
  coefficients %>%
  {.[paste("river", Info[["river"]] %>% rownames, sep="")]} %>%
  as.numeric
##
EnvMod[["d50"]] <- lm(log1p(d50)~river-1, data=Info[["site"]])
## EnvMod[["d50"]] %>% residuals %>% shapiro.test
## EnvMod[["d50"]] %>% residuals %>% density %>% plot
## EnvMod[["d50"]] %>% summary
Info[["river"]]["log1p_d50"] <-
  EnvMod[["d50"]] %>%
  coefficients %>%
  {.[paste("river", Info[["river"]] %>% rownames, sep="")]} %>%
  as.numeric
##
### Other descriptors by rivers
RD <- "../../Data/Descripteurs Riviere.csv" %>%
  read.csv(header = TRUE, sep=";", dec=",")
rownames(RD) <- RD[,1L]
RD <- RD[,-1L,drop=FALSE]
Info[["river"]]["degree-day"] <- RD[rownames(Info[["river"]]),"DD"]
Info[["river"]]["log_tot_phos"] <- log10(RD[rownames(Info[["river"]]),"TP"])
Info[["river"]]["macro_covr"] <- RD[rownames(Info[["river"]]),"MC"]
rm(RD)
##
Info[["density"]] <-
  "../../Data/densites_SPCL_sites.csv" %>%
  read.csv(as.is=TRUE, dec=",",sep=";")
rownames(Info[["density"]]) <- Info[["density"]][,"Site"]
Info[["density"]] %<>%
  {.[rownames(Info[["site"]]),-(1L:9L)]}
## Info[["density"]] %>% head
##
## Info[["density"]] %>% {100*sum(colSums(.)<=sqrt(.Machine$double.eps))/ncol(.)}
##
save(Info, file="../../Data/Info.rda")
##
disconnectFromDB(psqlCon)
rm(psqlCon)
##
