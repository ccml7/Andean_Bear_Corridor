library(raster)

setwd("/home/camilo/Documentos/Mariana/Corredor/Tremarctos_Corridor")

oc_gbif <- read.csv("./data/ocurrenciasGBIF.csv",
     header = T, sep = ";", dec = ".")

oc_cortolima <- read.csv("./data/registrososo_cortolima.csv",
     header = T, sep = ";", dec = ",")

head(oc_gbif)
