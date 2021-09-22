library(dismo)
library(raster)

## Occurrence points
ocurrence <- read.table("./data/ocurrencias_oso.txt",
    header = T, sep = "\t", dec = ",")

## Enviromental data 