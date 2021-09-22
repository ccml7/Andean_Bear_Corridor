library(dismo)
library(raster)
library(rgeos)

## Study Area
niche_m <- shapefile("./analysis_layers/shapes/buffer_oc_wgs84.shp")

## Occurrence points
ocurrences <- read.table("./data/ocurrencias_oso.txt",
    header = T, sep = "\t", dec = ",")

## Enviromental data
files <- list.files("./analysis_layers/raster/predictors",
            pattern = ".tif")

enviromental_data <- stack(paste("./analysis_layers/raster/predictors/",
                            files, sep = ""))


## Bias file
bias_layer <- raster("./analysis_layers/raster/bias_raster.tif")

## Background points
set.seed(1)
SpatialPoints()
background_points <- spsample(
                        x = niche_m,
                        n = 10000,
                        type = "random")

## Set the training and testing data

# Select the 50% of ocurrences for training
selected <- sample(seq(1, nrow(ocurrences)), nrow(ocurrences) * 0.5)

occ_training <- ocurrences[selected, ]
occ_test <- ocurrences[-selected, ]

## Extract climatic data for occurence and background points
env_training <- extract(enviromental_data, occ_training[, 2:3])
env_test <- extract(enviromental_data, occ_test[, 2:3])

env_background <- extract(enviromental_data, background_points)

presence_ausence <- c(rep(1, nrow(occ_training)),
                    rep(0, nrow(coordinates(background_points))))

unified_points <- data.frame(rbind(env_training, env_background))

## Training Maxent Model with tabular data
bear_maxent <- maxent(
                    x = enviromental_data,
                    p = occ_training[, c(3, 2)],
                    args = c("responsecurves")
)

env_test

dim(unified_points)
length(presence_ausence)

presence_ausence

install.packages("rJava")
install.packages("devtools")

library(rJava)
library(devtools)

devtools::install_version("rJava", version = "0.9-9", repos = "http://cran.us.r-project.org")
