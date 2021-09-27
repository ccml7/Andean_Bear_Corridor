library(dismo)
library(raster)
library(rgeos)
library(rJava)


## Study Area
niche_m <- shapefile("./analysis_layers/shapes/buffer_oc_wgs84.shp")

## Occurrence points
ocurrences <- read.table("./data/ocurrencias_oso.txt",
    header = T, sep = "\t", dec = ",")

ocurrences_sp <- shapefile("./analysis_layers/shapes/puntos_oso.shp")

## Enviromental data
files <- list.files("./analysis_layers/raster/predictors",
            pattern = ".tif")

enviromental_data <- stack(paste("./analysis_layers/raster/predictors/",
                            files, sep = ""))


## Bias file
bias_layer <- raster("./analysis_layers/raster/bias_raster.tif")

## Background points
# set.seed(1)
# SpatialPoints()
# background_points <- spsample(
#                         x = niche_m,
#                         n = 10000,
#                         type = "random")

# Select the 50% of ocurrences for training
selected <- sample(seq(1, nrow(ocurrences)), nrow(ocurrences) * 0.5)

env_bg <- extract(enviromental_data, gbif_bias_sp)
env_oc <- extract(enviromental_data, ocurrences_sp)

env_oc_training <- env_oc[selected, ]
env_oc_test <- env_oc[-selected, ]

## Set the training and testing data

occ_training <- ocurrences[selected, ]
occ_test <- ocurrences[-selected, ]

## Training Maxent Model with tabular data
bear_maxent <- maxent(
                    x = enviromental_data,
                    p = occ_training[, c(3, 2)],
                    a = gbif_bias_sp,
                    args = c("responsecurves"),
                    path = paste0("./outputs/maxent_output")
)

bear_maxent@results
pred1 <- predict(bear_maxent, enviromental_data)
plot(pred1)

writeRaster(pred1,
            "./outputs/Niche_Model_Prediction.tif",
            overwrite = T)

## Model Evaluation
mod_eval_train <- dismo::evaluate(p = env_oc_training,
                                  a = env_bg,
                                  model = bear_maxent)
print(mod_eval_train)

mod_eval_test <- dismo::evaluate(p = env_oc_test,
                                  a = env_bg,
                                  model = bear_maxent)
print(mod_eval_test)


## training AUC may be higher than testing AUC

thd1 <- threshold(mod_eval_train, "no_omission")
thd2 <- threshold(mod_eval_train, "spec_sens")
