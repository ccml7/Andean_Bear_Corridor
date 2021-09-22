library(raster)
setwd("/home/camilo/Documentos/Capas/")

wgs_84 <- "+proj=longlat +datum=WGS84 +no_defs"
origen_bogota <- "+proj=tmerc +lat_0=4.596200416666666 +lon_0=-74.07750791666666 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" # nolint


worldclim <- getData("worldclim", var = "bio", res = 0.5,
     lon = -75.8, lat = 3.7)
earth_env <- raster("MODCF_meanannual.tif")


setwd("/home/camilo/Documentos/Mariana/Corredor/Tremarctos_Corridor")

tolima_shape <- shapefile("./analysis_layers/shapes/tolima_wgs84.shp")
niche_m <- shapefile("./analysis_layers/shapes/buffer_oc_wgs84.shp")

worldclim_cut <- mask(crop(worldclim, niche_m), niche_m)
earth_env_cut <- mask(crop(earth_env, niche_m), niche_m)

for (i in seq(1, length(names(worldclim_cut)))) {
    writeRaster(worldclim_cut[[i]],
        paste("./analysis_layers/raster/predictors/",
        names(worldclim_cut)[i],
        "_NicheM.tif", sep = ""),
        overwrite = T)
}

writeRaster(earth_env_cut,
            "./analysis_layers/raster/predictors/Mean_Nubosity_NicheM.tif")

## Bias Layer
gbif_bias <- read.csv("./data/Gbif_Bias_Data.csv", header = T, sep = "\t")
names(gbif_bias)

# gbif_bias <- data.frame(ID = seq(1, dim(gbif_bias)[1]),
#                         Order = gbif_bias$order,
#                         Family = gbif_bias$family,
#                         Species = gbif_bias$verbatimScientificName,
#                         Longitude = gbif_bias$decimalLongitude,
#                         Latitude = gbif_bias$decimalLatitude)
# write.table(gbif_bias, "./data/Gbif_Bias_Data.csv", sep = "\t", row.names = F) # nolint

gbif_bias_sp <- gbif_bias
coordinates(gbif_bias_sp) <- ~Longitude + Latitude
crs(gbif_bias_sp) <- wgs_84

plot(tolima_shape, axes = TRUE)
plot(niche_m, add = T)
plot(gbif_bias_sp, add = T, col = "blue")

inside_m <- over(gbif_bias_sp, niche_m)
inside_points <- which(inside_m == 1)
gbif_bias_sp <- gbif_bias_sp[inside_points, ]


plot(tolima_shape, axes = TRUE)
plot(niche_m, add = T)
plot(gbif_bias_sp, add = T, col = "blue")


dummy <- earth_env_cut
dummy[values(dummy) > 0] <- 0

bias_brick <- rasterize(gbif_bias_sp,library(raster)
setwd("/home/camilo/Documentos/Capas/")

wgs_84 <- "+proj=longlat +datum=WGS84 +no_defs"
origen_bogota <- "+proj=tmerc +lat_0=4.596200416666666 +lon_0=-74.07750791666666 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" # nolint


worldclim <- getData("worldclim", var = "bio", res = 0.5,
     lon = -75.8, lat = 3.7)
earth_env <- raster("MODCF_meanannual.tif")


setwd("/home/camilo/Documentos/Mariana/Corredor/Tremarctos_Corridor")

tolima_shape <- shapefile("./analysis_layers/shapes/tolima_wgs84.shp")
niche_m <- shapefile("./analysis_layers/shapes/buffer_oc_wgs84.shp")

worldclim_cut <- mask(crop(worldclim, niche_m), niche_m)
earth_env_cut <- mask(crop(earth_env, niche_m), niche_m)

for (i in seq(1, length(names(worldclim_cut)))) {
    writeRaster(worldclim_cut[[i]],
        paste("./analysis_layers/raster/predictors/",
        names(worldclim_cut)[i],
        "_NicheM.tif", sep = ""),
        overwrite = T)
}

writeRaster(earth_env_cut,
            "./analysis_layers/raster/predictors/Mean_Nubosity_NicheM.tif")

## Bias Layer
gbif_bias <- read.csv("./data/Gbif_Bias_Data.csv", header = T, sep = "\t")
names(gbif_bias)

# gbif_bias <- data.frame(ID = seq(1, dim(gbif_bias)[1]),
#                         Order = gbif_bias$order,
#                         Family = gbif_bias$family,
#                         Species = gbif_bias$verbatimScientificName,
#                         Longitude = gbif_bias$decimalLongitude,
#                         Latitude = gbif_bias$decimalLatitude)
# write.table(gbif_bias, "./data/Gbif_Bias_Data.csv", sep = "\t", row.names = F) # nolint

gbif_bias_sp <- gbif_bias
coordinates(gbif_bias_sp) <- ~Longitude + Latitude
crs(gbif_bias_sp) <- wgs_84

plot(tolima_shape, axes = TRUE)
plot(niche_m, add = T)
plot(gbif_bias_sp, add = T, col = "blue")

inside_m <- over(gbif_bias_sp, niche_m)
inside_points <- which(inside_m == 1)
gbif_bias_sp <- gbif_bias_sp[inside_points, ]


plot(tolima_shape, axes = TRUE)
plot(niche_m, add = T)
plot(gbif_bias_sp, add = T, col = "blue")


dummy <- earth_env_cut
dummy[values(dummy) > 0] <- 0

bias_brick <- rasterize(gbif_bias_sp,
         dummy,
         fun = function(x, ...) {
            length(x)
         }
)

bias_count_raster <- bias_brick[[1]]
max(bias_count_raster)

max_count_bias <- cellStats(bias_count_raster,
                            stat = "max")


bias_count_raster <- bias_count_raster / max_count_bias
plot(bias_count_raster)

writeRaster(bias_count_raster, 
            "./analysis_layers/raster/bias_raster.tif",
            overwrite = TRUE
)

         dummy,
         fun = function(x, ...) {
            length(x)
         }
)

bias_count_raster <- bias_brick[[1]]
max(bias_count_raster)

max_count_bias <- cellStats(bias_count_raster,
                            stat = "max")


bias_count_raster <- bias_count_raster / max_count_bias
plot(bias_count_raster)

writeRaster(bias_count_raster, 
            "./analysis_layers/raster/bias_raster.tif",
            overwrite = TRUE
)
