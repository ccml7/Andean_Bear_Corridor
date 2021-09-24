source("./Data_Preparation.R")

## Upload Worldclim Tiles and Nubosity Layers
setwd("/home/camilo/Documentos/Capas/")

wgs_84 <- "+proj=longlat +datum=WGS84 +no_defs"
origen_bogota <- "+proj=tmerc +lat_0=4.596200416666666 +lon_0=-74.07750791666666 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" # nolint

worldclim <- getData("worldclim", var = "bio", res = 0.5,
     lon = -75.8, lat = 3.7)

earth_env <- raster("MODCF_meanannual.tif")

setwd("/home/camilo/Documentos/Mariana/Corredor/Tremarctos_Corridor")

#############################################################################
## IDEAM Forest No-Forest Layer
forest <- raster("./layers/bosque_tol.tif")
forest_wgs <- projectRaster(forest, worldclim)
forest_cut <- mask(crop(forest_wgs, oc_def_buffer_wgs), oc_def_buffer_wgs)

worldclim_cut <- mask(crop(worldclim, oc_def_buffer_wgs), oc_def_buffer_wgs)
earth_env_cut <- mask(crop(earth_env, oc_def_buffer_wgs), oc_def_buffer_wgs)

dummy <- earth_env_cut
dummy[values(dummy) > 0] <- 0

for (i in seq(1, length(names(worldclim_cut)))) {
    writeRaster(worldclim_cut[[i]],
        paste("./analysis_layers/raster/predictors/",
        names(worldclim_cut)[i],
        "_NicheM.tif", sep = ""),
        overwrite = T)
}

writeRaster(earth_env_cut,
            "./analysis_layers/raster/predictors/Mean_Nubosity_NicheM.tif",
            overwrite = TRUE)

writeRaster(forest_cut, "./analysis_layers/raster/predictors/Forest.tif",
            overwrite = TRUE)


writeRaster(dummy,
            "./analysis_layers/raster/dummy.tif",
            overwrite = TRUE)

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

plot(tolima_shape_wgs84, axes = TRUE)
plot(oc_def_buffer_wgs, add = T)
plot(gbif_bias_sp, add = T, col = "blue")

inside_m <- over(gbif_bias_sp, oc_def_buffer_wgs)
inside_points <- which(inside_m == 1)
gbif_bias_sp <- gbif_bias_sp[inside_points, ]


plot(tolima_shape_wgs84, axes = TRUE)
plot(oc_def_buffer_wgs, add = T)
plot(gbif_bias_sp, add = T, col = "blue")

gbif_bias_sp[, 2]

bias_brick <- rasterize(gbif_bias_sp,
         dummy,
         1
)

bias_raster <- bias_brick[[1]]

## Add roads to Bias file
# roads <- shapefile("./layers/red_vial_dptal.shp")
# roads <- spTransform(roads, wgs_84)
# roads <- crop(roads, oc_def_buffer_wgs)

# roads_raster <- rasterize(roads, 
#          dummy,
#          1
# )

# bias_raster <- bias_raster + roads_raster
bias_raster[values(bias_raster) > 1] <- 1

presences <- which(values(bias_raster) == 1)
pres_locs <- coordinates(bias_raster)[presences, ]

dens <- kde2d(pres_locs[, 1], pres_locs[, 2],
             n = c(nrow(bias_raster), ncol(bias_raster)))

bias_raster <- raster(dens)
plot(bias_raster)

max_count_bias <- cellStats(bias_raster,
                            stat = "max")

values(bias_raster) <- values(bias_raster) / max_count_bias

plot(bias_raster)

writeRaster(bias_raster, 
            "./analysis_layers/raster/bias_raster.tif",
            overwrite = TRUE
)
