library(raster)
library(dismo)
library(rJava)
library(rgeos)
library(MASS)

setwd("/home/camilo/Documentos/Mariana/Corredor/Tremarctos_Corridor")

## Coordinate systems
origen_bogota <- "+proj=tmerc +lat_0=4.596200416666666 +lon_0=-74.07750791666666 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" # nolint
wgs_84 <- "+proj=longlat +datum=WGS84 +no_defs"

tolima_shape_planar <- shapefile("./analysis_layers/shapes/tolima_planar.shp")
tolima_shape_wgs84 <- shapefile("./analysis_layers/shapes/tolima_wgs84.shp")

oc_gbif <- read.csv("./data/gbif.csv",
     header = T, sep = ";", dec = ",")

oc_cortolima <- read.csv("./data/registrososo_cortolima.csv",
     header = T, sep = ";", dec = ",")

names(oc_gbif) <- c("Latitude", "Longitude")
names(oc_cortolima) <- c("Latitude", "Longitude")

oc_def <- rbind(oc_cortolima, oc_gbif)
id_punto <- seq(1, dim(oc_def)[1], 1)

oc_def <- data.frame(id_punto, oc_def)

## From oriental cordillera
oc_cordillera_central <- which(oc_def$Longitude < -74.7)
oc_def <- oc_def[oc_cordillera_central, ]

## From less than 2000 m.a.s.l
oc_tierras_bajas <- which(oc_def$id_punto == 49)
oc_def <- oc_def[-oc_tierras_bajas, ]

oc_def_sp <- oc_def
coordinates(oc_def_sp) <- ~Longitude + Latitude
crs(oc_def_sp) <- wgs_84

## Filter anomalous points
## Outside Tolima
inside_tolima <- over(oc_def_sp, tolima_shape_wgs84)$COD_DEPART
outside <- which(is.na(inside_tolima) == TRUE)

oc_def <- oc_def[-outside, ]

## Transform again ocurrences to spatial points
oc_def_sp <- oc_def
coordinates(oc_def_sp) <- ~Longitude + Latitude
crs(oc_def_sp) <- wgs_84

plot(tolima_shape_wgs84, axes = T)
points(oc_def_sp)

shapefile(oc_def_sp, "./analysis_layers/shapes/puntos_oso.shp",
          overwrite = TRUE)

write.table(oc_def, "./data/ocurrencias_oso.txt",
     sep = "\t", row.names = FALSE, dec = ",")

## Reproject ocurrences to planar coordinates
oc_def_projected <- spTransform(oc_def_sp, origen_bogota)
shapefile(oc_def_projected, "./analysis_layers/shapes/puntos_oso.shp",
          overwrite = TRUE)

oc_def_buffer <- buffer(oc_def_projected, width = 20000, dissolve = T)
plot(tolima_shape_planar)
plot(oc_def_buffer, add = T)

### Crop buffer by Tolima limits
crop_buffer_oc <- crop(oc_def_buffer, tolima_shape_planar)
plot(crop_buffer_oc, col = "red", add = T)
points(oc_def_projected, col = "black", pch = 16)

shapefile(crop_buffer_oc,
      "./analysis_layers/shapes/buffer_oc_planar.shp",
       overwrite = T)

oc_def_buffer_wgs <- spTransform(crop_buffer_oc, wgs_84)

shapefile(oc_def_buffer_wgs,
          "./analysis_layers/shapes/buffer_oc_wgs84.shp",
          overwrite = T)
