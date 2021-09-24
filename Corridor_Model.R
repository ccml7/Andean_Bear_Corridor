source("ENM.R")
library(classInt)
library(dichromat)
library(leastcostpath)
library(gdistance)

?TransitionLayer

# values_t <- values(enviromental_data[[1]])
# p <- classIntervals(values_t, n = 10, style = "jenks")

# image(enviromental_data[[1]], breaks = p$brks, col = colorRampPalette(c("red", "yellow"))(10))

# p$brks
# p$var

# enviromental_data[[1]]

origin_points <- shapefile("./analysis_layers/shapes/Origin_Destination_Corridor.shp") #nolint

enm_projection <- raster("./outputs/Niche_Model_Prediction.tif")

# inv_projection <- (enm_projection - 1) * -1
# inv_projection <- inv_projection / 10
# plot(inv_projection)

transition_bear <- transition(enm_projection,
                             transitionFunction = mean,
                             directions = 16)

plot(raster(transition_bear))

corridor <- create_cost_corridor(cost_surface = transition_bear,
                    origin = origin_points[1, ],
                    destination = origin_points[2, ],
                    rescale = TRUE
)

lcp <- create_lcp(cost_surface = transition_bear,
                    origin = origin_points[1, ],
                    destination = origin_points[2, ],
                    directional =  FALSE
)

plot(enm_projection)
plot(lcp, add = T)
points(origin_points)

lcp_planar <- spTransform(lcp, origen_bogota)
corridor_buffer <- buffer(lcp_planar, width = 2500)

shapefile(lcp, "./outputs/LCP_wgs84.shp")
shapefile(lcp_planar, "./outputs/LCP_planar.shp")
shapefile(corridor_buffer, "./outputs/corridor_2500m.shp")
