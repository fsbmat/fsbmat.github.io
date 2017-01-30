rm(list=ls())
graphics.off()
library(ggmap)
library("rgdal")
library(ggplot2)
library(scales)
pb_poligonos_rgdal <- readOGR(dsn=path.expand("~/GitHub/webpage/fsbmat.github.io/blog_posts/28-01-2017/MG"), layer="31MUE250GC_SIR", verbose=FALSE, stringsAsFactors=FALSE);
pb_dados <- slot(object=pb_poligonos_rgdal, name="data");






lista_municipios <- c("FLORESTAL", "LAJINHA", "VIÇOSA");
indice_numerico <- which( pb_dados$NM_MUNICIP %in% lista_municipios );
print(indice_numerico);
indice_numerico[1]
mun_florestal <- pb_poligonos_rgdal[indice_numerico[1], ];
lat <- mun_florestal@polygons[[1]]@Polygons[[1]]@coords[,1]
lon <- mun_florestal@polygons[[1]]@Polygons[[1]]@coords[,2]
library(sp)
coordinates(mun_florestal)


mapImageData1 <- get_map(location = c(lon = -44.4321, lat = -19.8866),
                         color = "color",
                         source = "google",
                         maptype = "satellite",
                         zoom = 14)

ggmap(mapImageData1,
      extent = "device",
      ylab = "Latitude",
      xlab = "Longitude")


mapImageData4 <- get_map(location = c(lon = -44.4321, lat = -19.8866),
                         color = "color",
                         source = "google",
                         maptype = "hybrid",
                         zoom = 15)

ggmap(mapImageData4,
      extent = "panel",
      ylab = "Latitude",
      xlab = "Longitude")

# contour overlay
ggmap(get_map(maptype = "satellite"), extent = "device") +
  stat_density2d(aes(x = lon, y = lat, colour = class), data = chkpts, bins = 5)


# CARGA LIBRERIA Y DATOS
# -----------------------------------------
library(ggplot2)
library(ggmap)
geodata <- data.frame(
  lat = c(-19.6555,-19.6559,-19.7085,-19.5877,-19.6995,-19.6037,-19.7678,-19.7133, -19.7966),
  lon = c(-44.4321,-44.6167,-44.5859,-44.532,-44.3921,-44.3816,-44.3792,-44.3711, -44.276)
)


# MAPA 1
#-------------------------------------
qmplot(lon, lat, data = geodata, 
       colour = I("red"), 
       size   = I(6), 
       darken = .2,
       alpha  = 0.2,
       source = 'google')


# MAPA 2
#---------------------------------------
ciudad <- get_map("Florestal, Brasil", zoom=10)
p     <- ggmap(ciudad)
p + geom_point(data = geodata, 
               aes(x = lon, y = lat), 
               color = "red", 
               size  = 11, 
               alpha = 0.5)



# MAPA 3
# ----------------------------------------
ciudad <- get_map("Florestal, Brasil", zoom=15,maptype = "hybrid")
p      <- ggmap(ciudad)
p      <- p+stat_density2d(aes(x = lon, y = lat, fill=..level..), 
                           data=geodata,geom="polygon", alpha=0.2)
p + scale_fill_gradient(low = "yellow", high = "red")







