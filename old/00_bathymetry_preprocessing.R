#Preprocessing DHN bathymetric model

# 1. Loading Packages and Data ----

packages<-c('terra','sf')

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})



# 2. Creating Projection ----
# Using SIRGAS 2000 / Brazil Albers Equal Area (EPSG:10857) [Image of Albers Equal Area Conic projection]
# Crucial for Natural Neighbors interpolation and area calculations, 
# as it preserves area proportions across the continental scale.
crs_albers_brasil <- "+proj=aea +lat_0=-12 +lon_0=-54 +lat_1=-2 +lat_2=-22 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"


## Exclusive Economic Zone (EEZ) Mask ----
eez <- read_sf('input_data/study_area.shp')
plot(eez,add=T)
eez<-st_buffer(eez,200000)
eez <- st_transform(eez, crs_albers_brasil)

# Creating a base raster grid with 1km (1000m) resolution
mask <- rast(eez, res = 1000)
mask <- rasterize(eez, mask, field = 'id')


# 3. Converting DHN bathymetry----------
bat<-rast('input_data/DTM_BR_LEPLAC_MAR21_1000m_surfer_V7_0/DTM_BR_LEPLAC_MAR21_1000m_surfer_V7.GRD')
plot(bat)
crs(bat)
set.crs(bat,"EPSG:3395")
res(bat)
bat<-subst(bat, from=NA, to=0)

bat<-project(bat, mask,method='cubic')
plot(bat)
bat<-mask(bat,mask)
plot(bat)
writeRaster(bat,'input_data/bat_dhn_v2_mask.tif',overwrite=T)

