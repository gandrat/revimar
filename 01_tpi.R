
# 1. Loading Packages and Data ----

packages<-c('terra','sf','spatialEco')

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

crs_albers_brasil <- "+proj=aea +lat_0=-12 +lon_0=-54 +lat_1=-2 +lat_2=-22 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"


bat<-rast('input_data/gebco/gebco_2026_n10.0_s-38.0_w-55.0_e-24.0_geotiff.tif')
bat<-project(bat,crs_albers_brasil,method='bilinear')
ss<-read_sf('input_data/study_area.shp')
ss<-st_transform(ss,crs_albers_brasil)

res(bat)
plot(bat)
bat<-crop(bat,ss)

bat<-clamp(bat, upper = 0)

plot(bat)


# 3. Calculate Broad-scale and Fine-scale BPI----------
# The 'normalize = TRUE' argument standardizes the output (mean = 0, sd = 1).
# The 'scale' argument defines the moving window size in pixels. 


# TPI all scales
bpi_9 <- tpi(bat, scale = 9, win = "rectangle", normalize = TRUE)
bpi_21 <- tpi(bat, scale = 21, win = "rectangle", normalize = TRUE)
bpi_51 <- tpi(bat, scale = 51, win = "rectangle", normalize = TRUE)
bpi_501<-tpi(bat, scale = 501, win = "rectangle", normalize = TRUE)
bpi_101 <- tpi(bat, scale = 101, win = "rectangle", normalize = TRUE)

bpi<-c(bpi_9,bpi_21, bpi_51, bpi_101,bpi_501)
names(bpi)<-c('s9','s21','s51','s101','s501')

plot(bpi)

writeRaster(bpi,'output_data/bpi_v2.tif',overwrite=T)


#4. Reclassifying------
bpi<-rast('output_data/bpi_v2.tif')
res(bpi)

##Low TPI (Canyons)------

can<-bpi$s501<=-0.5

plot(can)

focal_window <- matrix(1, nrow = 21, ncol = 21)
canf <- focal(can, w = focal_window, fun = 'modal', na.rm = TRUE)
plot(canf)

can <- terra::sieve(canf,
                    threshold = 1000, 
                    directions = 8)
plot(can)

can<-mask(can,ss)
writeRaster(can,'output_data/l1_tpi_low.tif',overwrite=T)

##Hight TPI (Mountains)
mount<-bpi$s501>=.5

plot(mount)

focal_window <- matrix(1, nrow = 21, ncol = 21)
mountf <- focal(mount, w = focal_window, fun = 'modal', na.rm = TRUE)
plot(mountf)

mount <- terra::sieve(mountf,
                    threshold = 1000, 
                    directions = 8)
plot(mount)

mount<-mask(mount,ss)
writeRaster(mount,'output_data/l1_tpi_high.tif')

