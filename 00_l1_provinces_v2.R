#Defining oceanic provinces

# 1. Loading Packages and Data ----

packages<-c('terra','sf')

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

crs_albers_brasil <- "+proj=aea +lat_0=-12 +lon_0=-54 +lat_1=-2 +lat_2=-22 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

bat_dhn<-rast('input_data/bat_dhn_v2.tif')
bat_gebco<-rast('input_data/gebco/gebco_2026_n10.0_s-38.0_w-55.0_e-24.0_geotiff.tif')
ss<-read_sf('input_data/study_area.shp')
ss<-st_transform(ss,crs_albers_brasil)

bat<-project(bat_gebco,bat_dhn)
plot(bat)
bat<-crop(bat,ss)
# bat<-subst(bat, from=0, to=NA)

bat<-clamp(bat, upper = 0)

plot(bat)




# 2. Calculate Slope (in degrees)-----------
# Essential to distinguish flat plains (basins/shelves) from steep structural zones (slopes/talus)
slope_deg <- terrain(bat, v = "slope", unit = "degrees")
plot(slope_deg)
## Focal filter---------

focal_window <- matrix(1, nrow = 11, ncol = 11)

slope_degf <- focal(slope_deg, 
                    w = focal_window,
                    fun = 'max', 
                    na.rm = TRUE)

# slope_degf <- focal(slope_deg, 
#                     w = focal_window,
#                     fun = function(x, ...) quantile(x, probs = 0.90, na.rm = TRUE))
                    


# slope_deg<-aggregate(slope_deg,fact=10,method='max')
plot(slope_degf)
slope_deg<-slope_degf
hist(slope_deg)
# bat<-aggregate(bat,fact=10,method='mean')


plot(bat)


# 3. Apply the Classification Logic (Sequential Map Algebra)-----
benthic_zones <- rast(bat)
values(benthic_zones) <- NA

## Class 4: Abyssal Plain / Ocean Basin ----
# Everything deeper than the basin threshold, regardless of minor slope variations
benthic_zones[is.na(benthic_zones) & 
                bat <= -4500] <- 4


## Class 1: Continental Shelf (Plataforma Continental)------
# Everything shallower than the shelf break
benthic_zones[bat >= -250 & 
                slope_deg<0.5] <- 1
benthic_zones[bat >= -90] <- 1


## Class 2: Continental Slope (Talude Continental)--------
# Intermediate depths with a steep gradient
benthic_zones[is.na(benthic_zones) & 
                bat >= -2100] <- 2

benthic_zones[is.na(benthic_zones) &
                slope_deg>=1.5 &
                bat >= -3000] <- 2

## Class 3: Continental Rise (Sopé Continental)---------
# The remaining unclassified areas (intermediate to deep waters with a gentle gradient)
benthic_zones[is.na(benthic_zones) &
                slope_deg>=0.2] <- 3


benthic_zones[is.na(benthic_zones) & 
                bat <= -3500] <- 4

plot(benthic_zones)

# 4. Adding seamounts ----
# Load the seamounts shapefile  generated in QGIS
seamounts_vec <- read_sf('input_data/seamounts_v2.shp')
seamounts_vec<-st_transform(seamounts_vec,crs_albers_brasil)
seamounts_rast <- rasterize(seamounts_vec, benthic_zones, field = 5, background = NA)
plot(seamounts_rast)

benthic_zones<-cover(seamounts_rast,benthic_zones)
plot(benthic_zones)


#5. Noize Reduction-----------

## Sieve: Filter small polygons--------
pixels_tresholds <- 5000
benthic_zones_m <- terra::sieve(benthic_zones, 
                                threshold = pixels_tresholds, 
                                directions = 8)
plot(benthic_zones_m)

## Focal filter---------

focal_window <- matrix(1, nrow = 21, ncol = 21)
benthic_zones_m <- focal(benthic_zones_m, w = focal_window, fun = 'modal', na.rm = TRUE)
plot(benthic_zones_m)


benthic_zones_m<-cover(seamounts_rast,benthic_zones_m)

#6. Color table -------------- 
color_table_geo <- data.frame(
  value = c(1, 2, 3, 4, 5),
  color = c(
    "#00FFFF",  # 1: Continental Shelf (Cyan - representing shallow water)
    "#FFA500",  # 2: Continental Slope (Orange - representing steep transition)
    "#FFD700",   # 3: Continental Rise (Gold - representing gentle sediment accumulation at the base of the slope)
    "#000050",  # 4: Oceanic Basin (Navy - representing deep abyssal areas)
    "#FF4500"  # 5: Seamounts (Red - highlighting isolated topographic peaks)
  )
)

# The coltab() function embeds these colors directly into the raster's metadata
coltab(benthic_zones_m) <- color_table_geo

# Create the new, simplified Raster Attribute Table (RAT)
habitat_table <- data.frame(
  ID = c(1, 2, 3, 4, 5),
  Province = c(
    "A. Continental Shelf",
    "B. Slope",
    "C. Continental Rise",
    "D. Oceanic Basin",
    "E. Seamount"
  )
)

levels(benthic_zones_m) <- habitat_table

benthic_zones_m<-mask(benthic_zones_m,ss)
plot(benthic_zones_m)

# 7. Export the final classified product
writeRaster(benthic_zones_m, "output_data/l1_benthic_provinces_v2.tif", overwrite = TRUE)

