#Defining oceanic provinces

# 1. Loading Packages and Data ----

packages<-c('terra','sf','dplyr')

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

# 3. Creating Spatial Masks ----

## 3.1 Exclusive Economic Zone (EEZ) Mask ----
eez <- read_sf('input_data/study_area.shp')
eez <- st_transform(eez, crs_albers_brasil)

# Creating a base raster grid with 1km (1000m) resolution
mask <- rast(eez, res = 1000)
mask <- rasterize(eez, mask, field = 'id')

# Ensuring the mask strictly adheres to the Albers projection and resolution
mask <- project(mask, crs_albers_brasil, res = 1000)

# writeRaster(mask, 'input_data/eez_albers.tif', overwrite = TRUE)

## 3.2 Bathymetric Model Mask ----
bat <- rast('input_data/batimetria_dhn.tif')

# Aligning bathymetry to the EEZ grid
bat <- project(bat, mask,method='bilinear')

bat<-subst(bat, from=NA, to=0)
bat <- mask(bat, mask)
res(bat)
plot(bat, main = "Bathymetry")

writeRaster(bat, 'input_data/batimetria_dhn_albers.tif', overwrite = TRUE)


mask_bat<-bat>=-1000
plot(mask_bat)
mask_bat<-subst(mask_bat, from=0, to=NA)
writeRaster(mask_bat,'input_data/mask_bat1000.tif',overwrite=T)

# 4. Shelf, slope and basin ----
## 4.1. Continental shelf ----

shelf_raw <- ifel(bat >= -200, 1, NA)
plot(shelf_raw)

### 4.1.1 Contiguity Analysis (Patching) ----
# The patches() function identifies isolated groups of contiguous pixels.
# directions = 8 means pixels touching diagonally are considered connected.
shelf_patches <- patches(shelf_raw, directions = 8)
plot(shelf_patches)
# Calculate the frequency (number of pixels) of each unique patch
# This returns a data.frame with columns: layer, value (the Patch ID), and count
patch_counts <- freq(shelf_patches)

# Find the ID of the patch with the absolute highest pixel count.
# Topologically, this massive continuous polygon is the continental shelf.
largest_patch_id <- patch_counts$value[which.max(patch_counts$count)]

shelf_clean <- ifel(shelf_patches == largest_patch_id, 1, NA)

plot(bat, main = "Clean Continental Shelf (Seamounts Removed)", col = gray.colors(100))
plot(shelf_clean, col = "cyan", add = TRUE, legend = FALSE)

not_shelf<-ifel(is.na(shelf_clean),1,NA)
not_shelf<-mask(not_shelf,mask)
plot(not_shelf)

writeRaster(shelf_clean,'output_data/shelf.tif')

## 4.2 Slope ----

### 4.2.1. Thresholding the Slope ----

# oceanic basin and the flat continental shelf.
slope_raw <- ifel(bat >= -2000, 1, NA)
slope_raw<-mask(slope_raw,not_shelf)
plot(slope_raw)

slope_patches <- patches(slope_raw, directions = 8)
plot(slope_patches)

# Calculate the frequency (number of pixels) of each unique patch
# This returns a data.frame with columns: layer, value (the Patch ID), and count
patch_counts <- freq(slope_patches)

# Find the ID of the patch with the absolute highest pixel count.
# Topologically, this massive continuous polygon is the continental slope.
largest_patch_id <- patch_counts$value[which.max(patch_counts$count)]

slope_clean <- ifel(slope_patches == largest_patch_id, 1, NA)

plot(slope_clean)

### 4.2.2 Vectorization (SpatRaster to SpatVector) ----
# The as.polygons() function converts the raster cells into vector geometry.
# dissolve = TRUE is critical here: it merges all adjacent pixels with the 
# same value into a single, continuous multipart polygon.
slope_vect <- as.polygons(slope_clean, dissolve = TRUE)
slope_sf <- st_as_sf(slope_vect)
plot(slope_sf)

write_sf(slope_sf,'output_data/slope2000.shp')

### 4.2.3. Manual editing slope polygon ----
#Slope polygon splited on QGIS to remove Vitoria-Trindade mountains from slope

slope_sf<-read_sf('output_data/slope.shp')
plot(slope_sf)

slope<-rasterize(slope_sf, mask, field = 'patches')
plot(slope)

## 4.3. Basin ----
basin_raw <- ifel(!is.na(bat) & is.na(shelf_clean) & is.na(slope), 1, NA)
plot(basin_raw)
plot(shelf_clean)

slope<-slope*2
basin<-basin_raw*3

## 4.4. Merging rasters ----
geomorphology<-cover(shelf_clean,slope)
geomorphology<-cover(geomorphology,basin)

plot(geomorphology)


# 5. Exporting the Unified Map ----
writeRaster(geomorphology, 'output_data/benthic_geomorphology.tif', overwrite = TRUE)