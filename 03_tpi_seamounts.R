## Topographic Position Index (TPI) and Seamount Detection Script
## Using terra::terrain() function

# 1. Loading Packages and Data ----

packages<-c('terra')

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# Loading the full bathymetry in Albers Equal Area (EPSG:10857)
bat <- rast('input_data/batimetria_dhn_albers.tif')

## Topographic Position Index (TPI) and Seamount Detection Script
## Requires the unmasked bathymetry raster (Albers projection)

# 1. Loading Packages and Data ----
library(terra)

# Loading the full bathymetry in Albers Equal Area (EPSG:10857)
# Note: We use the unmasked version because seamounts often rise from the deep ocean floor
bat <- rast('input_data/batimetria_dhn_albers.tif')


# 2. Defining the Neighborhood Scale ----
# The radius defines what size of feature you are looking for.
# Since Albers is measured in meters, a 25000 radius equals 25 km.
# You will likely need to tweak this value based on the scale of your target seamounts.
search_radius_m <- 25000 

# Create a circular weight matrix for the focal function
# Pixels inside the circle get a weight of 1, outside get 0
focal_matrix <- focalMat(bat, d = search_radius_m, type = "circle")

# Convert the weights to a Boolean matrix (1 and NA) to calculate a true mean
# without shrinking the values by the matrix proportions
focal_matrix[focal_matrix > 0] <- 1
focal_matrix[focal_matrix == 0] <- NA


# 3. Calculating TPI ----
# Calculate the mean bathymetry of the surrounding neighborhood for each pixel
# na.rm = TRUE ensures edges and coastal boundaries don't return NA
neighborhood_mean <- focal(bat, w = focal_matrix, fun = mean, na.rm = TRUE)
plot(neighborhood_mean)
# TPI is simply the central pixel's depth minus the neighborhood mean depth
tpi <- bat - neighborhood_mean

plot(tpi, main = "Topographic Position Index (TPI)")


# 4. Identifying Seamount Tops (Thresholding) ----
# A standard geomorphological approach is to classify TPI based on its Standard Deviation (SD).
# Peaks/Seamounts are usually defined as areas where TPI is > 1 Standard Deviation above the mean.

# Calculate the global standard deviation of the TPI raster
tpi_sd <- global(tpi, fun = "sd", na.rm = TRUE)[1, 1]

# Create a Boolean mask: 1 for Seamount Tops, NA for everything else
seamount_tops <- ifel(tpi > 200, 1, NA)

# Ensure it is treated as categorical for mapping
seamount_tops <- as.factor(seamount_tops)
levels(seamount_tops) <- data.frame(ID = 1, Feature = "Seamount Top")

# Visualization: Overlaying the peaks on top of the original bathymetry
plot(bat, main = "Identified Seamount Tops", col = gray.colors(100))
plot(seamount_tops, col = "red", add = TRUE, legend = FALSE)


# 5. Exporting Results ----
writeRaster(tpi, 'output_data/tpi_continuous.tif', overwrite = TRUE)
writeRaster(seamount_tops, 'output_data/seamount_tops_binary.tif', overwrite = TRUE)