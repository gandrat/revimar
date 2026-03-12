## Sediment Interpolation Script
## Data from BNDO (Banco Nacional de Dados Oceanográficos)

# 1. Loading Packages ----
packages <- c('sf', 'ggplot2', 'terra', 'dplyr', 'Rsagacmd')

# Function to dynamically check, install (if missing), and load packages
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# Initialize the SAGA GIS environment via Rsagacmd
saga <- saga_gis()

# 2. Creating Projection ----
# Using SIRGAS 2000 / Brazil Albers Equal Area (EPSG:10857) [Image of Albers Equal Area Conic projection]
# Crucial for Natural Neighbors interpolation and area calculations, 
# as it preserves area proportions across the continental scale.
crs_albers_brasil <- "+proj=aea +lat_0=-12 +lon_0=-54 +lat_1=-2 +lat_2=-22 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# 3. Reading and Preprocessing Data ----

## 3.1 Points Data Setup ----
pts <- read_sf('input_data/bndo_eez/bndo_eez.shp')

# Reprojecting to Albers Equal Area
pts <- st_transform(pts, crs = crs_albers_brasil)

# Selecting, renaming, and mutating variables into a clean architecture
pts <- pts %>% transmute(
  id = num_amostr,
  lon = lon,
  lat = lat,
  time = as.POSIXct(data_hora, tz = "America/Sao_Paulo", format = '%Y/%m/%d %H:%M:%S'),
  sand = areia,
  gravel = cascalho,
  mud = argila + silte,          # Grouping fines into 'mud' for Folk classification
  tot = areia + cascalho + argila + silte
)

# Filtering out inconsistent samples (ensuring granulometric fractions sum to exactly 100%)
pts <- pts %>% filter(tot == 100)

## 3.2 Splitting Data: Interpolation vs. Accuracy Assessment ----
# Setting aside a 5% holdout dataset for validation
n_acc <- round(nrow(pts) * 0.05)

# Creating a control column
pts <- pts %>% mutate(use = 'interp')

# Randomly sampling indices for the accuracy subset
index <- sample(seq_len(nrow(pts)), size = n_acc)
pts$use[index] <- 'accuracy'

# Exporting the preprocessed and classified point data
write_sf(pts, 'input_data/bndo.shp')

# Separating the datasets based on the control column
pts_interp <- pts %>% filter(use == 'interp')
pts_acc <- pts %>% filter(use == 'accuracy')

# 4. Creating Spatial Masks ----

## 4.1 Exclusive Economic Zone (EEZ) Mask ----
eez <- read_sf('input_data/study_area.shp')
eez <- st_transform(eez, crs_albers_brasil)

# Creating a base raster grid with 1km (1000m) resolution
mask <- rast(eez, res = 1000)
mask <- rasterize(eez, mask, field = 'id')

# Ensuring the mask strictly adheres to the Albers projection and resolution
mask <- project(mask, crs_albers_brasil, res = 1000)

# writeRaster(mask, 'input_data/eez_albers.tif', overwrite = TRUE)

## 4.2 Bathymetric Model Mask ----
bat <- rast('input_data/batimetria_dhn.tif')

# Aligning bathymetry to the EEZ grid
bat <- project(bat, mask)
bat <- mask(bat, mask)
res(bat)
plot(bat, main = "Bathymetry")

# writeRaster(bat, 'input_data/batimetria_dhn_albers.tif', overwrite = TRUE)

# Creating a Boolean mask to restrict modeling to depths shallower than 1000m
# (e.g., continental shelf and upper slope)
mask_bat <- bat >= -1000
plot(mask_bat, main = "Bathymetric Mask (>= -1000m)")

# Converting 0s to NA to optimize processing in subsequent steps
mask_bat <- subst(mask_bat, from = 0, to = NA)
# writeRaster(mask_bat, 'input_data/mask_bat.tif', overwrite = TRUE)

# 5. Interpolating Sediments ----
# Using Natural Neighbors (Sibson's method) via SAGA GIS [Image of Natural Neighbor interpolation Voronoi area stealing]

## 5.1 Sand Interpolation ----
sand <- saga$grid_gridding$natural_neighbour(
  points = pts_interp,
  field = "sand",
  target_definition = 1,     # Tells SAGA to use a user-defined grid
  target_template = mask     # Provides the grid template (our 1km EEZ mask)
)

# Applying the bathymetric depth limit
sand <- mask(sand, mask_bat)
plot(sand, main = "Interpolated Sand (%)")
# writeRaster(sand, 'output_data/nn_sand.tif', overwrite = TRUE)

## 5.2 Mud Interpolation ----
mud <- saga$grid_gridding$natural_neighbour(
  points = pts_interp,         
  field = "mud",
  target_definition = 1,
  target_template = mask
)

mud <- mask(mud, mask_bat)
plot(mud, main = "Interpolated Mud (%)")
# writeRaster(mud, 'output_data/nn_mud.tif', overwrite = TRUE)

## 5.3 Gravel Interpolation ----
gravel <- saga$grid_gridding$natural_neighbour(
  points = pts_interp,         
  field = "gravel",
  target_definition = 1,
  target_template = mask
)

gravel <- mask(gravel, mask_bat)
plot(gravel, main = "Interpolated Gravel (%)")
# writeRaster(gravel, 'output_data/nn_gravel.tif', overwrite = TRUE)

# 6. Accuracy Assessment ----
# Extracting predicted values at the locations of the holdout points
acc_sand <- extract(sand, pts_acc)
names(acc_sand) <- c('id', 'sand_interp')

acc_mud <- extract(mud, pts_acc)
names(acc_mud) <- c('id', 'mud_interp')

acc_gravel <- extract(gravel, pts_acc)
names(acc_gravel) <- c('id', 'gravel_interp')

# Joining interpolated values back to the accuracy sf object
pts_acc <- pts_acc %>% mutate(
  sand_interp = acc_sand$sand_interp,
  mud_interp = acc_mud$mud_interp,
  gravel_interp = acc_gravel$gravel_interp
)

# Visualizing Error Distributions (Residuals) [Image of histogram showing error distribution]
ggplot(pts_acc, aes(x = sand - sand_interp)) + geom_histogram(bins = 30) + ggtitle("Sand Error Distribution")
ggplot(pts_acc, aes(x = mud - mud_interp)) + geom_histogram(bins = 30) + ggtitle("Mud Error Distribution")
ggplot(pts_acc, aes(x = gravel - gravel_interp)) + geom_histogram(bins = 30) + ggtitle("Gravel Error Distribution")

# Calculating statistical error metrics
accuracy <- pts_acc %>%
  st_drop_geometry() %>%
  filter(!is.na(sand_interp)) %>%
  summarise(
    # ME (Mean Error): Indicates bias (overestimation or underestimation)
    sand_me = mean(sand - sand_interp),
    # MAE (Mean Absolute Error): Average magnitude of the errors
    sand_mae = mean(abs(sand - sand_interp)),
    # RMSE (Root Mean Square Error): Penalizes larger errors
    sand_rmse = sqrt(mean((sand - sand_interp)^2)),
    
    mud_me = mean(mud - mud_interp),
    mud_mae = mean(abs(mud - mud_interp)),
    mud_rmse = sqrt(mean((mud - mud_interp)^2)),
    
    gravel_me = mean(gravel - gravel_interp),
    gravel_mae = mean(abs(gravel - gravel_interp)),
    gravel_rmse = sqrt(mean((gravel - gravel_interp)^2))
  )

print(accuracy)

# 7. Final Outputs & Workspace Saving ----

# Stacking the three matrices
sed <- c(sand, mud, gravel)
names(sed) <- c('sand', 'mud', 'gravel')

# Checking the sum of the fractions (ideally close to 100)
t <- sum(sed)
hist(t, main = "Histogram of Total Sediment Sum")
plot(t, main = "Spatial Distribution of Sediment Sum")

# Exporting the multi-band raster
writeRaster(sed, 'output_data/sediment.tif', overwrite = TRUE)

# Saving the R workspace objects for easy retrieval during the Folk classification step
save(sed, accuracy, pts_acc, pts_interp, file = 'rda/01_interp_sediment.Rda')