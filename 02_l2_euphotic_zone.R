## Light Penetration (Euphotic Zone) Script
## Deriving benthic light availability from CMEMS NetCDF data

# 1. Loading Packages and Data ----
packages <- c('terra', 'sf', 'dplyr','ggplot2','scales')

# Dynamically check, install, and load required packages
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

#1. Load the base rasters (Provinces and Bathymetry)---------

crs_albers_brasil <- "+proj=aea +lat_0=-12 +lon_0=-54 +lat_1=-2 +lat_2=-22 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

l1 <- rast('output_data/l1_benthic_provinces_v2.tif')
plot(l1)
res(l1)
bat <- rast('output_data/bathymetry_gebco.tif')
plot(bat)
res(bat)

bat<-project(bat,l1)

ss<-read_sf('input_data/study_area.shp')
ss<-st_transform(ss,crs_albers_brasil)


bat<-mask(bat,ss)

plot(bat)


# Load the Copernicus CMEMS Euphotic Depth (Zeu) NetCDF file.
# This file likely contains a time series (multiple layers/bands) of daily or monthly data.
zeu <- rast('input_data/cmems_mod_glo_bgc_my_0.083deg-lmtl_P1D-i_1773341382551.nc')


# 2. Temporal Aggregation ----
# Because the NetCDF 'zeu' has multiple layers (time steps), using median() here 
# does not calculate one single global number. Instead, it calculates the cell-wise 
# median across all time layers. This gives you a single raster layer representing 
# the median climatological light penetration for each pixel over the entire period.
zeu_med <- median(zeu, na.rm = TRUE)
plot(zeu_med)

# 3. Spatial Alignment (Resampling and Cropping) ----
# The CMEMS NetCDF is in unprojected geographic coordinates (Lat/Lon) with a ~9km resolution.
# We must project it to match our Albers Equal Area 1km grid.
zeu_med <- project(zeu_med, l1)

# Clip the projected Zeu raster exactly to the boundaries of our mapped geomorphology
zeu_med <- mask(zeu_med, bat)


writeRaster(zeu_med,'input_data/zeu_cmems_2024_v2.tif',overwrite=T)


# 4. Defining the Euphotic Benthos ----
# We create a boolean mask to identify where sunlight reaches the seafloor.
# Logic: If the depth light can penetrate (Zeu) is greater than or equal to the 
# depth of the seafloor (bat), then the bottom is illuminated (Euphotic = 1).
# We multiply 'bat' by -1 to turn negative depths (e.g., -40m) into positive values (40m) 
# so they can be directly compared to the positive Zeu values.


zeu_med<-rast('input_data/zeu_cmems_2024.tif')
plot(zeu_med, main = "Median Euphotic Depth (m)")
res(zeu_med)
zeu_med<-project(zeu_med,bat)

zeu_bat <- zeu_med >= (bat * -1)
plot(bat * -1)
plot(zeu_med)
plot(zeu_bat, main = "Light Reaches Bottom)")

writeRaster(zeu_bat,'output_data/l0_euphotic_zone.tif',overwrite=T)

# 5. Overlapping with Geomorphology (Habitat Zones) ----
# First, convert the TRUE/FALSE boolean raster into numeric 1 (Euphotic) and 0 (Aphotic)
zeu_bat_num <- ifel(zeu_bat, 1, 0)
plot(zeu_bat_num)

## 5.1 Map Algebra ----
# Combine the rasters to create unique two-digit IDs.
# Geomorphology L1 (1 to 5) * 10 + Photic Zone (0 or 1)
# 1=Shelf, 2=Slope, 3=Rise, 4=Basin, 5=Seamount
geom_photic <- (l1 * 10) + zeu_bat_num
plot(geom_photic, main = "Raw Two-Digit Combinations")

## 6. Reclassification Matrix ----
# We create a 2-column matrix: [Old Value, New Value].
# We force Slope (21), Rise (31), and Basin (41) to drop the euphotic distinction, 
# as these morphological features occur below the photic zone.
reclass_matrix <- matrix(c(
  10, 10,  # Aphotic Shelf    -> Aphotic Shelf
  11, 11,  # Euphotic Shelf   -> Euphotic Shelf
  20, 20,  # Aphotic Slope    -> Continental Slope
  21, 20,  # Euphotic Slope   -> Continental Slope (Removing photic distinction)
  30, 30,  # Aphotic Rise     -> Continental Rise
  31, 30,  # Euphotic Rise    -> Continental Rise (Removing photic distinction)
  40, 40,  # Aphotic Basin    -> Oceanic Basin
  41, 40,  # Euphotic Basin   -> Oceanic Basin (Removing photic distinction)
  50, 50,  # Aphotic Seamount -> Aphotic Seamount
  51, 51   # Euphotic Seamount-> Euphotic Seamount
), ncol = 2, byrow = TRUE)

# Apply the initial reclassification
geom_photic <- classify(geom_photic, reclass_matrix)

## 6.1 Isolating Deep Seamounts (Depth < -200m) ----
# Creates a 3rd Seamount class for isolated peaks in deep bathyal/abyssal waters
geom_photic <- ifel(geom_photic %in% c(50, 51) & bat < -200, 52, geom_photic)

## 7. Sequential ID Reclassification Matrix ----
# Map the two-digit logical IDs to a sequential 1 to 8 ID for the final legend
alpha_matrix <- matrix(c(
  11, 1,   # 1 = A1. Euphotic Shelf
  10, 2,   # 2 = A2. Mesophotic Shelf
  20, 3,   # 3 = B3. Bathyal Continental Slope
  30, 4,   # 4 = C4. Continental Rise
  40, 5,   # 5 = C5. Oceanic Basin
  51, 6,   # 6 = D1. Euphotic Seamount
  50, 7,   # 7 = D2. Mesophotic Seamount
  52, 8    # 8 = D3. Bathyal Seamount
), ncol = 2, byrow = TRUE)

geom_photic_id <- classify(geom_photic, alpha_matrix)

## 8. Raster Attribute Table (RAT) ----
# Convert the raster to a categorical factor
geom_photic_id <- as.factor(geom_photic_id)

# Create the Raster Attribute Table for the combined habitats
habitat_table <- data.frame(
  ID = 1:8,
  Habitat_Zone = c(
    "A1. Euphotic Shelf",
    "A2. Mesophotic Shelf",
    "B3. Bathyal Continental Slope",
    "C4. Continental Rise",
    "C5. Oceanic Basin",
    "D1. Euphotic Seamount",
    "D2. Mesophotic Seamount",
    "D3. Bathyal Seamount"
  )
)

# Apply the exact names to the raster levels
levels(geom_photic_id) <- habitat_table

## 9. Defining the Color Palette ----
# We map the exact colors to the specific IDs (1 to 8) in the categorical raster.
color_table <- data.frame(
  value = 1:8,
  color = c(
    "#00FFFF",  # 1: A1. Euphotic Shelf (Cyan)
    "#00688B",  # 2: A2. Mesophotic Shelf (Deep Sky Blue 4)
    "#FFA500",  # 3: B3. Bathyal Continental Slope (Orange)
    "#FFD700",  # 4: C4. Continental Rise (Gold - representing transitional sediment)
    "#000050",  # 5: C5. Oceanic Basin (Navy)
    "#FF4500",  # 6: D1. Euphotic Seamount (Orange Red)
    "#8B0000",  # 7: D2. Mesophotic Seamount (Dark Red)
    "#4A4000"   # 8: D3. Bathyal Seamount (Indigo/Dark Brown)
  )
)

coltab(geom_photic_id) <- color_table

# Final visual check
plot(geom_photic_id, main = "Benthic Habitats of the Brazilian Margin")

writeRaster(geom_photic, 'output_data/l2_province_photic_zones.tif')

# Visualization----------
geom_photic<-rast('output_data/l2_province_photic_zones.tif')
plot(geom_photic)
jpeg(filename = "figures/l2_map_geomorphology_photic_v2.jpg", 
     width = 40,       # Width of the image
     height = 50,       # Height of the image
     units = "cm",     # Units for width/height (inches)
     res = 300)        # Resolution in pixels per inch
plot(geom_photic, main = "Photic Zone (L2)",
     mar = c(3, 3, 3, 25),
     cex.main = 3,
     plg = list(cex = 2))
dev.off()

# 7. Statistical Summary and Export (km²) ----

# Use freq() on the final combined raster 
photic_summary <- freq(geom_photic)
photic_summary$value<-row_number(photic_summary)

print(photic_summary)

## 7.1. Plot areas ----
# Create a reference table with our exact names and colors
# Create a reference table with our exact names and colors
habitat_reference <- data.frame(
  value = 1:8,
  Habitat_Zone = c(
    "A1. Euphotic Shelf",
    "A2. Mesophotic Shelf",
    "B3. Bathyal Continental Slope",
    "C4. Continental Rise",
    "C5. Oceanic Basin",
    "D1. Euphotic Seamount",
    "D2. Mesophotic Seamount",
    "D3. Bathyal Seamount"
  ),
  color = c(
    "#00FFFF",  # 1: A1. Euphotic Shelf (Cyan)
    "#00688B",  # 2: A2. Mesophotic Shelf (Deep Sky Blue 4)
    "#FFA500",  # 3: B3. Bathyal Continental Slope (Orange)
    "#FFD700",  # 4: C4. Continental Rise (Gold)
    "#000050",  # 5: C5. Oceanic Basin (Navy)
    "#FF4500",  # 6: D1. Euphotic Seamount (Orange Red)
    "#8B0000",  # 7: D2. Mesophotic Seamount (Dark Red)
    "#4A4000"   # 8: D3. Bathyal Seamount (Indigo)
  )
)

# Merge the pixel counts with our reference table
photic_summary$value<-as.integer(photic_summary$value)

plot_data <- photic_summary %>%
  left_join(habitat_reference, by = "value") %>%
  select(value, count = count, color, Habitat_Zone)

# Create a named vector for ggplot's manual color scale
custom_colors <- setNames(plot_data$color, plot_data$Habitat_Zone)


## 7.2. Generating the Bar Chart (ggplot2) ----
ggplot(plot_data, aes(x = Habitat_Zone, y = count / sum(count), fill = Habitat_Zone, color = Habitat_Zone)) +
  geom_bar(stat = "identity") +
  coord_flip()+
  # Apply our exact cartographic colors
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  
    # Add labels and titles
  labs(
    title = "Photic Zone (L2)",
    y = NULL,
    x=NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0,1)) +
  
  # Apply a clean, professional theme
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none", # Hide legend since X-axis already has the names
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.y = element_blank(), # Remove vertical grid lines for a cleaner look
    panel.grid.major.x = element_blank()
  )
ggsave('figures/l2_plot_area_v2.jpg',width = 20, height = 12, dpi=150, units = 'cm',bg='white')



# 8. Export data ----
# Export the final combined raster
writeRaster(geom_photic_id, 
            'output_data/l2_province_photic_zones.tif', 
            datatype = "INT1U",
            overwrite = TRUE)

# Export the summary table as a CSV for your technical report
write.csv(photic_summary, 'output_data/habitat_area_summary_v2.csv', row.names = FALSE)


