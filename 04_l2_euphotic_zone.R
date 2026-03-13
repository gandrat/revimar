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

# Load the base Albers rasters (Geomorphology and Bathymetry)
geomorphology <- rast('output_data/benthic_geomorphology.tif')
bat <- rast('input_data/batimetria_dhn_albers.tif')

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
zeu_med <- project(zeu_med, geomorphology)

# Clip the projected Zeu raster exactly to the boundaries of our mapped geomorphology
zeu_med <- mask(zeu_med, geomorphology)


writeRaster(zeu_med,'input_data/zeu_cmems_2024.tif',overwrite=T)

plot(zeu_med, main = "Median Euphotic Depth (m)")


# 4. Defining the Euphotic Benthos ----
# We create a boolean mask to identify where sunlight reaches the seafloor.
# Logic: If the depth light can penetrate (Zeu) is greater than or equal to the 
# depth of the seafloor (bat), then the bottom is illuminated (Euphotic = 1).
# We multiply 'bat' by -1 to turn negative depths (e.g., -40m) into positive values (40m) 
# so they can be directly compared to the positive Zeu values.
zeu_bat <- zeu_med >= (bat * -1)

plot(zeu_bat, main = "Euphotic Benthos (1 = Light Reaches Bottom)")

writeRaster(zeu_bat,'output_data/euphotic_zone.tif',overwrite=T)

# 5. Overlapping with Geomorphology (Habitat Zones) ----
# First, convert the TRUE/FALSE boolean raster into numeric 1 (Euphotic) and 0 (Aphotic)
zeu_bat_num <- ifel(zeu_bat, 1, 0)
plot(zeu_bat_num)

## 5.1 Map Algebra ----
# Combine the rasters to create unique two-digit IDs.
# Geomorphology (1 to 4) * 10 + Photic (0 or 1)
# Example: Shelf (1) * 10 + Euphotic (1) = 11 (Euphotic Shelf)
geom_photic <- (geomorphology * 10) + zeu_bat_num
plot(geom_photic)


# 6. Reclassification Matrix ----
# We create a 2-column matrix: [Old Value, New Value].
# We force all Slope pixels (20, 21) to become 20.
# We force all Basin pixels (30, 31) to become 30.
# Shelf (10, 11) and Seamounts (40, 41) remain unchanged.

reclass_matrix <- matrix(c(
  10, 10,  # Aphotic Shelf    -> Aphotic Shelf
  11, 11,  # Euphotic Shelf   -> Euphotic Shelf
  20, 20,  # Aphotic Slope    -> Continental Slope
  21, 20,  # Euphotic Slope   -> Continental Slope (Removing photic distinction)
  30, 30,  # Aphotic Basin    -> Oceanic Basin
  31, 30,  # Euphotic Basin   -> Oceanic Basin (Removing photic distinction)
  40, 40,  # Aphotic Seamount -> Aphotic Seamount
  41, 41   # Euphotic Seamount-> Euphotic Seamount
), ncol = 2, byrow = TRUE)

# Apply the reclassification
geom_photic <- classify(geom_photic, reclass_matrix)

# 3. Isolating Abyssal Seamounts (Depth < -200m) ----
geom_photic <- ifel(geom_photic %in% c(40, 41) & bat < -200, 42, geom_photic)



alpha_matrix <- matrix(c(
  11, 1,   # 1 = A1. Euphotic Shelf
  10, 2,   # 2 = A2. Mesophotic Shelf
  20, 3,   # 3 = B3. Bathial Continental Slope
  30, 4,   # 4 = C4. Oceanic Basin
  41, 5,   # 5 = D1. Euphotic Seamount
  40, 6,   # 6 = D2. Mesophotic Seamount
  42, 7    # 7 = D3. Bathial Seamount
), ncol = 2, byrow = TRUE)

geom_photic_id <- classify(geom_photic, alpha_matrix)
plot(geom_photic_id)

# Convert to categorical raster
geom_photic_id <- as.factor(geom_photic_id)

# Create the Raster Attribute Table (RAT) for the combined habitats
# Note: Euphotic Basin (31) is highly unlikely due to depth, but included for mathematical completeness.
# Create the new, simplified Raster Attribute Table (RAT)
habitat_table <- data.frame(
  ID = 1:7,
  Habitat_Zone = c(
    "A1. Euphotic Shelf",
    "A2. Mesophotic Shelf",
    "B3. Bathial Slope",
    "C4. Oceanic Basin",
    "D1. Euphotic Seamount",
    "D2. Mesophotic Seamount",
    "D3. Bathial Seamount"
  )
)

# Apply the exact names to the new raster levels
geom_photic<-geom_photic_id
levels(geom_photic) <- habitat_table
plot(geom_photic)
## 6.1. Defining the Color Palette ----
# We map the exact colors to the specific IDs in the raster.
# ID 10 = Aphotic Shelf, 11 = Euphotic Shelf
# ID 20 = Continental Slope, 30 = Oceanic Basin
# ID 40 = Aphotic Seamount, 41 = Euphotic Seamount

# We create a data.frame where column 1 is the pixel ID and column 2 is the color (Hex or R name)
color_table <- data.frame(
  value = 1:7,
  color = c(
    "#00FFFF",  # 1: A1. Euphotic Shelf (Cyan)
    "#00688B",  # 2: A2. Mesophotic Shelf (Deep Sky Blue 4)
    "#FFA500",  # 3: B3. Bathial Continental Slope (Orange)
    "#000050",  # 4: C4. Oceanic Basin (Navy)
    "#FF4500",  # 5: D1. Euphotic Seamount (Orange Red)
    "#8B0000",  # 6: D2. Mesophotic Seamount (Dark Red)
    "#4A4000"   # 7: D3. Bathial Seamount (Indigo)
  )
)

coltab(geom_photic) <- color_table


# Visualization
geom_photic<-rast('output_data/geomorphology_photic_zones_v2.tif')

jpeg(filename = "figures/benthic_geomorphology_l2.jpg", 
     width = 18,       # Width of the image
     height = 17,       # Height of the image
     units = "cm",     # Units for width/height (inches)
     res = 300)        # Resolution in pixels per inch
plot(geom_photic, main = "Photic Zone (L2)",
     mar = c(3, 3, 3, 17))
dev.off()

# 7. Statistical Summary and Export (km²) ----

# Use freq() on the final combined raster 
photic_summary <- freq(geom_photic)
photic_summary$value<-row_number(photic_summary)

print(photic_summary)

## 7.1. Plot areas ----
# Create a reference table with our exact names and colors
habitat_reference <- data.frame(
  value = 1:7,
  Habitat_Zone = c(
    "A1. Euphotic Shelf",
    "A2. Mesophotic Shelf",
    "B3. Bathial Slope",
    "C4. Oceanic Basin",
    "D1. Euphotic Seamount",
    "D2. Mesophotic Seamount",
    "D3. Bathial Seamount"
  ),
  color = c(
    "#00FFFF",  # 1: A1. Euphotic Shelf (Cyan)
    "#00688B",  # 2: A2. Mesophotic Shelf (Deep Sky Blue 4)
    "#FFA500",  # 3: B3. Bathial Continental Slope (Orange)
    "#000050",  # 4: C4. Oceanic Basin (Navy)
    "#FF4500",  # 5: D1. Euphotic Seamount (Orange Red)
    "#8B0000",  # 6: D2. Mesophotic Seamount (Dark Red)
    "#4A4000"   # 7: D3. Bathial Seamount (Indigo)
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
ggsave('figures/benthic_habitat_level2.jpg',width = 20, height = 12, dpi=150, units = 'cm')


# 8. Export data ----
# Export the final combined raster
writeRaster(geom_photic, 'output_data/geomorphology_photic_zones_v2.tif', overwrite = TRUE)

# Export the summary table as a CSV for your technical report
write.csv(photic_summary, 'output_data/habitat_area_summary.csv', row.names = FALSE)


