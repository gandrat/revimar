## Light Penetration (Euphotic Zone) Script
## Deriving benthic light availability from CMEMS NetCDF data and combining with Geomorphology

# 1. Loading Packages and Data ----
packages <- c('terra', 'sf','dplyr','ggplot2','scales')

# Dynamically check, install, and load required packages
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# 1. Load the base rasters (Provinces and Bathymetry) ----
crs_albers_brasil <- "+proj=aea +lat_0=-12 +lon_0=-54 +lat_1=-2 +lat_2=-22 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

l1 <- rast('output_data/l1_benthic_provinces_v01.tif')
plot(l1, main = "L1 Benthic Provinces")

bat <- rast('output_data/bathymetry_gebco_v01.tif')


# 2. Temporal Aggregation of Zeu ----
# # Load the Copernicus CMEMS Euphotic Depth (Zeu) NetCDF file.
# zeu <- rast('input_data/cmems_mod_glo_bgc_my_0.083deg-lmtl_P1D-i_1773341382551.nc')
# 
# # Calculate the cell-wise median across all time layers to get a climatological baseline
# zeu_med <- median(zeu, na.rm = TRUE)
# 
# # 3. Spatial Alignment (Resampling and Cropping) ----
# # Project to Albers Equal Area and mask to match bathymetry exactly
# zeu_med <- project(zeu_med, l1)
# zeu_med <- mask(zeu_med, bat)
# writeRaster(zeu_med, 'input_data/zeu_cmems_2024_v2.tif', overwrite = TRUE)
# 
# # 4. Defining the Euphotic Benthos ----
# # Create a boolean mask to identify where sunlight reaches the seafloor.
# zeu_bat <- zeu_med >= (bat * -1)
# plot(zeu_bat, main = "Light Reaches Bottom (Euphotic)")
# writeRaster(zeu_bat, 'output_data/l0_euphotic_zone.tif', overwrite = TRUE)


zeu_bat<-rast('output_data/l0_euphotic_zone.tif')
# Convert the TRUE/FALSE boolean raster into numeric 1 (Euphotic) and 0 (Aphotic)
zeu_bat_num <- ifel(zeu_bat, 1, 0)

# 5. Overlapping with Provinces (Habitat Zones) ----
# Combine the rasters to create unique two-digit IDs.
# Geomorphology L1 (1 to 5) * 10 + Photic Zone (0 or 1)
# 1=Shelf, 2=Slope, 3=Rise, 4=Basin, 5=Seamount
geom_photic <- (l1 * 10) + zeu_bat_num

# 6. Reclassification and Deep Water Stratification ----

## 6.1 Initial Photic Clean-up
# Force Slope (20/21), Rise (30/31), and Basin (40/41) to drop the euphotic distinction.
reclass_matrix <- matrix(c(
  10, 10,  # Aphotic Shelf    -> Aphotic Shelf
  11, 11,  # Euphotic Shelf   -> Euphotic Shelf
  20, 20,  # Aphotic Slope    -> Continental Slope
  21, 20,  # Euphotic Slope   -> Continental Slope
  30, 30,  # Aphotic Rise     -> Continental Rise
  31, 30,  # Euphotic Rise    -> Continental Rise
  40, 40,  # Aphotic Basin    -> Oceanic Basin
  41, 40,  # Euphotic Basin   -> Oceanic Basin
  50, 50,  # Aphotic Seamount -> Aphotic Seamount
  51, 51   # Euphotic Seamount-> Euphotic Seamount
), ncol = 2, byrow = TRUE)

geom_photic <- classify(geom_photic, reclass_matrix)

## 6.2 Isolating Abyssal vs Bathyal Zones (Threshold: -3500m)
# We separate the deep features into Bathyal (> -2000m) and Abyssal (<= -2000m)
# 20 = Bathyal Slope, 22 = Abyssal Slope
geom_photic <- ifel(geom_photic == 20 & bat <= -3500, 22, geom_photic)

# 30 = Bathyal Rise, 32 = Abyssal Rise
geom_photic <- ifel(geom_photic == 30 & bat <= -3500, 32, geom_photic)

# 40 = Bathyal Basin, 42 = Abyssal Basin
geom_photic <- ifel(geom_photic == 40 & bat <= -3500, 42, geom_photic)

## 6.3 Isolating Deep Seamounts (Depth < -200m)
geom_photic <- ifel(geom_photic %in% c(50, 51) & bat < -200, 52, geom_photic)


# 7. Sequential ID Reclassification Matrix ----
# Map the logical IDs to a sequential 1 to 11 ID for the final legend
alpha_matrix <- matrix(c(
  11, 1,   
  10, 2,  
  20, 3,  
  22, 4,  
  30, 5,   
  32, 6,   
  40, 7,   # Merge with D4. Abyssal Oceanic Basin
  42, 7,   
  51, 8,   
  50, 9,  
  52, 10   
), ncol = 2, byrow = TRUE)

geom_photic_id <- classify(geom_photic, alpha_matrix)

# 8. Raster Attribute Table (RAT) ----
geom_photic_id <- as.factor(geom_photic_id)

habitat_table <- data.frame(
  ID = 1:10,
  Habitat_Zone = c(
    "A1. Euphotic Shelf",
    "A2. Mesophotic Shelf",
    "B3. Bathyal Continental Slope",
    "B4. Abyssal Continental Slope",
    "C3. Bathyal Continental Rise",
    "C4. Abyssal Continental Rise",
    "D4. Abyssal Oceanic Basin",
    "E1. Euphotic Seamount",
    "E2. Mesophotic Seamount",
    "E3. Bathyal Seamount"
  )
)

levels(geom_photic_id) <- habitat_table

# 9. Defining the Color Palette ----
color_table <- data.frame(
  value = 1:10,
  color = c(
    "#00FFFF",  # 1: A1. Euphotic Shelf (Cyan)
    "#00688B",  # 2: A2. Mesophotic Shelf (Deep Sky Blue)
    "#FFA500",  # 3: B3. Bathyal Slope (Orange)
    "#FF8C00",  # 4: B4. Abyssal Slope (Dark Orange)
    "#FFD700",  # 5: C5. Bathyal Rise (Gold)
    "#B8860B",  # 6: C6. Abyssal Rise (Dark Goldenrod)
    "#000050",  # 8: C8. Abyssal Basin (Navy)
    "#FF4500",  # 9: D1. Euphotic Seamount (Orange Red)
    "#b50000",  # 10: D2. Mesophotic Seamount (Dark Red)
    "#8c3100"   # 11: D3. Bathyal/Abyssal Seamount (Indigo)
  )
)

coltab(geom_photic_id) <- color_table
plot(geom_photic_id, main = "Benthic Habitats of the Brazilian Margin (L2)")

# 10. Statistical Summary and Export (km²) ----
photic_summary <- as.data.frame(freq(geom_photic_id))
# Ensure the value column is treated as integer for joining
photic_summary$value <- 1:10

## 10.1. Plot areas ----
habitat_reference <- data.frame(
  value = 1:10,
  Habitat_Zone = habitat_table$Habitat_Zone,
  color = color_table$color
)

plot_data <- photic_summary %>%
  left_join(habitat_reference, by = "value") %>%
  select(value, count = count, color, Habitat_Zone)

# Preserve order for the chart (reversed so A1 is on top when flipped)
plot_data$Habitat_Zone <- factor(plot_data$Habitat_Zone, levels = rev(habitat_reference$Habitat_Zone))
custom_colors <- setNames(plot_data$color, plot_data$Habitat_Zone)

## 10.2. Generating the Bar Chart (ggplot2) ----
ggplot(plot_data, aes(x = Habitat_Zone, y = count / sum(count), fill = Habitat_Zone, color = Habitat_Zone)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Photic and Bathymetric Zones (L2)", y = "Relative Area", x = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave('figures/l2_plot_area_v01.jpg', width = 20, height = 15, dpi = 150, units = 'cm', bg = 'white')

# 11. Export data ----
writeRaster(geom_photic_id, 
            'output_data/l2_province_photic_zones_v01.tif', 
            datatype = "INT1U",
            overwrite = TRUE)

write.csv(photic_summary, 'output_data/habitat_area_summary_v01.csv', row.names = FALSE)
