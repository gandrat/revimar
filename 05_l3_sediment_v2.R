## Final Benthic Habitat Integration & Cartography Script
## Integrates 4-Class Sediment (with Biogenic) and Excludes Deep Water Sediment
# 1. Loading Packages and Data ----
packages <- c('terra', 'dplyr','ggplot2','scales')

# Dynamically check, install, and load required packages
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# Load the photic geomorphology zones
geom_photic <- rast('output_data/l2_geomorphology_photic_zones_v2.tif')
plot(geom_photic)
# Load the 4-class sediment raster (Ensure IDs: 1=Muddy, 2=Sandy, 3=Gravel, 4=Biogenic)
sediments <- rast('output_data/l0_folk_classification_biogenic4.tif') 
plot(sediments)

# 2. Map Algebra (The Intersection) ----
# Shelf (1, 2) and shallow Seamounts (5, 6) get sediment IDs added (1 to 4).
# Slope (3), Basin (4), and Bathyal Seamount (7) get multiplied by 100 only.
habitats <- ifel(
  geom_photic %in% c(1, 2, 5, 6), 
  (geom_photic * 100) + sediments, 
  geom_photic * 100
)

plot(habitats)
# 3. Generating the Attribute Tables ----
# Define the structural zones that will receive sediment data
base_zones <- data.frame(
  Zone_ID = c(1, 2, 5, 6),
  Prefix = c("A1", "A2", "D1", "D2"),
  Zone_Name = c("Euphotic Shelf", "Mesophotic Shelf", "Euphotic Seamount", "Mesophotic Seamount")
)

# Define the new 4-class sediment order (a, b, c, d)
folk_classes <- data.frame(
  Sed_ID = 1:4,
  Letter = letters[1:4], 
  Sed_Name = c("Muddy", "Sandy", "Gravel", "Biogenic")
)

# Create combinations exclusively for zones with sediment data
combinations <- expand.grid(Zone_ID = c(1, 2, 5, 6), Sed_ID = 1:4) %>%
  left_join(base_zones, by = "Zone_ID") %>%
  left_join(folk_classes, by = "Sed_ID") %>%
  mutate(
    ID = (Zone_ID * 100) + Sed_ID,
    Habitat_Name = paste0(Prefix, Letter, ". ", Sed_Name, " ", Zone_Name)
  ) %>%
  select(ID, Habitat_Name)

# Define the deep-water zones without sediment data
zones_without_sediment <- data.frame(
  ID = c(300, 400, 700),
  Habitat_Name = c("B3. Bathial Continental Slope", "C4. Oceanic Basin", "D3. Bathial Seamount")
)

# Bind everything together and sort alphabetically to build the theoretical table
habitat_table <- bind_rows(combinations, zones_without_sediment) %>%
  arrange(Habitat_Name)


# 4. Generating the Color Palette ----
# Base colors for the structural zones
geo_colors <- data.frame(
  Zone_ID = c(1, 2, 5, 6), 
  Geo_Hex = c("#00FFFF", "#00688B", "#FF4500", "#8B0000")
)

# Base colors for the new 4 sediment classes
sed_colors <- data.frame(
  Sed_ID = 1:4, 
  Sed_Hex = c(
    "#4A7A40", # 1: Muddy (Dark Green)
    "#FCE47F", # 2: Sandy (Yellow)
    "#B95246", # 3: Gravel (Reddish Brown)
    "#DDA0DD"  # 4: Biogenic (Plum / Light Purple)
  )
)

# Function to mathematically mix two hex colors
mix_colors <- function(c1, c2, weight = 0.6) {
  rgb1 <- col2rgb(c1)
  rgb2 <- col2rgb(c2)
  mixed <- round(rgb1 * (1 - weight) + rgb2 * weight)
  rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

# Blend the colors for the 16 sediment-bearing classes
habitat_colors <- expand.grid(Zone_ID = c(1, 2, 5, 6), Sed_ID = 1:4) %>%
  left_join(geo_colors, by = "Zone_ID") %>%
  left_join(sed_colors, by = "Sed_ID") %>%
  rowwise() %>%
  mutate(
    value = (Zone_ID * 100) + Sed_ID, 
    color = mix_colors(Geo_Hex, Sed_Hex, weight = 0.6)
  ) %>%
  ungroup() %>%
  select(value, color)

# Add the standalone deep-water colors
deep_water_colors <- data.frame(
  value = c(300, 400, 700),
  color = c("#FFA500", "#000050", "#4A4000")
)

# The full theoretical color table
final_color_table <- bind_rows(habitat_colors, deep_water_colors) %>%
  arrange(value)


# 5. Filtering and Applying Levels ----
# Convert to factor first so terra builds its internal index of existing values
habitats <- as.factor(habitats)

# Extract the exact active Raster Attribute Table (RAT) that terra just created
active_rat <- levels(habitats)[[1]]
colnames(active_rat)[1] <- "ID"

# Merge our Habitat Names ONLY onto the true IDs that exist in the map
updated_rat <- active_rat %>%
  left_join(habitat_table, by = "ID") %>%
  as.data.frame()
# Re-apply the perfectly matched names

# Re-apply the perfectly matched names
levels(habitats) <- updated_rat

# Do the exact same foolproof join for the Color Table
safe_color_table <- data.frame(value = active_rat$ID) %>%
  left_join(final_color_table, by = "value") %>%
  as.data.frame()

# Apply the perfectly matched colors
coltab(habitats) <- safe_color_table

plot(habitats)


# 6. Visualization & Export ----
jpeg(filename = "figures/l3_benthic_habitats.jpg", 
     width = 40,        # Width of the image in cm
     height = 50,       # Height of the image in cm
     units = "cm",      
     res = 300)         # Resolution in pixels per inch

# Plotting with layout adjustments for long legend names
# Plotting with layout adjustments and explicitly calling the labels
plot(habitats, 
     main = "Benthic Habitats (L3)",
     type = "classes",  
     col = final_rat$color,  
     levels = updated_rat$Habitat_Name,          # Injeta os nomes corretos na legenda
     mar = c(3, 3, 3, 25),
     cex.main = 3,
     plg = list(cex = 2))                     

dev.off()

# Export the true, filtered maps and legends
# Using datatype INT2U to support IDs greater than 255 (prevents truncation warnings)
writeRaster(habitats, 'output_data/l3_benthic_habitats_v1.tif', 
            overwrite = TRUE)
write.csv(updated_rat, 'output_data/l3_benthic_habitats_legend.csv', row.names = FALSE)




# 6. Safe Exporting (8-bit GeoTIFF with Embedded Color Table) ----

# Step A: Convert the factor raster back to raw numbers (stripping the categorical metadata)
hab_numeric <- as.numeric(habitats)

# Step B: Create a strict reclassification matrix
# Column 1 = Old high IDs (101, 300, 700...)
# Column 2 = New sequential 8-bit IDs (1, 2, 3...)
final_rat$Export_ID <- 1:nrow(final_rat)
rcl_matrix <- as.matrix(final_rat[, c("ID", "Export_ID")])

# Step C: Reclassify! 
# The 'others = NA' argument is the safety net. It guarantees that any stray pixel 
# not in our master table is converted to NA, physically preventing the INT1U warning.
habitats_export <- classify(hab_numeric, rcl_matrix, others = NA)

# Step D: Convert back to a categorical factor and apply the NEW metadata
habitats_export <- as.factor(habitats_export)

# Apply the names mapped to the new Export_ID
levels(habitats_export) <- final_rat[, c("Export_ID", "Habitat_Name")]

# Apply the colors mapped to the new Export_ID
coltab(habitats_export) <- final_rat[, c("Export_ID", "color")]

# Step E: Final Export
# With absolutely no values above 255 remaining, INT1U will execute perfectly,
# and the color palette will be permanently embedded in the GeoTIFF.
writeRaster(habitats_export, 
            'output_data/l3_benthic_habitats_v1.tif', 
            datatype = "INT1U", 
            overwrite = TRUE)
