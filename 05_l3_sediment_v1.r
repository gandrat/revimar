## Final Benthic Habitat Integration & Cartography Script
## Integrates 4-Class Sediment (with Biogenic) and Excludes Deep Water Sediment

# 1. Loading Packages and Data ----
packages <- c('terra', 'dplyr', 'ggplot2', 'scales')

# Dynamically check, install, and load required packages
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# Load the photic geomorphology zones (L2) - 11 classes
geom_photic <- rast('output_data/l2_province_photic_zones.tif')
plot(geom_photic, main = "Input: Geomorphology + Photic (L2)")

# Load the 4-class sediment raster (Ensure IDs: 1=Muddy, 2=Sandy, 3=Gravel, 4=Biogenic)
sediments <- rast('output_data/l0_folk_classification_biogenic4_v1.tif') 
sediments <- project(sediments, geom_photic)

# 2. Map Algebra (The Intersection) ----
# Shelf (1, 2) get sediment IDs added (1 to 4).
# Others (3 to 11) get multiplied by 100 only, receiving no sediment suffix.
habitats <- ifel(
  geom_photic %in% c(1, 2), 
  (geom_photic * 100) + sediments, 
  geom_photic * 100
)

plot(habitats, main = "Raw Map Algebra L3")

# 3. Generating the Attribute Tables ----
# Define the structural zones that will receive sediment data
base_zones <- data.frame(
  Zone_ID = c(1, 2),
  Prefix = c("A1", "A2"),
  Zone_Name = c("Euphotic Shelf", "Mesophotic Shelf")
)

# Define the new 4-class sediment order (a, b, c, d)
folk_classes <- data.frame(
  Sed_ID = 1:4,
  Letter = letters[1:4], 
  Sed_Name = c("Muddy", "Sandy", "Gravel", "Biogenic")
)

# Create combinations exclusively for zones with sediment data
combinations <- expand.grid(Zone_ID = c(1, 2), Sed_ID = 1:4) %>%
  left_join(base_zones, by = "Zone_ID") %>%
  left_join(folk_classes, by = "Sed_ID") %>%
  mutate(
    ID = (Zone_ID * 100) + Sed_ID,
    Habitat_Name = paste0(Prefix, Letter, ". ", Sed_Name, " ", Zone_Name)
  ) %>%
  select(ID, Habitat_Name)

# Define the 9 deep-water zones without sediment data (IDs 300 to 1100)
zones_without_sediment <- data.frame(
  ID = c(300, 400, 500, 600, 700, 800, 900, 1000, 1100),
  Habitat_Name = c(
    "B3. Bathyal Continental Slope",
    "B4. Abyssal Continental Slope",
    "C3. Bathyal Continental Rise",
    "C4. Abyssal Continental Rise",
    "D3. Bathyal Oceanic Basin",
    "D4. Abyssal Oceanic Basin",
    "E1. Euphotic Seamount",
    "E2. Mesophotic Seamount",
    "E3. Bathyal Seamount"
  )
)

# Bind everything together and sort to build the theoretical table
habitat_table <- bind_rows(combinations, zones_without_sediment) %>%
  arrange(ID)

# 4. Generating the Color Palette ----
# Base colors for the structural zones
geo_colors <- data.frame(
  Zone_ID = c(1, 2), 
  Geo_Hex = c("#00FFFF", "#00688B")
)

# Base colors for the 4 sediment classes
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

# Blend the colors for the 8 sediment-bearing classes (2 zones * 4 sediments)
habitat_colors <- expand.grid(Zone_ID = c(1, 2), Sed_ID = 1:4) %>%
  left_join(geo_colors, by = "Zone_ID") %>%
  left_join(sed_colors, by = "Sed_ID") %>%
  rowwise() %>%
  mutate(
    ID = (Zone_ID * 100) + Sed_ID, 
    color = mix_colors(Geo_Hex, Sed_Hex, weight = 0.6)
  ) %>%
  ungroup() %>%
  select(ID, color)

# Add the standalone deep-water colors
deep_water_colors <- data.frame(
  ID = c(300, 400, 500, 600, 700, 800, 900, 1000, 1100),
  color = c(
    "#FFA500",  # B3. Bathyal Slope (Orange)
    "#FF8C00",  # B4. Abyssal Slope (Dark Orange)
    "#FFD700",  # C3. Bathyal Rise (Gold)
    "#B8860B",  # C4. Abyssal Rise (Dark Goldenrod)
    "#4169E1",  # D3. Bathyal Basin (Royal Blue)
    "#000050",  # D4. Abyssal Basin (Navy)
    "#FF4500",  # E1. Euphotic Seamount (Orange Red)
    "#8B0000",  # E2. Mesophotic Seamount (Dark Red)
    "#4A4000"   # E3. Bathyal Seamount (Indigo)
  )
)

# The full theoretical color table
final_color_table <- bind_rows(habitat_colors, deep_water_colors) %>%
  arrange(ID)

# Create the master reference table containing both Names and Colors
master_rat <- habitat_table %>%
  left_join(final_color_table, by = "ID")

# 5. Filtering and Applying Levels ----
# Convert to factor first so terra builds its internal index of existing values
habitats <- as.factor(habitats)

# Extract the active Raster Attribute Table (RAT) that terra just created based on pixels that actually exist
active_rat <- levels(habitats)[[1]]
colnames(active_rat)[1] <- "ID"

# Merge our Master Table ONLY onto the true IDs that exist in the map
updated_rat <- active_rat %>%
  left_join(master_rat, by = "ID") %>%
  as.data.frame()

# Re-apply the perfectly matched names and colors
levels(habitats) <- updated_rat[, c("ID", "Habitat_Name")]
coltab(habitats) <- updated_rat[, c("ID", "color")]
freq(habitats)
# 6. Visualization ----
plot(habitats)
jpeg(filename = "figures/l3_benthic_habitats_v1.jpg", 
     width = 40,        # Width of the image in cm
     height = 50,       # Height of the image in cm
     units = "cm",  
     res = 300)         # Resolution in pixels per inch

# Because coltab and levels are set, terra automatically handles the legend gracefully
plot(habitats, 
     main = "Benthic Habitats (L3) - Brazilian Margin",
     mar = c(3, 3, 3, 25),
     cex.main = 3,
     plg = list(cex = 1.5)) # Reduced cex slightly to accommodate 17 potential classes

dev.off()

# 7. Safe Exporting (8-bit GeoTIFF with Embedded Color Table) ----

# Step A: Convert the factor raster back to raw numbers (stripping the categorical metadata)
hab_numeric <- as.numeric(habitats)

# Step B: Create a strict sequential reclassification matrix
# We assign a sequential ID (1 to N) to the existing classes
updated_rat$Export_ID <- 1:nrow(updated_rat)

# Column 1 = Old high IDs (101, 300, 1100...) | Column 2 = New sequential 8-bit IDs (1, 2, 3...)
rcl_matrix <- as.matrix(updated_rat[, c("ID", "Export_ID")])

# Step C: Reclassify! 
# The 'others = NA' argument guarantees that any stray pixel is converted to NA, preventing INT1U warning.
habitats_export <- classify(hab_numeric, rcl_matrix, others = NA)

# Step D: Convert back to a categorical factor and apply the NEW metadata
habitats_export <- as.factor(habitats_export)

# Apply the names and colors mapped to the new sequential Export_ID
levels(habitats_export) <- updated_rat[, c("Export_ID", "Habitat_Name")]
coltab(habitats_export) <- updated_rat[, c("Export_ID", "color")]

# Step E: Final Export
# With absolutely no values above 255 remaining, INT1U executes perfectly.
writeRaster(habitats_export, 
            'output_data/l3_benthic_habitats_v1.tif', 
            datatype = "INT1U", 
            overwrite = TRUE)

# Export the legend for QGIS or reporting
write.csv(updated_rat, 'output_data/l3_benthic_habitats_legend_v1.csv', row.names = FALSE)
