## Final Benthic Habitat Integration & Cartography Script (L3)
## Integrates EUNIS 4-class sediment model + Biogenic layer with L2 biological zones

library(terra)
library(dplyr)
library(scales)

crs_albers_brasil <- "+proj=aea +lat_0=-12 +lon_0=-54 +lat_1=-2 +lat_2=-22 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# 1. Load Data and Align Spatial Grids ----
geom_photic <- rast('outputs/l2_biological_zones.tif')
sediments   <- rast('data/processed/sediment_folk_classification_5classes.tif') 
sediments   <- project(sediments, geom_photic, method = 'near')

# 2. Map Algebra (L2 Zones x Substrate Intersection) ----
# Shelf (1, 2) and shallow Seamounts (7, 8) receive sediment IDs (1 to 5).
# Deep-water features (3, 4, 5, 6, 9) remain multiplied by 100 without sediment overlay.
habitats_raw <- ifel(
  geom_photic %in% c(1, 2, 7, 8), 
  (geom_photic * 100) + sediments, 
  geom_photic * 100
)

# 3. Construct Master Reference Table (Habitat Names & Color Scheme) ----
# Structural photic zones that receive sediment data
base_zones <- data.frame(
  Zone_ID   = c(1, 2, 7, 8),
  Prefix    = c("A1", "A2", "E1", "E2"),
  Zone_Name = c("Euphotic Shelf", "Mesophotic Shelf", "Euphotic Seamount", "Mesophotic Seamount"),
  Geo_Hex   = c("#00FFFF", "#00688B", "#FF4500", "#8B0000")
)

# EUNIS modified sediment classes + Biogenic
folk_classes <- data.frame(
  Sed_ID   = 1:5,
  Letter   = letters[1:5],
  Sed_Name = c("Coarse", "Mixed", "Sand", "Mud", "Biogenic"),
  Sed_Hex  = c("#B95246", "#C68F8A", "#F6D173", "#7F9953", "#FF69B4")
)

# Mathematical color blending function
mix_colors <- function(c1, c2, weight = 0.6) {
  rgb1 <- col2rgb(c1)
  rgb2 <- col2rgb(c2)
  mixed <- round(rgb1 * (1 - weight) + rgb2 * weight)
  rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

# Shallow combinations with substrate (IDs: 101 to 805)
shallow_table <- expand.grid(Zone_ID = c(1, 2, 7, 8), Sed_ID = 1:5) %>%
  left_join(base_zones, by = "Zone_ID") %>%
  left_join(folk_classes, by = "Sed_ID") %>%
  rowwise() %>%
  mutate(
    Raw_ID       = (Zone_ID * 100) + Sed_ID,
    Habitat_Name = paste0(Prefix, Letter, ". ", Sed_Name, " ", Zone_Name),
    color        = mix_colors(Geo_Hex, Sed_Hex, weight = 0.6)
  ) %>%
  ungroup() %>%
  select(Raw_ID, Habitat_Name, color)

# Deep-water features without substrate data (IDs: 300 to 900)
deep_table <- data.frame(
  Raw_ID       = c(300, 400, 500, 600, 900),
  Habitat_Name = c("B3. Bathyal Slope",
                   "C3. Bathyal Continental Rise",
                   "C4. Abyssal Continental Rise",
                   "D4. Abyssal Oceanic Basin",
                   "E3. Bathyal Seamount"),
  color        = c("#FFA500", "#FFD700", "#B8860B", "#000050", "#8C3100")
)

# Unified theoretical master table
master_table <- bind_rows(shallow_table, deep_table) %>%
  arrange(Raw_ID)

# 4. Filter Active Map IDs and Map to 8-bit Sequential IDs (1 to N) ----
present_ids <- freq(habitats_raw)$value

final_rat <- master_table %>%
  filter(Raw_ID %in% present_ids) %>%
  mutate(ID_8bit = row_number())

# 5. Direct 8-bit Reclassification and Metadata Attachment ----
# Reclassification matrix: Column 1 = Raw_ID, Column 2 = ID_8bit
rcl_mat <- as.matrix(final_rat[, c("Raw_ID", "ID_8bit")])
habitats_8bit <- classify(habitats_raw, rcl_mat, others = NA)

# Convert to categorical SpatRaster and attach RAT and Color Palette
habitats_8bit <- as.factor(habitats_8bit)
levels(habitats_8bit) <- data.frame(ID = final_rat$ID_8bit, Habitat = final_rat$Habitat_Name)
coltab(habitats_8bit) <- data.frame(value = final_rat$ID_8bit, color = final_rat$color)
plot(habitats_8bit)
# 6. Data Export and Cartographic Output ----
# Export native 8-bit GeoTIFF with embedded color palette
writeRaster(
  habitats_8bit, 
  'outputs/l3_benthic_habitats_substrate.tif', 
  datatype  = "INT1U", 
  overwrite = TRUE
)


# Generate final cartographic visualization
jpeg(filename = "figures/l3_benthic_habitats.jpg", 
     width = 40, height = 45, units = "cm", res = 300)

plot(habitats_8bit, 
     main     = "Benthic Habitats of the Brazilian Margin (L3)",
     col      = final_rat$color,
     mar      = c(3, 3, 3, 26),
     cex.main = 2,
     plg      = list(cex = 1.3))

dev.off()
