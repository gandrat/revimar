## Final Benthic Habitat Integration & Cartography Script
## Excludes sediment from Bathyal Slope, Oceanic Basin, and Bathyal Seamounts

# 1. Loading Packages and Data ----
library(terra)
library(dplyr)
library(tidyr)

geom_photic <- rast('output_data/geomorphology_photic_zones_v2.tif')
sediments <- rast('output_data/folk_classification_5classes.tif') 

plot(geom_photic)

# 2. Map Algebra (The Intersection) ----
# Shelf (1, 2) and shallow Seamounts (5, 6) get sediment IDs added.
# Slope (3), Basin (4), and Bathyal Seamount (7) get multiplied by 100 only.
habitats <- ifel(
  geom_photic %in% c(1, 2, 5, 6), 
  (geom_photic * 100) + sediments, 
  geom_photic * 100
)


# 3. Generating the Attribute Tables (As before) ----
base_zones <- data.frame(
  Zone_ID = c(1, 2, 5, 6),
  Prefix = c("A1", "A2", "D1", "D2"),
  Zone_Name = c("Euphotic Shelf", "Mesophotic Shelf", "Euphotic Seamount", "Mesophotic Seamount")
)

folk_classes <- data.frame(
  Sed_ID = 1:5,
  Letter = letters[1:5],
  Sed_Name = c("Gravel", "Gravelly Sediment", "Sand", "Mixed Sand-Mud", "Mud")
)

combinations <- expand.grid(Zone_ID = c(1, 2, 5, 6), Sed_ID = 1:5) %>%
  left_join(base_zones, by = "Zone_ID") %>%
  left_join(folk_classes, by = "Sed_ID") %>%
  mutate(
    ID = (Zone_ID * 100) + Sed_ID,
    Habitat_Name = paste0(Prefix, Letter, ". ", Sed_Name, " ", Zone_Name)
  ) %>%
  select(ID, Habitat_Name)

zones_without_sediment <- data.frame(
  ID = c(300, 400, 700),
  Habitat_Name = c("B3. Bathial Continental Slope", "C4. Oceanic Basin", "D3. Bathial Seamount")
)

# The full theoretical table
habitat_table <- bind_rows(combinations, zones_without_sediment) %>%
  arrange(Habitat_Name)


# 4. Generating the Color Palette (As before) ----
geo_colors <- data.frame(Zone_ID = c(1, 2, 5, 6), Geo_Hex = c("#00FFFF", "#00688B", "#FF4500", "#8B0000"))
sed_colors <- data.frame(Sed_ID = 1:5, Sed_Hex = c("#B95246", "#C68F8A", "#FCE47F", "#9EB25C", "#4A7A40"))

mix_colors <- function(c1, c2, weight = 0.6) {
  rgb1 <- col2rgb(c1)
  rgb2 <- col2rgb(c2)
  mixed <- round(rgb1 * (1 - weight) + rgb2 * weight)
  rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

habitat_colors <- expand.grid(Zone_ID = c(1, 2, 5, 6), Sed_ID = 1:5) %>%
  left_join(geo_colors, by = "Zone_ID") %>%
  left_join(sed_colors, by = "Sed_ID") %>%
  rowwise() %>%
  mutate(value = (Zone_ID * 100) + Sed_ID, color = mix_colors(Geo_Hex, Sed_Hex, weight = 0.6)) %>%
  ungroup() %>%
  select(value, color)

deep_water_colors <- data.frame(
  value = c(300, 400, 700),
  color = c("#FFA500", "#000050", "#4A4000")
)

# The full theoretical color table
final_color_table <- bind_rows(habitat_colors, deep_water_colors) %>%
  arrange(value)


# 5. THE FIX: Filtering and Applying Levels ----
# Find out exactly which IDs actually exist in the raster using freq()
existing_ids <- freq(habitats)$value

# Filter our theoretical tables to ONLY include classes that actually occur
habitat_table_actual <- habitat_table %>% 
  filter(ID %in% existing_ids)

final_color_table_actual <- final_color_table %>% 
  filter(value %in% existing_ids) %>%
  as.data.frame() # Keeping the tibble-to-dataframe fix for the colors

# Ensure the raster is categorical and apply the perfectly matched tables
habitats <- as.factor(habitats)

levels(habitats) <- habitat_table_actual
coltab(habitats) <- final_color_table_actual

# Visualization
jpeg(filename = "figures/benthic_habitats_l3.jpg", 
     width = 18,       # Width of the image
     height = 17,       # Height of the image
     units = "cm",     # Units for width/height (inches)
     res = 300)        # Resolution in pixels per inch

plot(habitats, main = "Benthic Habitats (L3)",
     mar = c(3, 3, 3, 17))
dev.off()
# Exporting the true, filtered maps and legends
writeRaster(habitats, 'output_data/benthic_habitats_v1.tif', overwrite = TRUE)
write.csv(habitat_table_actual, 'output_data/benthic_habitats_legend.csv', row.names = FALSE)