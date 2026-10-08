## Benthic Habitat Classification Script (EUNIS modified Folk + Biogenic)
## Implements the EUNIS Folk sediment classification (Long, 2006)

packages <- c('terra', 'sf', 'dplyr')

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# 1. Loading Environment and Data ----
sed <- rast('data/processed/sediment.tif')
plot(sed)

sand   <- sed$sand
mud    <- sed$mud
gravel <- sed$gravel

# 2. Raster Math Normalization (100% Assurance) ----
sand   <- clamp(sand, lower = 0)
mud    <- clamp(mud, lower = 0)
gravel <- clamp(gravel, lower = 0)

total_sum <- sand + mud + gravel
plot(total_sum, main = "Total Sum Before Normalization")

sand_norm   <- (sand / total_sum) * 100
mud_norm    <- (mud / total_sum) * 100
gravel_norm <- (gravel / total_sum) * 100

sed_norm <- c(sand_norm, mud_norm, gravel_norm)
writeRaster(sed_norm, 'data/processed/sediment_norm.tif', overwrite = TRUE)

# 3. Calculate Sand/Mud Ratio (Diagram Axis) ----
sum_sand_mud <- sand_norm + mud_norm
sand_ratio   <- ifel(sum_sand_mud == 0, 0, (sand_norm / sum_sand_mud) * 100)
plot(sand_ratio, main = "Sand Ratio (%)")

# 4. Apply EUNIS Folk Boundaries (4 Classes) ----
folk_classes <- 
  # 1: Coarse sediment (Gravel >= 80% OU gravel >= 5% com Sand >= 90%)
  ifel(gravel_norm >= 80 | (gravel_norm >= 5 & sand_ratio >= 90), 1,
       
       # 2: Mixed sediment (Cascalho 5% a 80% com Sand < 90%)
       ifel(gravel_norm >= 5 & sand_ratio < 90, 2,
            
            # 3: Sand (gravel < 5% com Sand >= 90%)
            ifel(gravel_norm < 5 & sand_ratio >= 90, 3,
                 
                 # 4: Mud to muddy sand (gravel < 5% com Sand < 90%)
                 ifel(gravel_norm < 5 & sand_ratio < 90, 4, NA))))

plot(folk_classes)
# 5. Add Biogenic Substrate (Class 5) ----
bio <- read_sf('gis/gis_revimar.gpkg', layer = 'biogenico')
bio_vect <- project(vect(bio), folk_classes)
bio_rast <- rasterize(bio_vect, folk_classes, touches = TRUE)

# Substitui pela classe 5 (Biogenic) onde houver polígono biogênico
sediment_classes <- ifel(!is.na(bio_rast), 5, folk_classes)

# 6. Convert to Categorical Raster (Metadata RAT) ----
sediment_classes <- as.factor(sediment_classes)

folk_table <- data.frame(
  ID = 1:5,
  Habitat = c(
    "Coarse sediment", 
    "Mixed sediment", 
    "Sand", 
    "Mud to muddy sand", 
    "Biogenic"
  )
)

levels(sediment_classes) <- folk_table

# 7. Defining Color Palette (EUNIS Diagram Palette + Hot Pink) ----
folk_color_table <- data.frame(
  value = 1:5,
  color = c(
    "#B95246",  # 1: Coarse sediment (Terracotta / Reddish Brown)
    "#C68F8A",  # 2: Mixed sediment (Mauve / Dusty Pink)
    "#F6D173",  # 3: Sand (Yellow)
    "#7F9953",  # 4: Mud to muddy sand (Olive Green)
    "#FF69B4"   # 5: Biogenic (Hot Pink)
  )
)

coltab(sediment_classes) <- folk_color_table

# 8. Visual Export ----
jpeg(filename = "figures/sediment_eunis_biogenic.jpg", 
     width = 18, 
     height = 17, 
     units = "cm", 
     res = 300)
plot(sediment_classes, main = "EUNIS Folk & Biogenic Substrates")
dev.off()

plot(sediment_classes, main = "EUNIS Folk & Biogenic Substrates")

# 9. Export Final GeoTIFF ----
writeRaster(
  sediment_classes, 
  'data/processed/sediment_folk_classification_5classes.tif', 
  datatype = "INT1U",
  overwrite = TRUE
)
