# Benthic Habitat Mapping of the Brazilian Margin

An automated, reproducible R-based geospatial pipeline for classifying, integrating, and visualizing benthic habitats along the Brazilian continental margin. 

Built to support historical oceanographic data recovery (**Revimar** project) and provide foundational ecological layers for Marine Spatial Planning (**PEM Sul** and broader Brazilian MSP initiatives), this repository processes geomorphological, bathymetric, and sedimentological data to generate publication-ready habitat maps and statistics.

## Project Overview

This pipeline translates raw spatial data into a standardized, ecologically meaningful benthic habitat classification system. It utilizes high-resolution raster math and map algebra via the `terra` package to seamlessly cross physical structure (geomorphology and light penetration) with substrate texture (Folk, 1954).

### Key Methodological Steps:
1. **Structural & Photic Framework:** Categorizes the seafloor into 7 primary zones based on depth and light availability:
   * **Shelf (A):** Euphotic and Mesophotic
   * **Slope (B):** Bathyal Continental Slope
   * **Basin (C):** Oceanic Basin
   * **Seamounts (D):** Euphotic, Mesophotic, and Bathyal Seamounts
2. **Sediment Classification:** Normalizes interpolated Sand, Mud, and Gravel fractions and classifies them into a robust, simplified 5-class Folk system (Gravel, Gravelly Sediment, Sand, Mixed Sand-Mud, Mud) to minimize interpolation noise in regional models.
3. **Habitat Integration:** Uses map algebra to cross the structural framework with the sediment classes, generating 23 highly specific habitat classes (e.g., `A1c. Sand Euphotic Shelf`). *Note: Deep-water features (Bathyal Slope, Oceanic Basin, Bathyal Seamount) are retained as purely structural zones due to sampling constraints.*
4. **Cartographic Blending:** Automatically generates customized `.tif` color tables by mathematically blending standard British Geological Survey (BGS) sediment colors with depth-appropriate structural hues (e.g., cyan for shallow, navy for deep).
