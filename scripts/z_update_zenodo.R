
library(zen4R)

# Retrieve Zenodo token and initialize manager
zenodo_token <- Sys.getenv("ZENODO_TOKEN")
zenodo <- ZenodoManager$new(token = zenodo_token)

# Retrieve the official record (e.g., creating a new version)
record_id <- "22834505"
my_deposit <- zenodo$depositRecordVersion(record_id)

# Update metadata for the new release
my_deposit$setVersion("0.5")
my_deposit$setLicense("cc-by-4.0")
my_deposit$setKeywords(c("R", "Marine Habitats", "Benthic", "Active Development"))

# 4. Attaching the final map exported by the 'terra' package
zenodo$uploadFile("outputs/l1_benthic_provinces.tif", record = my_deposit)
zenodo$uploadFile("outputs/l2_biological_zones", record = my_deposit)
zenodo$uploadFile("outputs/l3_substrate.tif", record = my_deposit)

# Synchronize metadata with Zenodo (Publishing step requires explicit confirmation)
my_deposit <- zenodo$depositRecord(my_deposit)
# zenodo$publishRecord(my_deposit$id)