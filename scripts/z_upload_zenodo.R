library(zen4R)

# 1. Authentication (Generate this token in your Zenodo dashboard)
mytoken <- Sys.getenv("ZENODO_TOKEN")

# Secure authentication
zenodo <- ZenodoManager$new(token = mytoken)

# 2. Create a new deposit (Draft)
my_deposit <- zenodo$createEmptyRecord()

# 3. Automated Basic Metadata completion (DataCite)
my_deposit$setTitle("REVIMAR-MAP: Benthic Marine Habitats Mapping of Brazil")
my_deposit$setDescription("Continuous spatial modeling integrating geomorphological, photic, and sedimentological variables.")
my_deposit$addCreator(name = "Gandra, Tiago", affiliation = "Federal Institute of Rio Grande do Sul")
my_deposit$setResourceType("dataset")

# 4. Attaching the final map exported by the 'terra' package
zenodo$uploadFile("outputs/l1_benthic_provinces.tif", record = my_deposit)
zenodo$uploadFile("outputs/l2_biological_zones", record = my_deposit)
zenodo$uploadFile("outputs/l3_substrate.tif", record = my_deposit)

# 5. Automatic publication (Generates the DOI for the data)
zenodo$publishRecord(my_deposit$id)