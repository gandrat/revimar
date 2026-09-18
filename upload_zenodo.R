
library(zen4R)

# 1. Autenticação (Você gera este token no painel do Zenodo)
mytoken <- Sys.getenv("ZENODO_TOKEN")

# Autenticação segura
zenodo <- ZenodoManager$new(token = mytoken)

# 2. Criação de um novo depósito (Draft)
meu_deposito <- zenodo$createEmptyRecord()

# 3. Preenchimento automatizado de Metadados Básicos (DataCite)
meu_deposito$setTitle("REVIMAR-MAP: Mapeamento de Habitats Marinhos Bentônicos do Brasil")
meu_deposito$setDescription("Modelagem espacial contínua integrando variáveis geomorfológicas, fóticas e sedimentológicas.")
meu_deposito$addCreator(name = "Gandra, Tiago", affiliation = "Federal Institute of Rio Grande do Sul")
meu_deposito$setResourceType("dataset")

# 4. Anexando o mapa final exportado pelo pacote 'terra'
zenodo$uploadFile("output_data/l2_province_photic_zones.tif", record = meu_deposito)
zenodo$uploadFile("output_data/l1_benthic_provinces_v2.tif", record = meu_deposito)

# 5. Publicação automática (Gera o DOI para os dados)
zenodo$publishRecord(meu_deposito$id)