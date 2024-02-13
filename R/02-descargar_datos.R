## -------------------- DESCARGAR LOS DATOS --------------------
## -------------------------------------------------------------

## Descargar los datos de RNA-seq del estudio de interés
rse_gene_SRP192782 <- create_rse_manual(
  project = "SRP192782",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)
## Explorar el objeto RSE
rse_gene_SRP192782
## Explorar los metadatos del objeto RSE SRP192782
metadata(rse_gene_SRP192782)
