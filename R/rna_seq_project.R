## Name: rna_seq_project
## Author : Ana Marisol García Mejía
## Description:
## Last update: 12/02/2024
## Version:
##
##
##
##
## -------------------- Librerías/Paquetes --------------------
## Cargar las librerías
suppressPackageStartupMessages(library(recount3))
suppressPackageStartupMessages(library(sessioninfo))
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
##
## -------------------- Descargar datos --------------------
##
## Descargar los datos de RNA-seq del estudio de interés
rse_gene_SRP192782 <- create_rse_manual(
  project = "SRP192782",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)
## Explorar los metadatos del objeto RSE SRP192782
metadata_SRP192782 <- metadata(rse_gene_SRP192782)
##
## -------------------- Pre procesamiento de los datos  --------------------

## Convertir las cuentas por nucleótido a cuentas por lectura
assay(rse_gene_SRP192782, "counts") <- compute_read_counts(rse_gene_SRP192782)
## Reemplazar guiones cortos que podrían ocasionar problemas en varios atributos
rse_gene_SRP192782$sra.sample_attributes[1:5]
rse_gene_SRP192782$sra.sample_attributes <- gsub("cell-id-assigned", "cell_id_assigned", rse_gene_SRP192782$sra.sample_attributes)
rse_gene_SRP192782$sra.sample_attributes <- gsub("cell-id-sorted", "cell_id_sorted", rse_gene_SRP192782$sra.sample_attributes)
rse_gene_SRP192782$sra.sample_attributes <- gsub("mouse-id", "mouse_id", rse_gene_SRP192782$sra.sample_attributes)
##
## Expandir los atributos SRA para facilitar el manejo de la información
rse_gene_SRP192782 <- expand_sra_attributes(rse_gene_SRP192782)
# Verificar que los atributes no generen algún problema
rse_gene_SRP192782$sra.sample_attributes [1]
colData(rse_gene_SRP192782)[,grepl("^sra_attribute", colnames(colData(rse_gene_SRP192782)))]
##
## Cambiar el tipo de character a numeric o factor
rse_gene_SRP192782$sra_attribute.cell_id_assigned <- as.numeric(rse_gene_SRP192782$sra_attribute.cell_id_assigned)
rse_gene_SRP192782$sra_attribute. <- as.numeric(rse_gene_SRP192782$sra_attribute.)

##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
