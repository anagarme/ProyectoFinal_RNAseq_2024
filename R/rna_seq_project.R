## Name: rna_seq_project
## Author : Ana Marisol García Mejía
## Last update: 12/02/2024
## Version: 2.0



## -------------------- PAQUETES ------------------------------
## ------------------------------------------------------------
## Cargar los paquetes necesarios

suppressPackageStartupMessages(library(recount3))
suppressPackageStartupMessages(library(sessioninfo))
suppressPackageStartupMessages(library(edgeR))
library(limma)
library(ggplot2)
library(pheatmap)

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

## ---------------- PROCESAMIENTO DE LOS DATOS ----------------
## ------------------------------------------------------------

## FORMATEAR LOS DATOS
## Convertir las cuentas por nucleótido a cuentas por lectura
assay(rse_gene_SRP192782, "counts") <- compute_read_counts(rse_gene_SRP192782)

## Mostrar los atributos del objeto
rse_gene_SRP192782$sra.sample_attributes[1:5]

## Reemplazar guiones cortos que podrían ocasionar problemas en varios atributos
rse_gene_SRP192782$sra.sample_attributes <- gsub("cell-id-assigned", "cell_id_assigned", rse_gene_SRP192782$sra.sample_attributes)
rse_gene_SRP192782$sra.sample_attributes <- gsub("cell-id-sorted", "cell_id_sorted", rse_gene_SRP192782$sra.sample_attributes)
rse_gene_SRP192782$sra.sample_attributes <- gsub("mouse-id", "mouse_id", rse_gene_SRP192782$sra.sample_attributes)


## Verificar que los atributos no generen algún problema
rse_gene_SRP192782$sra.sample_attributes [1]
## Hacer la información más fácil de usar
rse_gene_SRP192782 <- expand_sra_attributes(rse_gene_SRP192782)
colData(rse_gene_SRP192782)[,grepl("^sra_attribute", colnames(colData(rse_gene_SRP192782)))]

## Cambiar el tipo de character a factor para los análisis estadísticos
rse_gene_SRP192782$sra_attribute.cell_id_assigned <- factor(rse_gene_SRP192782$sra_attribute.cell_id_assigned)
rse_gene_SRP192782$sra_attribute.cell_id_sorted <- factor(rse_gene_SRP192782$sra_attribute.cell_id_sorted)
rse_gene_SRP192782$sra_attribute.source_name <- factor(rse_gene_SRP192782$sra_attribute.source_name)
rse_gene_SRP192782$sra_attribute.tissue <- factor(rse_gene_SRP192782$sra_attribute.tissue)
rse_gene_SRP192782$sra_attribute.mouse_id <- factor(rse_gene_SRP192782$sra_attribute.mouse_id)
## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP192782)[ ,grepl("^sra_attribute.[cell_id_assigned|cell_id_sorted|tissue]", colnames(colData(rse_gene_SRP192782)))]))


## --------------------- FILTRAR LOS DATOS --------------------
## ------------------------------------------------------------

## Proporción de lecturas asignadas a genes
rse_gene_SRP192782$assigned_gene_prop <- rse_gene_SRP192782$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP192782$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP192782$assigned_gene_prop)

## Resumen por grupo de tejido
with(colData(rse_gene_SRP192782), tapply(assigned_gene_prop, sra_attribute.tissue, summary))

## Eliminar muestras de baja calidad y genes con niveles de expresión muy bajos
## Copia de seguridad
rse_gene_SRP192782_unfiltered <- rse_gene_SRP192782

## Histograma
hist(rse_gene_SRP192782$assigned_gene_prop, col="plum2")
abline(v=0.35,col="steelblue", lwd=7, lty = "dashed")
## Eliminar las muestras con una proporción baja
table(rse_gene_SRP192782$assigned_gene_prop < 0.35)
## Eliminar los genes con proporciones baja
rse_gene_SRP192782 <- rse_gene_SRP192782[, rse_gene_SRP192782$assigned_gene_prop > 0.35]

## Eliminar genes con niveles bajos de expresión
gene_means <- rowMeans(assay(rse_gene_SRP192782, "counts"))
summary(gene_means)
rse_gene_SRP192782 <- rse_gene_SRP192782[gene_means > 0.01, ]
dim(rse_gene_SRP192782)

## Porcentaje de genes que retuvimos
round(nrow(rse_gene_SRP192782) / nrow(rse_gene_SRP192782_unfiltered) * 100, 2)

## --------------------- NORMALIZACIÓN  --------------------
## ------------------------------------------------------------
dge <- DGEList(
  counts = assay(rse_gene_SRP192782, "counts"),
  genes = rowData(rse_gene_SRP192782)
)
dge <- calcNormFactors(dge)

## --------------------- EXPRESIÓN DIFERENCIAL  ---------------
## ------------------------------------------------------------

## Tejido
ggplot(as.data.frame(colData(rse_gene_SRP192782)), aes(y = assigned_gene_prop, x = sra_attribute.tissue)) +
  geom_boxplot() +
  theme_classic(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Tissue group")

## Células T
ggplot(as.data.frame(colData(rse_gene_SRP192782)), aes(y = assigned_gene_prop, x = sra_attribute.cell_id_sorted)) +
  geom_boxplot() +
  theme_classic(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("T cells group")

## --------------------- Modelo estadístico  ---------------
## ------------------------------------------------------------
mod<- model.matrix(~ 0 + sra_attribute.tissue+ sra_attribute.cell_id_sorted + assigned_gene_prop,
                    data = colData(rse_gene_SRP192782)
)
colnames(mod)

vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP192782),
  sort.by = "none"
)
dim(de_results)


## Visualicemos los resultados estadísticos
plotMA(eb_results, coef = 2)

volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)
volcanoplot(eb_results, coef = 4, highlight = 5, names = de_results$gene_name, col = "cornflowerblue", hl.col="darkorange")

de_results[de_results$gene_name %in% c("Xist", "Foxp3", "Klrk1","Ctla4","Flicr"), ]
## --------------------- VISUALIZACIÓN  --------------------
## ------------------------------------------------------------


## Creemos una tabla con información de las muestras
## y con nombres de columnas más amigables

exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 30, ]
## Crear un data frame con la información de las muestras
df <- as.data.frame(colData(rse_gene_SRP192782)[, c("sra_attribute.cell_id_assigned", "sra_attribute.tissue", "sra_attribute.cell_id_sorted")])
colnames(df) <- c("Cell id", "Tissue", "T Cells")
## Cambiar el nombre de los identificadores
nombres <- rownames(de_results)
rownames(exprs_heatmap) <- de_results$gene_name[match(rownames(exprs_heatmap), nombres)]

## Heatmap

pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)

