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
