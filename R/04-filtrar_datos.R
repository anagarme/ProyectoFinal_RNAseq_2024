## --------------------- FILTRAR LOS DATOS --------------------
## ------------------------------------------------------------

## Proporción de lecturas asignadas a genes
rse_gene_SRP192782$assigned_gene_prop <- rse_gene_SRP192782$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP192782$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP192782$assigned_gene_prop)

## Resumen por grupo de tejido
with(colData(rse_gene_SRP192782), tapply(assigned_gene_prop, sra_attribute.tissue, summary))

## Eliminar muestras de baja calidad y genes con niveles de expresión muy bajos
rse_gene_SRP192782_unfiltered <- rse_gene_SRP192782
## Histograma
hist(rse_gene_SRP192782$assigned_gene_prop, col="plum2")
abline(v=0.4,col="steelblue", lwd=7, lty = "dashed")
## Eliminar las muestras con una proporción baja
table(rse_gene_SRP192782$assigned_gene_prop < 0.4)
## Eliminar los genes con proporciones bajas con base al histograma
rse_gene_SRP192782 <- rse_gene_SRP192782[, rse_gene_SRP192782$assigned_gene_prop > 0.4]

## Eliminar genes con niveles bajos de expresión
gene_means <- rowMeans(assay(rse_gene_SRP192782, "counts"))
summary(gene_means)
rse_gene_SRP192782 <- rse_gene_SRP192782[gene_means > 0.01, ]
dim(rse_gene_SRP192782)

## Porcentaje de genes que retuvimos
round(nrow(rse_gene_SRP192782) / nrow(rse_gene_SRP192782_unfiltered) * 100, 2)
