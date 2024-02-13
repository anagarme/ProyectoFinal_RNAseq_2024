## --------------------- NORMALIZACIÓN  --------------------
## ------------------------------------------------------------
dge <- DGEList(
  counts = assay(rse_gene_SRP192782, "counts"),
  genes = rowData(rse_gene_SRP192782)
)
dge <- calcNormFactors(dge)
