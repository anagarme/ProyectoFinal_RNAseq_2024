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

## --------------------- Modelo estadístico  ------------------
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


## Para visualizar los resultados estadísticos
plotMA(eb_results, coef = 2)
volcanoplot(eb_results, coef = 4, highlight = 5, names = de_results$gene_name, col = "cornflowerblue", hl.col="darkorange")

de_results[de_results$gene_name %in% c("Xist", "Foxp3", "Klrk1","Ctla4","Flicr"), ]
