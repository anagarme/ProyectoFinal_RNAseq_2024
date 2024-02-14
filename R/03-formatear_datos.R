## ---------------- PROCESAMIENTO DE LOS DATOS ----------------
## --------------------- FORMATEAR ----------------------------

## Convertir las cuentas por nucleótido a cuentas por lectura
assay(rse_gene_SRP192782, "counts") <- compute_read_counts(rse_gene_SRP192782)

## Mostrar los atributos del objeto
rse_gene_SRP192782$sra.sample_attributes[1:3]

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
rse_gene_SRP192782$sra_attribute.tissue <- factor(rse_gene_SRP192782$sra_attribute.tissue)
rse_gene_SRP192782$sra_attribute.mouse_id <- factor(rse_gene_SRP192782$sra_attribute.mouse_id)
## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP192782)[ ,grepl("^sra_attribute.[cell_id_assigned|cell_id_sorted|tissue]", colnames(colData(rse_gene_SRP192782)))]))
