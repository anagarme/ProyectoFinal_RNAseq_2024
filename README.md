# Proyecto final RNA-seq 2024-2
_Ana Marisol García Mejía_

Actualización 12/02/2024
## Sobre el proyecto
El objetivo de este proyecto es aplicar los conocimientos y habilidades adquiridas durante el módulo "Análisis de datos de secuenciación masiva"
impartido en el curso de Bioinformática y Estadística II en el semestre 2024-2. 
### Estudio de interés: SRP192782
Los datos utilizados para el análisis fueron obtenidos de [IL-33 Signaling Alters Regulatory T Cell Diversity in Support of Tumor Development](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6990979/).
**Abstract:** 
"Regulatory T cells (Tregs) can impair anti-tumor immune responses and are associated with poor prognosis in multiple cancer types. 
Tregs in human tumors span diverse transcriptional states distinct from those of peripheral Tregs, but their contribution to tumor
development remains unknown. Here, we use single-cell RNA sequencing (RNA-seq) to longitudinally profile dynamic shifts in the 
distribution of Tregs in a genetically engineered mouse model of lung adenocarcinoma. In this model, interferon-responsive Tregs
are more prevalent early in tumor development, whereas a specialized effector phenotype characterized by enhanced expression of 
the interleukin-33 receptor ST2 is predominant in advanced disease. Treg-specific deletion of ST2 alters the evolution of effector
Treg diversity, increases infiltration of CD8+ T cells into tumors, and decreases tumor burden. Our study shows that ST2 plays a 
critical role in Treg-mediated immunosuppression in cancer, highlighting potential paths for therapeutic intervention."
## Metodología 
Para realizar el siguiente análisis es fundamental correr los sripts contenidos en la carpeta [R](https://github.com/anagarme/ProyectoFinal_RNAseq_2024/tree/01f15f3fa0ef73dc2723b17332ab8acb0410caea/R) de manera secuencial.
### Librerías
En el primero de los scripts, [R/01-librerias.R](https://github.com/anagarme/ProyectoFinal_RNAseq_2024/blob/01f15f3fa0ef73dc2723b17332ab8acb0410caea/R/01-librerias.R), se encuentran las librerías necesarias para la descarga, manipulación, análisis y visualización de los datos.
### Descargar datos
Posteriormente, el script [02-descargar_datos.R](https://github.com/anagarme/ProyectoFinal_RNAseq_2024/blob/01f15f3fa0ef73dc2723b17332ab8acb0410caea/R/02-descargar_datos.R) descarga los datos de RNA-seq del estudio
de interés a través del paquete 'recount3'. Los cuales se guardan en el objeto _rse_gene_SRP192782_. 
```
rse_gene_SRP192782 <- create_rse_manual(
  project = "SRP192782",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)
```
### Formateo de los datos
El script [03-formatear_datos.R](https://github.com/anagarme/ProyectoFinal_RNAseq_2024/blob/01f15f3fa0ef73dc2723b17332ab8acb0410caea/R/03-formatear_datos.R) se encarga de procesar los datos
para que su manejo sea más sencillo, y de ser necesario hacer correcciones para que no existan complicaciones.

Primero convierte las cuentas sin procesar en recuentos de lectura.
```
## Convertir las cuentas por nucleótido a cuentas por lectura
assay(rse_gene_SRP192782, "counts") <- compute_read_counts(rse_gene_SRP192782)
```
Adicionalmente exploramos los atributos del objeto RSE.
```
rse_gene_SRP192782$sra.sample_attributes[1:3]
  [1] "cell-id-assigned;;Tconv|cell-id-sorted;;Tconv|mouse-id;;Young2|source_name;;mouse tissue|tissue;;Lung"
  [2] "cell-id-assigned;;Treg|cell-id-sorted;;Treg|mouse-id;;4400|source_name;;mouse tissue|tissue;;Lung"    
  [3] "cell-id-assigned;;Tconv|cell-id-sorted;;Tconv|mouse-id;;3642|source_name;;mouse tissue|tissue;;Lung"  
```
Sin embargo, podemos observar que algunos nombres de los atributos tienen guiones, lo cual puede traer consigo algunos problemas para el análisis. Para ello se modifican dichos nombres,
sustituyendo los guiones por guiones bajos.

```
rse_gene_SRP192782$sra.sample_attributes <- gsub("cell-id-assigned", "cell_id_assigned", rse_gene_SRP192782$sra.sample_attributes)
rse_gene_SRP192782$sra.sample_attributes <- gsub("cell-id-sorted", "cell_id_sorted", rse_gene_SRP192782$sra.sample_attributes)
rse_gene_SRP192782$sra.sample_attributes <- gsub("mouse-id", "mouse_id", rse_gene_SRP192782$sra.sample_attributes)
```
Por otro lado, con la función `expand_sra_attributes()` se expanden y separa los atributos, los cuales se incorporan a un _data frame_.
El script también sirve para cambiar el tipo de dato, de _character_ a _factor_, con la finalidad de facilitar los análisis estadísticos.
```
rse_gene_SRP192782 <- expand_sra_attributes(rse_gene_SRP192782)
colData(rse_gene_SRP192782)[,grepl("^sra_attribute", colnames(colData(rse_gene_SRP192782)))]

## Cambiar el tipo de character a factor para los análisis estadísticos
rse_gene_SRP192782$sra_attribute.cell_id_assigned <- factor(rse_gene_SRP192782$sra_attribute.cell_id_assigned)
rse_gene_SRP192782$sra_attribute.cell_id_sorted <- factor(rse_gene_SRP192782$sra_attribute.cell_id_sorted)
rse_gene_SRP192782$sra_attribute.source_name <- factor(rse_gene_SRP192782$sra_attribute.source_name)
rse_gene_SRP192782$sra_attribute.tissue <- factor(rse_gene_SRP192782$sra_attribute.tissue)
rse_gene_SRP192782$sra_attribute.mouse_id <- factor(rse_gene_SRP192782$sra_attribute.mouse_id)
```
### Filtrado de datos
Una vez que los datos tienen el tipo de formato adecuado, se procede a descartar aquellos datos que no tengan una buena calidad o aquellos que su nivel de expresión no se significativo, 
esto a través del script [04-filtrar_datos.R](https://github.com/anagarme/ProyectoFinal_RNAseq_2024/blob/86754e6acade9922268ecfe9e494c57669c3396f/R/04-filtrar_datos.R). 

Se calcula la proporción de lecturas asignadas a genes al dividir las cuentas de genes asignados por las cuentas totales.
```
rse_gene_SRP192782$assigned_gene_prop <- rse_gene_SRP192782$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP192782$recount_qc.gene_fc_count_all.total
```
A continuación se eliminan algunas muestras que se considerende baja calidad y genes con niveles de expresión muy bajos. Esto mediante la visualización de gráficos. 

```
## Copia de seguridad
rse_gene_SRP192782_unfiltered <- rse_gene_SRP192782
## Histograma
hist(rse_gene_SRP192782$assigned_gene_prop, col="plum2")
abline(v=0.35,col="steelblue", lwd=7, lty = "dashed")
```
![](/histogram_assigned_gene_prop.png) 
Al observar dicho histograma nos permite elegir un valor de corte.













![](/results/aa_count_plot_SARS_CoV_2_human_CHE_SARS_CoV_2.png) # Para impag
