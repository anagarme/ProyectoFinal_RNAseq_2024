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
En el primero de los scripts, [01-librerias.R](https://github.com/anagarme/ProyectoFinal_RNAseq_2024/blob/01f15f3fa0ef73dc2723b17332ab8acb0410caea/R/01-librerias.R), se encuentran las librerías necesarias para la descarga, manipulación, análisis y visualización de los datos.
### Descargar datos
Posteriormente, el script [02-descargar_datos.R](https://github.com/anagarme/ProyectoFinal_RNAseq_2024/blob/01f15f3fa0ef73dc2723b17332ab8acb0410caea/R/02-descargar_datos.R) descarga los datos de RNA-seq del estudio
de interés a través del paquete _recount3_. Los cuales se guardan en el objeto _rse_gene_SRP192782_. 
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
Por otro lado, con la función `expand_sra_attributes()` se expanden y separan los atributos, los cuales se incorporan a un _data frame_.
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
![](plots/histogram_assigned_gene_prop.png) 


Al observar dicho histograma nos permite elegir un valor de corte. En este caso 0.4, el cual se utiliza para eliminar las muestras con una proporción menor. 
```
table(rse_gene_SRP192782$assigned_gene_prop < 0.4)
## Eliminar los genes con proporciones bajas con base al histograma
rse_gene_SRP192782 <- rse_gene_SRP192782[, rse_gene_SRP192782$assigned_gene_prop > 0.4]
```
Posteriormene calcula los niveles medios de expresión de los genes en las muestras y elimina los genes con niveles muy bajos.
```
## Eliminar genes con niveles bajos de expresión
gene_means <- rowMeans(assay(rse_gene_SRP192782, "counts"))
summary(gene_means)
rse_gene_SRP192782 <- rse_gene_SRP192782[gene_means > 0.01, ]
```

Por último, se calcula el porcentaje de genes filtrados.
```
round(nrow(rse_gene_SRP192782) / nrow(rse_gene_SRP192782_unfiltered) * 100, 2)
[1] 56.59
```
Lo que significa que se conservó 56.59% de los genes. Es decir, las dimensiones originales del objeto fueron 55,421 genes y 2,942 muestras y ahora se tiene 31,361 genes y 2703 muestras.

### Normalización
Además, con el script [05-normalizacion.R](https://github.com/anagarme/ProyectoFinal_RNAseq_2024/blob/247b79c1e9a1a96188ea51860ac5c32dab675e08/R/05-normalizacion.R) se normalizaron los datos con el paquete `edgeR` con la finalidad de reducir la indicendia de falsos positivos y para que todos tengan la misma escala.

### Expresión diferencial
Adicionalmente, la finalidad del script [06-expresion_diferencial.R](https://github.com/anagarme/ProyectoFinal_RNAseq_2024/blob/d4b71a22c8e630d6003b0350f50cdb8d467ed29d/R/06-expresion_diferencial.R) es la visualización de los datos mediante gráficos comos _boxplots_ para analizar la diferencia en la expresión génica de las muestras en distintas
variables. Observamos que no hay diferencia significativa entre los grupos de tejidos. No obstante, en el grupo de células T reguladoras (_Treg_) se nota que la expresión de genes es 
mayor que en los otros tipos de células T.
```
## Boxplot por grupo de tejido
ggplot(as.data.frame(colData(rse_gene_SRP192782)), aes(y = assigned_gene_prop, x = sra_attribute.tissue)) +
  geom_boxplot() +
  theme_classic(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Tissue group")
```
![](plots/BoxplotTissue.png) 
```
## Boxplot por grupo de células T
ggplot(as.data.frame(colData(rse_gene_SRP192782)), aes(y = assigned_gene_prop, x = sra_attribute.cell_id_sorted)) +
  geom_boxplot() +
  theme_classic(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("T cells group")
```
![](plots/BoxplotTcells.png) 

Asimismo, el script permite crear un modelo estadístico con base en las variables antes mencionadas. Dicho modelo convierte las cuentas de las lecturas por millón a logaritmo base 2 (_log2_)
de las cuentas por millón.
```
mod<- model.matrix(~ 0 + sra_attribute.tissue+ sra_attribute.cell_id_sorted + assigned_gene_prop,
                   data = colData(rse_gene_SRP192782)
)
```
A continuación, se estima la relación-varianza promedio de los datos mediante la función `voom` del paquete `limma`.
```
vGene <- voom(dge, mod, plot = TRUE)
```
![](plots/voom.png)

Se ajustan los datos a la variable _tissue_. Este paso es crucial para identificar genes que muestran una expresión diferencial significativa entre diferentes condiciones.

```
eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP192782),
  sort.by = "none"
)
```
![](plots/tissuevoom.png)

Subsecuentemente se realiza un _volcano plot_ para observar los genes más diferencialmente expresados.
```
volcanoplot(eb_results, coef = 4, highlight = 5, names = de_results$gene_name, col = "cornflowerblue", hl.col="darkorange")
```
![](plots/volcano.png)
Con base al _volcano plot_ se observan los cinco genes con mayor expresión diferencial: Xist, Foxp3, Klrk1, Ctla4 y Flicr. Profundizando en **Foxp3** ya que es uno de los genes más expresados, se encontró que es un factor regulador de la transcripción, que participa directamente en la función de las células reguladoras T CD4+ humanas y murinas. Lo cual hace sentido puesto que es un modelo murino con adenocarcinoma pulmonar. 

###  Visualización 
Por último, el script [07-visualizacion.R](https://github.com/anagarme/ProyectoFinal_RNAseq_2024/blob/78761415c1a4c75370a8f6e6d43de62496aba442/R/07-visualizacion.R) permite visualizar los 50 genes con mayor expresión diferencial.
![](plots/heatmap.png)



## Referencias 
Li A, Herbst RH, Canner D, Schenkel JM, Smith OC, Kim JY, Hillman M, Bhutkar A, Cuoco MS, Rappazzo CG, Rogers P, Dang C, Jerby-Arnon L, Rozenblatt-Rosen O, Cong L, Birnbaum M, Regev A, Jacks T. IL-33 Signaling Alters Regulatory T Cell Diversity in Support of Tumor Development. Cell Rep. 2019 Dec 3;29(10):2998-3008.e8. doi: 10.1016/j.celrep.2019.10.120. PMID: 31801068; PMCID: PMC6990979.
