################################################################################
# 3. Integración, reducción dimensional y clustering.
# - Carga de objetos Seurat filtrados por QC.
# - Integración de muestras en un único objeto.
# - Normalización con SCTransform.
# - Reducción dimensional (PCA + Harmony).
# - Construcción de grado de vecinos y clustering (Leiden).
# - Cálculo de UMAP y generación de figuras.
################################################################################


library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)


################################################################################
# 1) Integración de muestras y normalización.
################################################################################


# Cargar todos los objetos QC ya filtrados.
qc_dir <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_QC"

HEALTHY_1_QC <- readRDS(file.path(qc_dir, "HEALTHY_1_QC.rds"))
HEALTHY_2_QC <- readRDS(file.path(qc_dir, "HEALTHY_2_QC.rds"))
HEALTHY_3_QC <- readRDS(file.path(qc_dir, "HEALTHY_3_QC.rds"))
PCA_QC       <- readRDS(file.path(qc_dir, "PCA_QC.rds"))
CRPC_1_QC    <- readRDS(file.path(qc_dir, "CRPC_1_QC.rds"))
CRPC_2_QC    <- readRDS(file.path(qc_dir, "CRPC_2_QC.rds"))
CRPC_3_QC    <- readRDS(file.path(qc_dir, "CRPC_3_QC.rds"))

# Merge en un único objeto integrado que contiene todas las muestras:
# - Se combinan las células de Healthy, PCa y CRPC en un mismo objeto Seurat.
# - Se preservan los metadatos de condición y sample_id para análisis posteriores.

merged_all <- merge(
  HEALTHY_1_QC,
  y = list(
    HEALTHY_2_QC,
    HEALTHY_3_QC,
    PCA_QC,
    CRPC_1_QC,
    CRPC_2_QC,
    CRPC_3_QC),
  project = "PCa_progression")

# Comprobar distribución por condición y muestra:
table(merged_all$condition)
table(merged_all$sample_id)
ncol(merged_all)  # Número total de células combinadas

# Normalización con SCTransform:
# - Modela los conteos con uan distribución binomial negativa regularizada.
# - Identifica genes altamente variables y corrige efectos de profundidad.
# - Sustituye la necesidad de NormalizeData + FindVariableFeatures.
merged_all <- SCTransform(
  merged_all,
  verbose = FALSE)

# No es necesario realizar ScaleData en SCT, pero se mantiene por consistencia:
DefaultAssay(merged_all) <- "SCT"
merged_all <- ScaleData(
  merged_all,
  assay = "SCT",
  verbose = TRUE)

# Comprobaciones de que SCT se realizó correctamente:
DefaultAssay(merged_all)
Assays(merged_all)
merged_all[["SCT"]]
length(VariableFeatures(merged_all))
head(VariableFeatures(merged_all))


################################################################################
# 2) Reducción dimensional (PCA + Harmony)
################################################################################


# PCA sobre SCT para reducir dimensionalidad antes de Harmony
# - npcs = 50: se calculan hasta 50 componentes principales.
# Estas PCs se usarán como entrada para la corrección de bath con Harmony.

merged_all <- RunPCA(merged_all, npcs = 50, verbose = FALSE)

# Corrección de batch con Harmony:
# - group.by.vars = "sample_id": se corrige el efecto de la muestra.
# - Se busca eliminar variabilidad técnica presenvando las diferencias biológicas 
# entre las condiciones (Healthy / PCa / CRPC).

merged_all <- RunHarmony(
  merged_all,
  group.by.vars = "sample_id",
  verbose = TRUE) # Harmony reporta convergencia en la 6ª iteración,alcanzando correctamente
# el criterio de convergencia.


################################################################################
# 3) Clustering y UMAP (Grafo de vecinos y visualización)
################################################################################


# Vecinos + clustering:
# - Se contruye un grafo de vecinos más cercanos a partir de los embbedings
# corregidos por Harmony.
# Parámetros principales:
# - dims = 1:30: PCs utilizadas para capturar variabilidad biológica.
# - k.param = 30: tamaño del vecindario (número de vecinos por célula).
# - algorithm = 4: Leiden, mas robusto que Louvain para datos scrNA-seq.
# - resolution = 0.5: resolución intermedia que produce 18 clústeres.

merged_all <- FindNeighbors(merged_all, 
                            reduction = "harmony", 
                            dims = 1:30, 
                            k.param = 30,
                            verbose = TRUE)
merged_all <- FindClusters(merged_all, 
                           algorithm = 4,
                           resolution = 0.5,
                           verbose = TRUE)

# Recuento de clusters:
# - Se establece la identidad activa como seurat_clusters.
# - Se calcula el número total de clústers detectados.
length(unique(Idents(merged_all)))
Idents(merged_all) <- merged_all$seurat_clusters
num_clusters <- length(unique(Idents(merged_all)))
num_clusters

# UMAP para visualización:
# - Se calculan coordenadas UMAP usando las embeddings de Harmony.
# Parámetros:
# - dims = 1:20: subconjunto de PCs para estabilizar UMAP.
# - min.dist = 0.3: controla la separación entre clústers (0 = más compacto, 
# valores mayores = más disperso).
# - n.neighbors = 100: da más peso a la estructura global del dataset.

merged_all <- RunUMAP(
  merged_all,
  reduction = "harmony",
  dims = 1:20,
  min.dist = 0.3,
  n.neighbors = 100)

# Plots clave:
# p_clusters: UMAP coloreado por clúster (seurat_clusters).
# p_condition: UMAP coloreado por condición biológica (Healthy /PCa / CRPC).
# p_sample: UMAP coloreado por sample_ID (donante).
p_clusters <- DimPlot(merged_all, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP - Clusters")

p_condition <- DimPlot(merged_all, reduction = "umap", group.by = "condition") +
  ggtitle("UMAP - Condition (Healthy / PCa / CRPC)")

p_sample <- DimPlot(merged_all, reduction = "umap", group.by = "sample_id") +
  ggtitle("UMAP - sample_id (donante / biopsia)")

# Ruta de salida del PDF
output_pdf <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/UMAP_separate_pages.pdf"

# Exportar a PDF con cada gráfico en una página separada
pdf(file = output_pdf, width = 8, height = 7)  
print(p_clusters)
print(p_condition)
print(p_sample)
dev.off()

# Guardado de los gráficos individuales en PNG:
ggsave(
  filename = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/UMAP_clusters.png",
  plot = p_clusters, width = 8, height = 7, dpi = 300)

ggsave(
  filename = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/UMAP_condition.png",
  plot = p_condition, width = 8, height = 7, dpi = 300)

ggsave(
  filename = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/UMAP_sample.png",
  plot = p_sample, width = 8, height = 7, dpi = 300)

# Guardar el objeto integrado, el cuál incluirá:
# - Datos normalizados (SCT).
# - Embbedings PCA y Harmony.
# - Clústers (seurat_clusters).
# - Coordenadas UMAP y metadatos.
saveRDS(
  merged_all,
  file = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_integrated.rds")