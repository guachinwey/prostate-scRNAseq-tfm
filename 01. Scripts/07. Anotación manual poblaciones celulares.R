################################################################################
# 7. Anotación manual de clústeres y UMAPs por condición.
# - Carga del objeto anotado automáticamente (ScType).
# - Definición del mapeo clúster → población celular (anotación manual).
# - Incorporación de la anotación manual al meta.data (CellType_Manual).
# - Guardado del objeto anotado manualmente.
# - UMAP global y estratificado por condición (Healthy, PCa, CRPC).
################################################################################


library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)  


################################################################################
# 1) Carga del objeto anotado automáticamente (ScType)
################################################################################


# Se parte del objeto anotado automáticamente (ScType), que contiene:
# - Normalización SCT.
# - Clustering (seurat_clusters).
# - Anotación automática preliminar (CellType).
# - Coordenadas UMAP y metadatos.
merged_all_annotated <- readRDS(
  "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_annotated.rds"
)


################################################################################
# 2) Definición de anotación manual clúster → tipo celular
################################################################################


# Mapeo clúster → nombre de población celular, definido a partir de:
# - Genes marcadores.
# - Anotación automática.
# - Revisión bibliográfica y conocimiento biológico.
cluster_to_celltype_manual <- c(
  "1"  = "Benign luminal epithelial",
  "2"  = "Luminal metabolic epithelial",
  "3"  = "Secretory luminal epithelial",
  "4"  = "Macrophages",
  "5"  = "T cells",
  "6"  = "Basal epithelial",
  "7"  = "Endothelial cells",
  "8"  = "CAFs", # Cancer-associated fibroblasts
  "9"  = "Pericytes / smooth muscle cells",
  "10" = "Cycling / proliferating cells",
  "11" = "Luminal KLK11/12+",
  "12" = "Luminal MALAT1+",
  "13" = "Mucinous luminal",
  "14" = "Neuroendocrine / EMT-like",
  "15" = "Mast cells",
  "16" = "B cells",
  "17" = "CRISP3-high luminal",
  "18" = "PCA3-high epithelial"
)

# Usamos los Idents (seurat_clusters) como etiquetas de clúster por célula.
clusters_vec <- as.character(merged_all_annotated$seurat_clusters)

# Mapeamos cada clúster a su tipo celular manual.
celltype_manual <- cluster_to_celltype_manual[clusters_vec]

# Cualquier clúster no mapeado se asigna como "Unknown".
celltype_manual[is.na(celltype_manual)] <- "Unknown"

# Aseguramos que el vector tenga nombres de célula (colnames)
names(celltype_manual) <- colnames(merged_all_annotated)

# Añadimos el vector como columna explícita en el meta.data:
# - Nueva columna: CellType_Manual.
merged_all_annotated <- AddMetaData(
  merged_all_annotated,
  metadata = celltype_manual,
  col.name = "CellType_Manual"
)

# Por seguridad, cualquier entrada NA en CellType_Manual se reemplaza por "Unknown".
merged_all_annotated$CellType_Manual[is.na(merged_all_annotated$CellType_Manual)] <- "Unknown"

# Guardar objeto con anotación manual incorporada.
saveRDS(
  merged_all_annotated,
  file = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_annotated_manual.rds"
)


################################################################################
# 3) UMAP con anotación manual (global y por condición)
################################################################################


merged_all_annotated <- readRDS(
  "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_annotated_manual.rds"
)

# UMAP global con anotación manual
p_umap_manual <- DimPlot(
  merged_all_annotated,
  reduction = "umap",
  group.by  = "CellType_Manual",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP – Poblaciones celulares (anotación manual)")


# UMAP sólo Healthy
obj_healthy <- subset(merged_all_annotated, subset = condition == "Healthy")

p_healthy <- DimPlot(
  obj_healthy,
  reduction = "umap",
  group.by  = "CellType_Manual",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP – Poblaciones celulares (Healthy)")


# UMAP sólo PCa.
obj_pca <- subset(merged_all_annotated, subset = condition == "PCa")

p_pca <- DimPlot(
  obj_pca,
  reduction = "umap",
  group.by  = "CellType_Manual",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP – Poblaciones celulares (PCa)")


# UMAP sólo CRPC.
obj_crpc <- subset(merged_all_annotated, subset = condition == "CRPC")

p_crpc <- DimPlot(
  obj_crpc,
  reduction = "umap",
  group.by  = "CellType_Manual",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP – Poblaciones celulares (CRPC)")


################################################################################
# 4) Exportar UMAPs a PDF
################################################################################


pdf("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/UMAP_population_manual.pdf",
    width = 10, height = 8)
print(p_umap_manual)
print(p_healthy)
print(p_pca)
print(p_crpc)
dev.off()