################################################################################
# 7.1. Revisión exploratoria de la población "Luminal MALAT1+"
# - Visualización de expresión de MALAT1 junto con métricas de QC
# - Comparación con otras poblaciones luminales diferenciadas
# - Análisis de asociación entre MALAT1 y complejidad/transcripción celular
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)

################################################################################
# 1) Cargar objeto anotado
################################################################################

merged_all_annotated <- readRDS(
  "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_annotated.rds"
)


################################################################################
# 2) Gráfico de violín: MALAT1 y métricas de QC por tipo celular
################################################################################
# - Se incluyen: expresión de MALAT1, número de genes detectados (nFeature),
#   número de transcritos (nCount), % mitocondrial, % ribosomal
# - Permite ver si MALAT1+ se asocia a perfiles de calidad celular atípicos

VlnPlot(
  merged_all_annotated,
  features = c("MALAT1", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  group.by = "CellType_Manual",
  pt.size = 0
)

################################################################################
# 3) Comparación centrada en poblaciones luminales diferenciadas
################################################################################
# - Subset de poblaciones luminales (incluye "Luminal MALAT1+")
# - Análisis más fino dentro del compartimento luminal
# - Objetivo: detectar si MALAT1+ representa una transición o artefacto

luminal_set <- c(
  "Benign luminal epithelial",
  "Luminal metabolic epithelial",
  "Secretory luminal epithelial",
  "Luminal KLK11/12+",
  "Luminal MALAT1+",
  "CRISP3-high luminal",
  "PCA3-high epithelial",
  "Mucinous luminal"
)

obj_lum <- subset(merged_all_annotated, subset = CellType_Manual %in% luminal_set)

# Violin plot dentro de luminales seleccionadas
VlnPlot(
  obj_lum,
  features = c("MALAT1", "nFeature_RNA", "percent.mt"),
  group.by = "CellType_Manual",
  pt.size = 0
)

################################################################################
# 4) Dispersión: relación entre MALAT1 y complejidad transcriptómica
################################################################################
# - Evaluar si hay correlación entre expresión de MALAT1 y:
#   - Nº genes detectados (nFeature_RNA)
#   - Nº total de UMI (nCount_RNA)
# - Puede revelar si MALAT1 se asocia a hipertranscripción o artefactos

FeatureScatter(obj_lum, feature1 = "MALAT1", feature2 = "nFeature_RNA")
FeatureScatter(obj_lum, feature1 = "MALAT1", feature2 = "nCount_RNA")