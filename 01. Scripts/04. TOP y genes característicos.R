################################################################################
# 04. Detección de marcadores y anotación de poblaciones celulares.
# Objetivo general del script:
# - Partir del objeto integrado (merged_all) con clustering ya realizado.
# - Detectar genes marcadores por clúster usando el assay SCT.
# - Construir tablas TOP2 y TOP10 por clúster para su uso en la memoria.
# - Generar heatmaps globales y por bloques de clústeres.
# - Preparar el material necesario para la anotación (automática y manual).
################################################################################


library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(writexl)
library(readxl)
library(openxlsx)
library(HGNChelper)

# Al usar ScType:
# - sctype_score_.R: función principal para scoring de tipos celulares.
# - gene_sets_prepare.R: preparación de sets de genes específicos por tejido.
source("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Funciones/sctype_score_.R")
source("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Funciones/gene_sets_prepare.R")


################################################################################
# 1) Cargar el objeto integrado con clustering ya realizado
################################################################################


merged_all <- readRDS(
  "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_integrated.rds"
)

# Asegurarnos de que el clustering está en Idents:
Idents(merged_all) <- merged_all$seurat_clusters
levels(Idents(merged_all))
length(unique(Idents(merged_all)))   # debería dar 17

# Trabar sobre el assay normalizado con SCT:
DefaultAssay(merged_all) <- "SCT"


################################################################################
# 2) Preparar objeto para detección de marcadores (SCT)
################################################################################


# Esta función ajusta internamente el modelo de varianza de SCT para que 
# FindAllMarkers funcione correctamente con este tipo de normalización.
merged_all <- PrepSCTFindMarkers(merged_all)


################################################################################
# 3) Detección de genes marcadores por clúster
################################################################################


# Encontrar marcadores de todos los clústeres frente al resto:
#   - only.pos = TRUE → solo genes sobreexpresados
#   - min.pct = 0.25 → al menos en el 25% de las células del clúster
#   - logfc.threshold = 0.25 → cambio mínimo de expresión (log2FC)
markers_all <- FindAllMarkers(
  merged_all,
  assay           = "SCT",
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

# Guardar la tabla completa de marcadores en Excel
output_markers_xlsx <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/FindAllMarkers_SCT.xlsx"
write_xlsx(markers_all, output_markers_xlsx)


################################################################################
# 4) TOP2 y TOP10 genes por clúster
################################################################################


# Top 2 genes por clúster (para resumen/tablas)
top2 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

output_top2 <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/TOP2_markers.xlsx"
write_xlsx(top2, output_top2)

# Top 10 genes por clúster (para heatmaps)
top10 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

output_top10 <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/TOP10_markers.xlsx"
write_xlsx(top10, output_top10)

cat("✅ TOP2 y TOP10 por clúster guardados en Excel.\n")


################################################################################
# 5. Heatmap global de marcadores TOP10
################################################################################


# Heatmap con todos los clústeres y sus top10 genes
pdf("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/heatmap_TOP10_all_clusters.pdf",
    width = 10, height = 10)

DoHeatmap(
  merged_all,
  features = unique(top10$gene),
  group.by = "seurat_clusters"
) + NoLegend()

dev.off()


################################################################################
# 6. Heatmaps por bloques de clústeres.
################################################################################


# Lista de identidades de clúster 
cluster_ids <- levels(Idents(merged_all))
cluster_ids

# Divide la lista de clústeres en bloques de tamaño 5
chunk_size <- 5
cluster_chunks <- split(cluster_ids, ceiling(seq_along(cluster_ids) / chunk_size))

pdf("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/heatmap_TOP10_by_cluster_blocks.pdf",
    width = 10, height = 10)

for (cl_group in cluster_chunks) {
  genes_group <- top10$gene[top10$cluster %in% cl_group]
  p <- DoHeatmap(
    merged_all,
    group.by = "seurat_clusters",
    features = genes_group,
    cells    = WhichCells(merged_all, idents = cl_group)
  ) + ggtitle(paste("Clústeres:", paste(cl_group, collapse = ", "))) +
    NoLegend()
  print(p)
}

dev.off()