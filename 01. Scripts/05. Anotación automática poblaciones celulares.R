################################################################################
# 5. Anotación automática de tipos celulares y UMAPs por condición.
# - Carga de librerías y funciones ScType.
# - Cálculo de scores de tipo celular con ScType sobre datos SCT.
# - Resumen de anotación por clúster (top tipos celulares).
# - Filtrado de anotaciones débiles → "Unknown".
# - Incorporación de la anotación al meta.data (columna CellType).
# - Visualización UMAP global y estratificada por condición (Healthy, PCa, CRPC).
################################################################################


library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)  
library(limma)


# Funciones auxiliares para ScType:
# - sctype_score_.R: calcula los scores de tipo celular por gen set.
# - gene_sets_prepare.R: prepara las listas de genes positivos/negativos
#   a partir de la base de datos ScTypeDB.
source("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Funciones/sctype_score_.R")
source("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Funciones/gene_sets_prepare.R")


################################################################################
# 1) Anotación de tipos celulares con ScType
################################################################################


# Base de datos ScType:
# - db_: ruta al fichero Excel con firmas de ScType reducidas (ScTypeDB_short.xlsx).
# - tissue: tejido de interés; en este caso "Prostate" para usar firmas específicas.
db_ <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/ScTypeDB_short.xlsx"
tissue <- "Prostate"

# Preparar listas de genes a partir de la base de datos:
# - gs_positive: lista de genes positivos por tipo celular para próstata.
# - gs_negative (no usado aquí, se pasa como NULL).
gs_list <- gene_sets_prepare(db_, tissue)

# Usamos la matriz escalada de SCT:
# - scale.data contiene la expresión normalizada y escalada por SCTransform.
sct_scaled <- GetAssayData(merged_all, assay = "SCT", slot = "scale.data")

# Cálculo de scores ScType:
# - scRNAseqData: matriz de expresión (genes x células).
# - scaled = TRUE: se indica que los datos ya están escalados.
# - gs: lista de firmas positivas por tipo celular.
# - gs2: firmas negativas (no usadas -> NULL).
es.max <- sctype_score(
  scRNAseqData = sct_scaled,
  scaled       = TRUE,
  gs           = gs_list$gs_positive,
  gs2          = NULL
)

# Resumen por clúster:
# - Para cada clúster se suman los scores de todas sus células
#   para cada tipo celular.
# - Se ordenan de mayor a menor y se guardan las 10 mejores etiquetas
#   para facilitar la interpretación.
sctype_results <- do.call(
  "rbind",
  lapply(unique(merged_all$seurat_clusters), function(cl) {
    cells_cl  <- rownames(merged_all@meta.data)[merged_all$seurat_clusters == cl]
    scores_cl <- rowSums(es.max[, cells_cl, drop = FALSE])
    scores_cl <- sort(scores_cl, decreasing = TRUE)
    head(
      data.frame(
        cluster = cl,
        type    = names(scores_cl),
        scores  = scores_cl,
        ncells  = length(cells_cl)
      ),
      10
    )
  })
)

# Seleccionamos la mejor etiqueta (máximo score) por clúster.
sctype_best <- sctype_results %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = scores)

# Filtrado de etiquetas débiles:
# - Si el score máximo es < ncells/4, se considera poco informativo.
# - En esos casos, se fuerza la etiqueta a "Unknown".
sctype_best$type[as.numeric(as.character(sctype_best$scores)) < sctype_best$ncells / 4] <- "Unknown"

# Guardar anotación ScType en Excel para revisión manual posterior.
output_annot_xlsx <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/AnnotatedClusters_scType.xlsx"
write_xlsx(sctype_best, output_annot_xlsx)

# Añadir anotación al meta.data:
# - Se crea una columna "CellType" inicializada como "Unknown".
# - Para cada clúster se asigna el tipo celular propuesto por ScType.
merged_all$CellType <- "Unknown"
for (j in unique(sctype_best$cluster)) {
  cl_type <- sctype_best[sctype_best$cluster == j, ]
  merged_all$CellType[merged_all$seurat_clusters == j] <- as.character(cl_type$type[1])
}

# Guardar objeto anotado con la columna CellType incorporada.
saveRDS(
  merged_all,
  file = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_annotated.rds"
)


################################################################################
# 2) Visualización UMAP con anotación de tipos celulares
################################################################################

# Se visualiza:
# - UMAP global con todas las poblaciones anotadas.
# - UMAP por condición (Healthy, PCa, CRPC) de forma independiente.

# UMAP global con anotación automática (ScType)
p1 <- DimPlot(
  merged_all,
  reduction = "umap",
  group.by  = "CellType",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP – Poblaciones celulares anotadas")

# UMAP solo Healthy
obj_healthy <- subset(merged_all, subset = condition == "Healthy")

p_healthy <- DimPlot(
  obj_healthy,
  reduction = "umap",
  group.by  = "CellType",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP – Poblaciones celulares (Healthy)")

# UMAP solo PCa
obj_pca <- subset(merged_all, subset = condition == "PCa")

p_pca <- DimPlot(
  obj_pca,
  reduction = "umap",
  group.by  = "CellType",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP – Poblaciones celulares (PCa)")

# UMAP solo CRPC
obj_crpc <- subset(merged_all, subset = condition == "CRPC")

p_crpc <- DimPlot(
  obj_crpc,
  reduction = "umap",
  group.by  = "CellType",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP – Poblaciones celulares (CRPC)")


# Guardar todas las figuras en un único PDF
pdf("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/UMAP_population_auto.pdf",
    width = 10, height = 8)
print(p1)
print(p_healthy)
print(p_pca)
print(p_crpc)
dev.off()