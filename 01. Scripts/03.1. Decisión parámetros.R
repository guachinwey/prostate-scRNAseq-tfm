################################################################################
# 3.1. Justificación de parámetros de reducción dimensional y clustering.
# - Selección del número de PCs (Elbow + varianza explicada).
# - Evaluación de estabilidad de clústers según dims (ARI / NMI).
# - Exploración de k.param en FindNeighbors.
# - Comparación de distintas resoluciones de clustering.
# - UMAP comparativo entre resoluciones y condiciones.
################################################################################

# Carga de librerías necesarias para todo el script
library(Seurat)
library(dplyr)
library(mclust)    # adjustedRandIndex (ARI)
library(aricode)   # NMI
library(ggplot2)


################################################################################
# 1) Selección del número de componentes principales (PCs)
################################################################################

# ElbowPlot para explorar la caída de varianza explicada por PC
ElbowPlot(merged_all, ndims = 50)

# Cálculo de la varianza explicada y varianza acumulada
var_exp <- merged_all[["pca"]]@stdev^2
cumvar <- cumsum(var_exp) / sum(var_exp)

# Curva de varianza explicada acumulada
plot(cumvar, type = "l")
abline(h = 0.8, col = "red", lty = 2)  # línea de referencia en el 80% de varianza


################################################################################
# 2) Evaluación de estabilidad de clustering según número de PCs (dims)
################################################################################


# Conjuntos de PCs a probar:
# - 1:20, 1:25, 1:30, 1:35
# Estos se usarán en FindNeighbors / FindClusters para comparar estabilidad.
dims_list <- list(
  d20 = 1:20,
  d25 = 1:25,
  d30 = 1:30,
  d35 = 1:35)

# Parámetros fijos para todas las pruebas de dims
k_param    <- 30
resolution <- 0.5
algo       <- 4   # Leiden

set.seed(123)  # reproducibilidad

# Repetir vecinos + clustering para cada dims
for (name in names(dims_list)) {
  dims_use <- dims_list[[name]]
  message("Calculando vecinos y clusters para: ", name,
          " (PCs = ", min(dims_use), ":", max(dims_use), ")")
  
  # Vecinos usando un subconjunto de PCs de Harmony
  merged_all <- FindNeighbors(
    merged_all,
    reduction = "harmony",
    dims      = dims_use,
    k.param   = k_param,
    verbose   = FALSE
  )
  
  # Clustering Leiden con resolución fija
  merged_all <- FindClusters(
    merged_all,
    algorithm  = algo,
    resolution = resolution,
    verbose    = FALSE
  )
  
  # Guardar los clústers resultantes en una columna nueva del meta.data
  # para poder compararlos posteriormente entre distintas dims.
  colname <- paste0("clusters_dims", max(dims_use))
  merged_all@meta.data[[colname]] <- Idents(merged_all)
}

# Ahora el objeto contiene columnas:
#  - clusters_dims20, clusters_dims25, clusters_dims30, clusters_dims35
# que representan las asignaciones de clúster para cada elección de PCs.

# Calcular matrices de ARI y NMI entre dimsets
cluster_cols <- c("clusters_dims20", "clusters_dims25",
                  "clusters_dims30", "clusters_dims35")

# Matrices vacías para ARI y NMI entre todas las combinaciones de dimsets
ari_mat <- matrix(NA,
                  nrow = length(cluster_cols),
                  ncol = length(cluster_cols),
                  dimnames = list(cluster_cols, cluster_cols))

nmi_mat <- ari_mat

for (i in seq_along(cluster_cols)) {
  for (j in seq_along(cluster_cols)) {
    c1 <- merged_all@meta.data[[cluster_cols[i]]]
    c2 <- merged_all@meta.data[[cluster_cols[j]]]
    
    # ARI: mide concordancia entre dos particiones de clústeres
    ari_mat[i, j] <- adjustedRandIndex(c1, c2)
    
    # NMI: mide información mutua normalizada entre particiones
    nmi_mat[i, j] <- NMI(c1, c2)
  }
}

# Tablas de estabilidad de clustering según dims:
ari_mat
nmi_mat


################################################################################
# 3) Exploración del parámetro k.param en FindNeighbors
################################################################################

# Aquí se evalúa cómo cambia el número de clústers si variamos k.param:
# - k = 20, 30, 40
# Manteniendo:
# - dims = 1:30, resolution = 0.5, algorithm = 4 (Leiden)

for (k in c(20, 30, 40)) {
  tmp <- FindNeighbors(merged_all, reduction = "harmony", dims = 1:30,
                       k.param = k)
  tmp <- FindClusters(tmp, resolution = 0.5, algorithm = 4)
  cat("k =", k, "→ clusters:", length(unique(Idents(tmp))), "\n")
}


################################################################################
# 4) Exploración de resoluciones de clustering y UMAP asociado
################################################################################

# Vecinos fijos usando Harmony con dims 1:30 y k.param = 30
merged_all <- FindNeighbors(
  merged_all,
  reduction = "harmony",
  dims = 1:30,        
  k.param = 30,       
  verbose = TRUE
)

# Se prueban varias resoluciones y se guardan en metadatos
resolutions_to_try <- c(0.2, 0.3, 0.5)

for (res in resolutions_to_try) {
  merged_all <- FindClusters(
    merged_all,
    algorithm = 4,     
    resolution = res,
    verbose = FALSE)
  colname <- paste0("clusters_res", res)
  merged_all@meta.data[[colname]] <- Idents(merged_all)
}

# Elección de resolución final para análisis (ejemplo: 0.5)
Idents(merged_all) <- merged_all$clusters_res0.5  # por ejemplo 0.3

num_clusters <- length(unique(Idents(merged_all)))
cat("Número de clusters con res=0.5:", num_clusters, "\n")

# Recalcular UMAP con la configuración elegida:
# - dims = 1:30 (más PCs que en el script principal, a modo de comparación).
# - n.neighbors = 30 (más local que 100, estructura más detallada).
merged_all <- RunUMAP(
  merged_all,
  reduction    = "harmony",
  dims         = 1:30,    # antes 1:20
  min.dist     = 0.3,
  n.neighbors  = 30       # antes 100, más local
)


################################################################################
# 5) Visualizaciones comparativas de resoluciones y metadatos
################################################################################

# UMAP coloreado por clústers con distintas resoluciones:
p_res02 <- DimPlot(
  merged_all,
  reduction = "umap",
  group.by  = "clusters_res0.2",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP - clusters res=0.2")

p_res03 <- DimPlot(
  merged_all,
  reduction = "umap",
  group.by  = "clusters_res0.3",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP - clusters res=0.3")

p_res05 <- DimPlot(
  merged_all,
  reduction = "umap",
  group.by  = "clusters_res0.5",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP - clusters res=0.5")

# UMAP coloreado por condición biológica (Healthy / PCa / CRPC)
p_condition <- DimPlot(
  merged_all,
  reduction = "umap",
  group.by  = "condition"
) + ggtitle("UMAP - Condition (Healthy / PCa / CRPC)")

# UMAP coloreado por sample_id (donante / biopsia)
p_sample <- DimPlot(
  merged_all,
  reduction = "umap",
  group.by  = "sample_id"
) + ggtitle("UMAP - sample_id (donante / biopsia)")


################################################################################
# 6) Exportar figuras a PDF (cada plot en una página)
################################################################################

# Guardar todas las UMAPs en un único PDF, cada una en una página.
output_pdf <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/UMAP_multires.pdf"
pdf(file = output_pdf, width = 8, height = 7)
print(p_res02)
print(p_res03)
print(p_res05)
print(p_condition)
print(p_sample)
dev.off()