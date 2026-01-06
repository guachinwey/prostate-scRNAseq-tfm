################################################################################
# 9. Análisis de genes diferencialmente expresados (DEGs) por tipo celular
# - Carga del objeto anotado con CellType_Manual.
# - Activación del assay SCT para análisis de expresión diferencial.
# - Comparaciones entre condiciones: Healthy vs PCa, Healthy vs CRPC, PCa vs CRPC.
# - Resultados exportados a Excel (una hoja por tipo celular).
################################################################################

library(Seurat)
library(dplyr)
library(openxlsx)

################################################################################
# 1) Definición de parámetros y rutas de entrada/salida
################################################################################

input_rds <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_annotated_manual.rds"

out_dir_deg <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/DEGs_by_CellType"
dir.create(out_dir_deg, showWarnings = FALSE, recursive = TRUE)

# Comparaciones entre condiciones a ejecutar
comparisons <- list(
  Healthy_vs_PCa  = c("Healthy", "PCa"),
  Healthy_vs_CRPC = c("Healthy", "CRPC"),
  PCa_vs_CRPC     = c("PCa", "CRPC")
)

# Umbrales para DE
min_cells_per_group <- 25   # mínimo de células por grupo
min_pct             <- 0.10 # expresión mínima del gen
logfc_thr           <- 0.25 # umbral de log2FC


################################################################################
# 2) Cargar objeto y realizar chequeos
################################################################################

obj <- readRDS(input_rds)

# Validación de columnas necesarias
stopifnot("CellType_Manual" %in% colnames(obj@meta.data))
stopifnot("condition" %in% colnames(obj@meta.data))

# DE se realizará sobre el assay "SCT" (debe existir)
stopifnot("SCT" %in% names(obj@assays))
DefaultAssay(obj) <- "SCT"

# Asignación de identidades por condición
Idents(obj) <- "condition"

################################################################################
# 3) Función para encontrar DEGs dentro de un tipo celular específico
################################################################################

run_deg_for_celltype <- function(object,
                                 celltype_name,
                                 group1,
                                 group2,
                                 min_cells = 25,
                                 min_pct = 0.10,
                                 logfc_thr = 0.25) {
  
  # Subset del objeto por tipo celular
  sub <- subset(object, subset = CellType_Manual == celltype_name)
  
  # Chequeo de tamaño de muestra por grupo
  n1 <- sum(sub$condition == group1)
  n2 <- sum(sub$condition == group2)
  
  if (n1 < min_cells || n2 < min_cells) {
    message("   (skip) N insuficiente: ", group1, "=", n1, " | ", group2, "=", n2)
    return(NULL)
  }
  
  # Asignación de identidades para DE
  Idents(sub) <- "condition"
  
  # Verificación y selección del assay "SCT"
  if (!("SCT" %in% names(sub@assays))) {
    stop("No existe el assay 'SCT' en el objeto. Revisa DefaultAssay/SCTransform.")
  }
  DefaultAssay(sub) <- "SCT"
  
  # Preparar el objeto para DE (si falla, se sigue igual)
  suppressWarnings({
    sub <- tryCatch(
      PrepSCTFindMarkers(sub, assay = "SCT", verbose = FALSE),
      error = function(e) sub
    )
  })
  
  # Intento principal: FindMarkers usando Wilcoxon
  deg <- tryCatch(
    {
      FindMarkers(
        object = sub,
        ident.1 = group1,
        ident.2 = group2,
        assay = "SCT",
        test.use = "wilcox",
        min.pct = min_pct,
        logfc.threshold = logfc_thr
      )
    },
    error = function(e) {
      # En caso de error por diferencias en librerías, reintenta sin recorrect_umi
      msg <- conditionMessage(e)
      if (grepl("multiple models with unequal library sizes", msg, ignore.case = TRUE)) {
        message("   (retry) SCT multi-model detectado -> usando recorrect_umi = FALSE")
        return(
          FindMarkers(
            object = sub,
            ident.1 = group1,
            ident.2 = group2,
            assay = "SCT",
            test.use = "wilcox",
            min.pct = min_pct,
            logfc.threshold = logfc_thr,
            recorrect_umi = FALSE
          )
        )
      } else {
        stop(e)
      }
    }
  )
  
  # Formateo final de los resultados
  deg <- deg %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::arrange(p_val_adj, desc(avg_log2FC))
  
  return(deg)
}

################################################################################
# 4) Loop por comparaciones y exportación de resultados
################################################################################

celltypes <- sort(unique(obj$CellType_Manual))

for (comp_name in names(comparisons)) {
  
  groups <- comparisons[[comp_name]]
  g1 <- groups[1]
  g2 <- groups[2]
  
  message("\n==============================")
  message("Comparación: ", comp_name, " (", g1, " vs ", g2, ")")
  message("==============================\n")
  
  deg_list <- list()
  
  for (ct in celltypes) {
    message(" - DEGs en: ", ct)
    
    deg_ct <- run_deg_for_celltype(
      object = obj,
      celltype_name = ct,
      group1 = g1,
      group2 = g2,
      min_cells = min_cells_per_group,
      min_pct = min_pct,
      logfc_thr = logfc_thr
    )
    
    if (!is.null(deg_ct) && nrow(deg_ct) > 0) {
      deg_list[[ct]] <- deg_ct
    }
  }
  
  if (length(deg_list) == 0) {
    next
  }
  
  # Renombrar hojas de Excel evitando caracteres problemáticos
  safe_names <- names(deg_list) %>%
    gsub("[\\[\\]\\*\\?\\:/\\\\]", "_", .) %>%
    substr(1, 28)
  names(deg_list) <- safe_names
  
  out_xlsx <- file.path(out_dir_deg, paste0("DEGs_", comp_name, "_by_CellType.xlsx"))
  openxlsx::write.xlsx(deg_list, file = out_xlsx)
}