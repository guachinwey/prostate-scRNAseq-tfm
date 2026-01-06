################################################################################
# 12. Heatmaps de DEGs (TOP UP/DOWN) por tipo celular y comparación
# - Lectura de resultados DEG desde Excel (uno por comparación).
# - Filtro de genes técnicos (MT-, RPL/RPS, HB).
# - Selección de genes más significativos (TOP UP/DOWN) por tipo celular.
# - Cálculo de expresión promedio (AverageExpression, datos SCT absolutos).
# - Generación de heatmaps con expresión por condición (2 columnas).
# - Exportación en PDF.
################################################################################

library(Seurat)
library(dplyr)
library(openxlsx)
library(pheatmap)
library(stringr)

################################################################################
# 1) Parámetros y rutas
################################################################################

# Objeto anotado con tipos celulares
obj <- readRDS("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_annotated_manual.rds")

# Validación de columnas requeridas
stopifnot("CellType_Manual" %in% colnames(obj@meta.data))
stopifnot("condition" %in% colnames(obj@meta.data))

# Archivos Excel con resultados DEG por comparación
deg_files <- list(
  Healthy_vs_PCa  = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/DEGs_by_CellType/DEGs_Healthy_vs_PCa_by_CellType.xlsx",
  Healthy_vs_CRPC = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/DEGs_by_CellType/DEGs_Healthy_vs_CRPC_by_CellType.xlsx",
  PCa_vs_CRPC     = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/DEGs_by_CellType/DEGs_PCa_vs_CRPC_by_CellType.xlsx"
)

# Carpeta de salida para los PDFs
out_dir <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/DEG_heatmaps"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Tipos celulares seleccionados
celltypes_to_plot <- c(
  "Benign luminal epithelial",
  "Macrophages",
  "T cells"
)

# Parámetros DEG
padj_thr  <- 0.05
logfc_thr <- 0.25

# Nº de genes TOP a mostrar
top_up   <- 15
top_down <- 15

# Activar el filtro de genes técnicos
filter_technical_genes <- TRUE


################################################################################
# 2) Funciones auxiliares
################################################################################

# Filtro para eliminar genes técnicos no informativos
filter_deg_genes <- function(df) {
  if (!"gene" %in% colnames(df)) df$gene <- rownames(df)
  
  df %>%
    filter(!is.na(gene), gene != "") %>%
    filter(!grepl("^MT-", gene)) %>%
    filter(!grepl("^RPL|^RPS", gene)) %>%
    filter(!grepl("^HB", gene))
}

# Selección de genes TOP UP y TOP DOWN (por log2FC y padj)
select_top_up_down <- function(df, top_up = 15, top_down = 15) {
  stopifnot("avg_log2FC" %in% colnames(df))
  stopifnot("p_val_adj" %in% colnames(df))
  if (!"gene" %in% colnames(df)) df$gene <- rownames(df)
  
  df_sig <- df %>%
    filter(!is.na(p_val_adj)) %>%
    filter(p_val_adj < padj_thr) %>%
    filter(abs(avg_log2FC) >= logfc_thr)
  
  up <- df_sig %>%
    filter(avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = top_up)
  
  down <- df_sig %>%
    filter(avg_log2FC < 0) %>%
    arrange(avg_log2FC) %>%
    slice_head(n = top_down)
  
  list(
    up_genes   = unique(up$gene),
    down_genes = unique(down$gene),
    all_genes  = unique(c(up$gene, down$gene))
  )
}

# Función para crear heatmap de expresión promedio (AverageExpression)
make_heatmap_avgexpr <- function(obj, celltype, group1, group2, genes, title, out_pdf) {
  
  obj_ct <- subset(obj, subset = CellType_Manual == celltype & condition %in% c(group1, group2))
  
  if (ncol(obj_ct) == 0) {
    return(invisible(NULL))
  }
  
  DefaultAssay(obj_ct) <- "SCT"
  
  genes <- genes[genes %in% rownames(obj_ct)]
  if (length(genes) < 2) {
     return(invisible(NULL))
  }
  
  avg <- AverageExpression(
    obj_ct,
    group.by = "condition",
    assays   = "SCT",
    slot     = "data",
    features = genes,
    verbose  = FALSE
  )$SCT
  
  # Asegurar orden de columnas
  avg <- avg[, intersect(c(group1, group2), colnames(avg)), drop = FALSE]
  
  # Exportar heatmap a PDF
  pdf(out_pdf, width = 7, height = 9)
  pheatmap(
    avg,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 8,
    main = title
  )
}

################################################################################
# 3) Bucle principal por comparación y tipo celular
################################################################################

for (comp in names(deg_files)) {
  
  xlsx <- deg_files[[comp]]
  if (!file.exists(xlsx)) {
    next
  }
  
  parts <- strsplit(comp, "_vs_")[[1]]
  group1 <- parts[1]
  group2 <- parts[2]
  
  sheets <- getSheetNames(xlsx)
  for (ct in celltypes_to_plot) {
     if (!ct %in% sheets) {
      next
    }
    
    # Leer hoja correspondiente
    df <- read.xlsx(xlsx, sheet = ct)
    if (!"gene" %in% colnames(df)) df$gene <- rownames(df)
    
    # Filtrar genes técnicos si aplica
    if (filter_technical_genes) {
      df <- filter_deg_genes(df)
    }
    
    # Obtener genes seleccionados
    top_lists <- select_top_up_down(df, top_up = top_up, top_down = top_down)
    genes_use <- top_lists$all_genes
    if (length(genes_use) < 2) {
      next
    }
    
    out_pdf <- file.path(
      out_dir,
      paste0("Heatmap_", comp, "_", str_replace_all(ct, "[^A-Za-z0-9]+", "_"), ".pdf")
    )
    
    # Generar heatmap
    make_heatmap_avgexpr(
      obj      = obj,
      celltype = ct,
      group1   = group1,
      group2   = group2,
      genes    = genes_use,
      title    = paste0(ct, " — ", group1, " vs ", group2, "\nAverageExpression (SCT data) | TOP UP/DOWN"),
      out_pdf  = out_pdf
    )
  }
}