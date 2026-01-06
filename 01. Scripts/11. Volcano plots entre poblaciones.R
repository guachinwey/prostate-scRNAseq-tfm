################################################################################
# Volcano plots
# - Comparaciones: Healthy vs PCa, Healthy vs CRPC, PCa vs CRPC.
# - Inversión del eje X: log2FC positivo → sobreexpresión en grupo 2.
# - Etiquetas de los genes más significativos.
# - Exportación en PDF.
################################################################################

library(openxlsx)
library(dplyr)
library(ggplot2)
library(stringr)

# Archivos Excel con DEGs por tipo celular
deg_files <- list(
  Healthy_vs_PCa  = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/DEGs_by_CellType/DEGs_Healthy_vs_PCa_by_CellType.xlsx",
  Healthy_vs_CRPC = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/DEGs_by_CellType/DEGs_Healthy_vs_CRPC_by_CellType.xlsx",
  PCa_vs_CRPC     = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/DEGs_by_CellType/DEGs_PCa_vs_CRPC_by_CellType.xlsx"
)

# Directorio de salida para los volcano plots
out_dir <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/Volcano_plots_DEGs_corrected"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Tipos celulares previamente seleccionados como robustos
celltypes_to_plot <- c(
  "Benign luminal epithelial",
  "Luminal metabolic epithelial",
  "Macrophages",
  "T cells",
  "Endothelial cells"
)

# Parámetros para DE y visualización
padj_thr  <- 0.05   # umbral de p-valor ajustado
logfc_thr <- 0.25   # umbral de log2FC
top_label <- 10     # nº máximo de genes a etiquetar

################################################################################
# Función para generar volcano plot con eje X invertido (grupo2 vs grupo1)
################################################################################

make_volcano <- function(df, title, subtitle = NULL) {
  
  # Asegurar columna gene
  if (!"gene" %in% colnames(df)) df$gene <- rownames(df)
  stopifnot("p_val_adj" %in% colnames(df))
  stopifnot("avg_log2FC" %in% colnames(df))
  
  # Invertir dirección del log2FC (interpretable como: mayor expresión en grupo 2)
  df <- df %>%
    mutate(
      avg_log2FC = -avg_log2FC,
      neglog10_padj = -log10(p_val_adj + 1e-300),
      sig = (p_val_adj < padj_thr) & (abs(avg_log2FC) >= logfc_thr)
    )
  
  # Selección de genes a rotular (más significativos)
  label_df <- df %>%
    filter(sig) %>%
    arrange(p_val_adj) %>%
    slice_head(n = top_label)
  
  # Crear gráfico
  ggplot(df, aes(x = avg_log2FC, y = neglog10_padj)) +
    geom_point(aes(shape = sig), alpha = 0.6) +
    geom_vline(xintercept = c(-logfc_thr, logfc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    geom_text(
      data = label_df,
      aes(label = gene),
      size = 3,
      vjust = 1.2,
      check_overlap = TRUE
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "avg_log2FC (positivo = mayor en grupo 2)",
      y = "-log10(p_val_adj)"
    )
}

################################################################################
# Loop por comparaciones y tipos celulares
################################################################################

for (comp in names(deg_files)) {
  
  xlsx <- deg_files[[comp]]
  sheets <- getSheetNames(xlsx)
  
  # Extraer nombres de las condiciones (grupo 1 vs grupo 2)
  conds <- unlist(strsplit(comp, "_vs_"))
  group1 <- conds[1]
  group2 <- conds[2]
  
  for (ct in celltypes_to_plot) {
    
    if (!ct %in% sheets) {
      next
    }
    
    # Cargar tabla de DEGs
    df <- read.xlsx(xlsx, sheet = ct)
    
    # Volcano plot con dirección corregida del log2FC
    p <- make_volcano(
      df,
      title = paste0(ct, " — ", comp),
      subtitle = paste0("Eje X: +log2FC → mayor en ", group2)
    )
    
    # Guardar como PDF con nombre limpio
    out_pdf <- file.path(
      out_dir,
    )
    
    ggsave(out_pdf, plot = p, width = 7, height = 5)
  }
}