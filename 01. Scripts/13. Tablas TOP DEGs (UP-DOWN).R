################################################################################
# 13. Tablas de TOP DEGs (UP/DOWN) por tipo celular y comparación
# - Una hoja por cada tipo celular en CellType_Manual, sin omisiones.
# - Si no hay hoja en el Excel original → se escribe mensaje (probable N insuficiente).
# - Si hay hoja pero no hay DEGs significativos → se indica en la hoja.
# - Si hay DEGs → se listan los TOP UP y DOWN según criterios definidos.
################################################################################

library(openxlsx)
library(dplyr)
library(stringr)

################################################################################
# 1) Parámetros y rutas
################################################################################

# Objeto anotado con tipos celulares (para extraer todos los CellType_Manual posibles)
obj <- readRDS("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_annotated_manual.rds")
stopifnot("CellType_Manual" %in% colnames(obj@meta.data))

# Lista completa de tipos celulares presentes
all_celltypes <- sort(unique(obj$CellType_Manual))

# Archivos Excel de entrada con DEGs por comparación
deg_files <- list(
  Healthy_vs_PCa  = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/DEGs_by_CellType/DEGs_Healthy_vs_PCa_by_CellType.xlsx",
  Healthy_vs_CRPC = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/DEGs_by_CellType/DEGs_Healthy_vs_CRPC_by_CellType.xlsx",
  PCa_vs_CRPC     = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/DEGs_by_CellType/DEGs_PCa_vs_CRPC_by_CellType.xlsx"
)

# Directorio de salida
out_dir <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/TOP_DEGs_tables_allCellTypes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Parámetros de filtrado
padj_thr  <- 0.05
logfc_thr <- 0.25
n_top     <- 15   # número máximo de genes UP/DOWN a mostrar

################################################################################
# 2) Función para obtener tabla de TOP UP y DOWN (o mensaje si no hay DEGs)
################################################################################

get_top_deg_table <- function(df, n = 15, padj_thr = 0.05, logfc_thr = 0.25) {
  
  if (!"gene" %in% colnames(df)) df$gene <- rownames(df)
  
  required_cols <- c("gene", "p_val_adj", "avg_log2FC")
  if (!all(required_cols %in% colnames(df))) {
    return(tibble(Message = "No hay DEGs significantes (faltan columnas requeridas)"))
  }
  
  df_f <- df %>%
    filter(!is.na(p_val_adj), !is.na(avg_log2FC)) %>%
    filter(p_val_adj < padj_thr, abs(avg_log2FC) >= logfc_thr)
  
  if (nrow(df_f) == 0) {
    return(tibble(Message = "No hay DEGs significantes bajo los umbrales establecidos"))
  }
  
  # Selección de TOP UP y TOP DOWN
  top_up <- df_f %>%
    filter(avg_log2FC > 0) %>%
    arrange(p_val_adj, desc(avg_log2FC)) %>%
    slice_head(n = n) %>%
    mutate(Direction = "UP")
  
  top_down <- df_f %>%
    filter(avg_log2FC < 0) %>%
    arrange(p_val_adj, avg_log2FC) %>%
    slice_head(n = n) %>%
    mutate(Direction = "DOWN")
  
  bind_rows(top_up, top_down) %>%
    arrange(factor(Direction, levels = c("UP", "DOWN")), p_val_adj)
}

################################################################################
# 3) Bucle por comparación: se genera un Excel con TODAS las poblaciones
################################################################################

for (comp in names(deg_files)) {
  xlsx_path <- deg_files[[comp]]
  if (!file.exists(xlsx_path)) {
    next
  }
  
  sheets <- getSheetNames(xlsx_path)
  
  out_list <- list()
  
  for (ct in all_celltypes) {
    
    # Asegurar nombre compatible con Excel (máx 31 caracteres, sin símbolos conflictivos)
    safe_ct <- ct %>%
      str_replace_all("[\\[\\]\\*\\?\\:/\\\\]", "_") %>%
      str_replace_all("[^A-Za-z0-9_ ]", "_") %>%
      substr(1, 31)
    
    # Si NO hay hoja en el Excel: registrar motivo
    if (!ct %in% sheets) {
      out_list[[safe_ct]] <- tibble(
        Message = "No se generaron DEGs para esta población (la hoja no existe en el Excel; posible N insuficiente o ausencia en alguna condición)."
      )
      next
    }
    
    # Si hay hoja: calcular tabla de TOP genes (o mensaje si no hay significativos)
    df <- read.xlsx(xlsx_path, sheet = ct)
    
    top_tbl <- get_top_deg_table(
      df = df,
      n = n_top,
      padj_thr = padj_thr,
      logfc_thr = logfc_thr
    )
    
    out_list[[safe_ct]] <- top_tbl
  }
  
  # Evitar nombres duplicados si truncamiento genera colisiones
  names(out_list) <- make.unique(names(out_list), sep = "_")
  
  out_xlsx <- file.path(out_dir, paste0("TOP15_UP_DOWN_DEGs_", comp, "_allCellTypes.xlsx"))
  write.xlsx(out_list, file = out_xlsx, overwrite = TRUE)
}