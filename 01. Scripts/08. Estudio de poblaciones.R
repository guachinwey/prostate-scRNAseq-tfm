################################################################################
# 8. Análisis de distribución celular por condición
# - Carga del objeto anotado con tipos celulares.
# - Conteo y cálculo de porcentajes por combinación CellType–condición.
# - Exportación de tablas en formatos largo y ancho a Excel.
# - Visualización con barplots apilados:
#   (A) Porcentajes por tipo celular (¿de qué condición proviene cada tipo?).
#   (B) Porcentajes por condición (¿cómo se compone cada condición?).
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(tidyr)
library(forcats)


################################################################################
# 1) Carga del objeto anotado y selección de columna de anotación
################################################################################

merged_all_annotated <- readRDS(
  "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_annotated_manual.rds"
)

# Seleccionar la columna que contiene la anotación celular.
celltype_col <- "CellType_Manual"    # Anotación manual

# Verificación de columnas requeridas en meta.data.
stopifnot(celltype_col %in% colnames(merged_all_annotated@meta.data))
stopifnot("condition" %in% colnames(merged_all_annotated@meta.data))

# Extraer dataframe con condición y tipo celular, renombrando columna como "CellType".
meta_df <- merged_all_annotated@meta.data %>%
  dplyr::select(condition, !!sym(celltype_col)) %>%
  dplyr::rename(CellType = !!sym(celltype_col))

# Reemplazo de valores faltantes por "Unknown".
meta_df$CellType[is.na(meta_df$CellType)] <- "Unknown"


################################################################################
# 2) Tabla resumen (formato largo): CellType–condición con n y %
################################################################################

# Objetivo: para cada tipo celular, calcular el % de células que provienen de cada condición.
celltype_condition_summary <- meta_df %>%
  dplyr::group_by(CellType, condition) %>%
  dplyr::summarise(n = n(), .groups = "drop_last") %>%
  dplyr::mutate(pct = n / sum(n) * 100) %>%
  dplyr::ungroup()

celltype_condition_summary


################################################################################
# 3) Tabla resumen en formato ancho
################################################################################

# Pivotear la tabla para obtener columnas separadas por condición (n y %).
celltype_condition_wide <- celltype_condition_summary %>%
  tidyr::pivot_wider(
    id_cols     = CellType,
    names_from  = condition,
    values_from = c(n, pct),
    values_fill = 0
  ) %>%
  dplyr::mutate(total_cells = n_CRPC + n_Healthy + n_PCa) %>%
  dplyr::arrange(desc(total_cells))

celltype_condition_wide


################################################################################
# 4) Exportación de tablas a archivo Excel
################################################################################

out_xlsx <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/CellType_condition_summary.xlsx"

openxlsx::write.xlsx(
  list(
    "long_format" = celltype_condition_summary,
    "wide_format" = celltype_condition_wide
  ),
  file = out_xlsx
)


################################################################################
# 5A) Barplot apilado: % de condición dentro de cada tipo celular
################################################################################

# Esta visualización muestra la distribución de condiciones dentro de cada tipo celular.

# Ordenar los tipos celulares por abundancia total para mejorar la visualización.
celltype_order <- celltype_condition_wide$CellType

celltype_condition_summary <- celltype_condition_summary %>%
  dplyr::mutate(CellType = factor(CellType, levels = celltype_order))

p_bar_celltype_within <- ggplot(celltype_condition_summary,
                                aes(x = CellType, y = pct, fill = condition)) +
  geom_col(color = "black", width = 0.85) +
  scale_y_continuous(expand = c(0, 0), name = "Porcentaje dentro del tipo celular (%)") +
  xlab("Tipo celular") +
  ggtitle("Distribución por condición dentro de cada tipo celular") +
  theme_bw(base_size = 12) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_bar_celltype_within)

# Guardar gráfico en PNG y PDF
out_png <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/CellType_condition_barplot_withinCellType.png"
out_pdf <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/CellType_condition_barplot_withinCellType.pdf"

ggsave(out_png, plot = p_bar_celltype_within, width = 10, height = 5, dpi = 300)
ggsave(out_pdf, plot = p_bar_celltype_within, width = 10, height = 5)


################################################################################
# 5B) Barplot apilado: composición celular por condición
################################################################################

# Cómo se compone cada condición en términos de tipos celulares.

condition_celltype_summary <- meta_df %>%
  dplyr::group_by(condition, CellType) %>%
  dplyr::summarise(n = n(), .groups = "drop_last") %>%
  dplyr::mutate(pct = n / sum(n) * 100) %>%
  dplyr::ungroup()

# Ordenar tipos celulares por abundancia global (opcional)
global_order <- meta_df %>%
  dplyr::count(CellType, sort = TRUE) %>%
  dplyr::pull(CellType)

condition_celltype_summary <- condition_celltype_summary %>%
  dplyr::mutate(CellType = factor(CellType, levels = global_order))

p_bar_condition <- ggplot(condition_celltype_summary,
                          aes(x = condition, y = pct, fill = CellType)) +
  geom_col(color = "black", width = 0.75) +
  scale_y_continuous(expand = c(0, 0), name = "Porcentaje dentro de la condición (%)") +
  xlab("Condición") +
  ggtitle("Composición celular por condición") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_bar_condition)

# Guardar gráfico en PNG y PDF
out_png2 <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/Condition_celltype_composition_barplot.png"
out_pdf2 <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/Condition_celltype_composition_barplot.pdf"

ggsave(out_png2, plot = p_bar_condition, width = 8, height = 5, dpi = 300)
ggsave(out_pdf2, plot = p_bar_condition, width = 8, height = 5)