################################################################################
# 6. Estudio de poblaciones: distribución clúster–condición.
# - Carga del objeto anotado automáticamente (scType).
# - Cálculo de nº de células y porcentaje por clúster–condición (formato largo).
# - Generación de tabla “ancha” con n_ y pct_ por condición.
# - Exportación de ambas tablas a Excel.
# - Barplot apilado de composición de cada clúster por condición (Healthy, PCa, CRPC).
################################################################################

# Estudio de poblaciones

library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)  


################################################################################
# 1) Cargar objeto anotado y preparar tabla clúster–condición 
################################################################################


merged_all_annotated <- readRDS("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_annotated.rds")

# Se construye una tabla resumen clúster–condición:
# - seurat_clusters: identidad de clúster
# - condition: estado biológico (Healthy / PCa / CRPC).
# - n: nº de células por combinación clúster–condición.
# - pct: porcentaje que representan dentro de cada clúster.
cluster_condition_summary <- merged_all_annotated@meta.data %>%
  dplyr::select(seurat_clusters, condition) %>%          # nos quedamos solo con lo que necesitamos
  dplyr::group_by(seurat_clusters, condition) %>%        # agrupamos por clúster y condición
  dplyr::summarise(n = n(), .groups = "drop_last") %>%   # contamos nº de células (n)
  dplyr::mutate(
    pct = n / sum(n) * 100                               # porcentaje dentro de cada clúster
  ) %>%
  dplyr::ungroup()

# Resumen
cluster_condition_summary


################################################################################
# 2) Tabla “ancha” (n_ y pct_ para cada condición)
################################################################################


cluster_condition_wide <- cluster_condition_summary %>%
  # pasamos a formato ancho creando columnas n_Healthy, n_PCa, ... y pct_Healthy, pct_PCa, ...
  tidyr::pivot_wider(
    id_cols      = seurat_clusters,
    names_from   = condition,
    values_from  = c(n, pct),
    values_fill  = 0
  ) %>%
  # opcional: ordenar por número total de células
  dplyr::mutate(
    total_cells = n_CRPC + n_Healthy + n_PCa
  ) %>%
  dplyr::arrange(desc(total_cells))

# Vista rápida
cluster_condition_wide


################################################################################
# 3) Exportar tablas a Excel (formato largo y ancho)
################################################################################


out_xlsx <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/Cluster_condition_summary.xlsx"

openxlsx::write.xlsx(
  list(
    "long_format"  = cluster_condition_summary,  # tabla larga (una fila por cluster-condición)
    "wide_format"  = cluster_condition_wide      # tabla ancha (n_ y pct_ por condición)
  ),
  file = out_xlsx
)


################################################################################
# 4) Barplot apilado: composición de cada clúster por condición
################################################################################


# Aseguramos que los clústeres estén ordenados de forma razonable
cluster_condition_summary <- cluster_condition_summary %>%
  dplyr::mutate(
    seurat_clusters = forcats::fct_inorder(seurat_clusters)
  )

# Gráfico: barras apiladas por clúster, altura = porcentaje, color = condición.
# - Permite ver visualmente qué condiciones dominan en cada clúster.
p_bar <- ggplot(cluster_condition_summary,
                aes(x = seurat_clusters,
                    y = pct,
                    fill = condition)) +
  geom_col(color = "black", width = 0.8) +
  scale_y_continuous(expand = c(0, 0), name = "Porcentaje de células (%)") +
  xlab("Clúster") +
  ggtitle("Composición de cada clúster por condición") +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Mostrar en pantalla
print(p_bar)

# Guardar en archivo (PNG y PDF)
out_png <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/Cluster_condition_barplot.png"
out_pdf <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/Cluster_condition_barplot.pdf"

ggsave(out_png, plot = p_bar, width = 8, height = 5, dpi = 300)
ggsave(out_pdf, plot = p_bar, width = 8, height = 5)