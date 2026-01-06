################################################################################
# 2. Control de calidad (QC) por muestra.
# - Cálculo de métricas (nFeature, nCount, %MT, %Ribo, MALAT1)
# - Fijación de umbrales adaptativos (percentiles).
# - Filtrado de células de baja calidad.
# - Tabla resumen de QC por muestra.
# - Boxplots e histogramas de distribución.
################################################################################


library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(patchwork)


################################################################################
# 1) Función para calcular métricas de QC por muestra
# Función: compete_qc_metrics()

# Objetivo:
# - Calcular y añadir al objeto Seurat las métricas básicas de QC por célula.

# Formato de entrada:
# - seu: objeto Seurat crudo (sin filtrar).
# - sample_name: nombre de la muestra, usado para etiquetar la tabla.

# Métricas calculadas:
# - percent.mt: porcentaja de genes mitocondriales (^MT-)
# - percent.ribo: porcentaje de genes ribosomales (RPL*/RPS*)
# - MALAT1: conteos crudos del gen MALAT1 por célula.

# Salida:
# - lista con:
# - $seu: el objeto Seurat actualizado con columnas percent.mt y percent.ribo.
# - $qc: un data.frame con las métricas de QC por célula.
################################################################################


compute_qc_metrics <- function(seu, sample_name) {
  
  # Porcentaje mitocondrial (solo si no existe en metadata):
  if (!"percent.mt" %in% colnames(seu@meta.data)) {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  }
  
  # Porcentaje de genes ribosomales
  # Si no se detectan, el resultado es NA
  ribo_genes <- grep("^RPL|^RPS", rownames(seu), value = TRUE)
  if (length(ribo_genes) > 0) {
    seu[["percent.ribo"]] <- PercentageFeatureSet(seu, features = ribo_genes)
  } else {
    seu$percent.ribo <- NA_real_
  }
  
  # Expresión de MALAT1 (indicador de fracción nuclear / calidad)
  # Si no está presente, se informa con NA.
  if ("MALAT1" %in% rownames(seu)) {
    mal_counts <- as.numeric(
      Matrix::colSums(
        GetAssayData(seu, slot = "counts")["MALAT1", , drop = FALSE]
      )
    )
  } else {
    mal_counts <- rep(NA_real_, ncol(seu))
  }
  
  # Tabla de QC por célula con todas las métricas
  qc_df <- data.frame(
    cell          = colnames(seu),
    sample        = sample_name,
    nCount_RNA    = seu$nCount_RNA,
    nFeature_RNA  = seu$nFeature_RNA,
    percent.mt    = seu$percent.mt,
    percent.ribo  = seu$percent.ribo,
    MALAT1_expr   = mal_counts,
    stringsAsFactors = FALSE
  )
  
  return(list(seu = seu, qc = qc_df))
}


################################################################################
# 2) Rutas de entrada (objetos Seurat crudos)
################################################################################


# Lista con rutas a los objetos Seurat generados en el Script 1.
# Cada entranda corresponde a una muestra antes de QC.

input_paths <- list(
  HEALTHY_1 = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_HEALTH/HEALTHY_1.rds",
  HEALTHY_2 = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_HEALTH/HEALTHY_2.rds",
  HEALTHY_3 = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_HEALTH/HEALTHY_3.rds",
  PCA       = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_PCA/PCA.rds",
  CRPC_1    = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_CRPC/CRPC_1.rds",
  CRPC_2    = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_CRPC/CRPC_2.rds",
  CRPC_3    = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_CRPC/CRPC_3.rds"
)

# Carpeta de salida:
# - output_dir_qc: objetos Seurat filtrados por QC.
# - output_dir_figs: boxplots de métricas de QC por muestra.

output_dir_qc   <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_QC"
output_dir_figs <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/QC_figures"

dir.create(output_dir_qc,   showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_figs, showWarnings = FALSE, recursive = TRUE)

# Lista para guardar los objetos QC
seurat_qc_objects <- list()

# Lista para guardar las tablas de QC por muestra
qc_tables <- list()


################################################################################
# 3) Bucle principal: QC data-driven por muestra
# En este bucle:
# - Se carga cada objeto Seurat crudo.
# - Se calculan métricas de QC.
# - Se generan boxplots recortados a [p1, p99].
# - Se calculan umbrales adaptativos por percentiles.
# - Se filtran células de baja calidad.
# - Se guardan los objetos filtrados y se reporta el número de células retenidas.
################################################################################


for (sample in names(input_paths)) {
  message("Procesando QC de: ", sample)
  
  # Cargar objeto Seurat crudo
  seu_raw <- readRDS(input_paths[[sample]])
  
  # Calcular métricas de QC (mito, ribo, MALAT1, etc.)
  qc_list <- compute_qc_metrics(seu_raw, sample_name = sample)
  seu_raw <- qc_list$seu
  qc_df   <- qc_list$qc
  qc_tables[[sample]] <- qc_df
  
  # Pasar de formato ancho a largo para facilitar la generación de boxplots
  qc_long <- qc_df |>
    select(nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, MALAT1_expr) |>
    pivot_longer(cols = everything(),
                 names_to = "metric",
                 values_to = "value")
  
  # Recortar valores al rango [p1, p99] para evitar colas extremas
  qc_long <- qc_long |>
    group_by(metric) |>
    mutate(
      q1  = quantile(value, 0.01, na.rm = TRUE),
      q99 = quantile(value, 0.99, na.rm = TRUE),
      value_trunc = pmin(pmax(value, q1), q99)
    ) |>
    ungroup()
  
  # Boxplot de métricas QC por muestra (nUMI, nFeatura, %mito, %ribo, MALAT1)
  p_box <- ggplot(qc_long, aes(x = metric, y = value_trunc)) +
    geom_boxplot(outlier.shape = NA, fill = "lightblue") +
    theme_bw(base_size = 12) +
    coord_flip() +
    labs(
      title = paste("Métricas de QC -", sample),
      x     = "Métrica",
      y     = "Valor (recortado al percentil 1–99)"
    )

  # Guardar boxplots como PDF por muestra  
  ggsave(
    filename = file.path(output_dir_figs, paste0(sample, "_QC_boxplots.pdf")),
    plot     = p_box,
    width    = 7,
    height   = 5
  )
  
  # Definir umbrales (percentiles) por muestra:
  
  # nFeature_RNA: entre p1 y p99, acotando a [200, 6000]
  feat_quant   <- quantile(qc_df$nFeature_RNA, c(0.01, 0.99), na.rm = TRUE)
  min_features <- max(200, floor(feat_quant[1]))
  max_features <- min(6000, ceiling(feat_quant[2]))
  
  # percent.mt: usar p99, pero no más estricto que 15% (consistente con literatura)
  mito_99 <- quantile(qc_df$percent.mt, 0.99, na.rm = TRUE)
  if (is.na(mito_99)) {
    max_percent_mt <- 15
  } else {
    max_percent_mt <- min(15, mito_99)
  }
  
  # MALAT1: solo se registra el p99 (no se usa como filtro estricto)
  malat1_99 <- quantile(qc_df$MALAT1_expr, 0.99, na.rm = TRUE)
  
  # Reporte de los umbrales calculados para esa muestra
  message("Umbrales propuestos para ", sample, " :",
          "\n  nFeature_RNA entre ", min_features, " y ", max_features,
          "\n  percent.mt <= ", round(max_percent_mt, 2),
          "\n  MALAT1_expr p99 = ", ifelse(is.na(malat1_99), "NA", round(malat1_99, 2)),
          "\n")
  
  # Aplicar filtrado con esos umbrales
  seu_qc <- subset(
    seu_raw,
    subset = nFeature_RNA >= min_features &
      nFeature_RNA <= max_features &
      (is.na(percent.mt) | percent.mt <= max_percent_mt)
  )
  
  # Guardar objeto filtrado en la lista y en disco y resumen
  seurat_qc_objects[[sample]] <- seu_qc
  
  saveRDS(
    seu_qc,
    file = file.path(
      output_dir_qc,
      paste0(sample, "_QC.rds")
    )
  )
  
  message(
    "Células en ", sample,
    " | antes: ", ncol(seu_raw),
    " | después QC: ", ncol(seu_qc), "\n"
  )
}


################################################################################
# 4) Tabla resumen final de QC por muestra.
# Se construye una tabla resumen con:
# - Número de células antes y después de QC.
# - Porcentaje retenido.
# - Percentiles p1/p99 de nFeature_RNA.
# - p99 de percent.mt, MALAT1 y percent.ribo
# La tabla se exporta como CSV para documentar el proceso de filtrado.
################################################################################


# Lista donde almacenaremos los datos
qc_summary <- data.frame()

for (sample in names(seurat_qc_objects)) {
  
  # Objeto crudo (antes de QC) y filtrado (después de QC)
  seu_raw <- readRDS(input_paths[[sample]])
  seu_qc  <- seurat_qc_objects[[sample]]
  
  # Métricas por célula del objeto sin filtrar (para cálculo de percentiles)
  nfeature_raw <- seu_raw$nFeature_RNA
  percent_mt_raw <- PercentageFeatureSet(seu_raw, pattern = "^MT-")
  malat1_raw <- GetAssayData(seu_raw)["MALAT1", ]
  ribo_raw <- PercentageFeatureSet(seu_raw, pattern = "^RPL|^RPS")
  
  # Percentiles principales apra resumir la distribución de cada métrica
  p1_nFeature  <- quantile(nfeature_raw, 0.01, na.rm = TRUE)
  p99_nFeature <- quantile(nfeature_raw, 0.99, na.rm = TRUE)
  p99_mt       <- quantile(percent_mt_raw, 0.99, na.rm = TRUE)
  p99_malat1   <- quantile(malat1_raw, 0.99, na.rm = TRUE)
  p99_ribo     <- quantile(ribo_raw, 0.99, na.rm = TRUE)
  
  # Células antes/después del QC y porcentaje retenido
  before  <- ncol(seu_raw)
  after   <- ncol(seu_qc)
  kept_pct <- round((after / before) * 100, 2)
  
  # Añadir fila resumen para esa muestra
  qc_summary <- rbind(qc_summary, data.frame(
    Sample = sample,
    Cells_before = before,
    Cells_after = after,
    Percent_kept = kept_pct,
    nFeature_p1 = round(p1_nFeature, 2),
    nFeature_p99 = round(p99_nFeature, 2),
    max_percent_mt_used = min(15, round(p99_mt, 2)),
    MALAT1_p99 = round(as.numeric(p99_malat1), 2),
    Ribo_p99 = round(as.numeric(p99_ribo), 2)
  ))
}

# Guardar tabla resumen de QC como CSV
write_csv(qc_summary,
          "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/QC_table.csv")

qc_summary


################################################################################
# 5) Histograma.
# Se generan histogramas con densidad:
# - Distribución de nFeature_RNA.
# - Distribución de percent.mt
# Por cada muestra tras QC, se guardan en PDFs independientes.
################################################################################


# Carpeta donde guardar los PDFs de histogramas
out_dir <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/Results/QC_histograms"

# Función: genera histogramas + densidad para una muestra

# Entrada:
# - seu: objeto Seurat ya filtrado por QC.
# - sample_name: nombre de la muestra (se usa en los títulos)
# - out_dir: carpeta de salida para los PDFs.

# Salida:
# - Un PDF por muestra con dos paneles:
# 1) Histograma + densidad de nFeature_RNA.
# 2) Histograma + densidad de percent.mt.

plot_qc_histograms <- function(seu, sample_name, out_dir) {
  
  # Data.frame con las métricas que se van a representar
  df <- data.frame(
    nFeature_RNA = seu$nFeature_RNA,
    percent.mt   = seu$percent.mt
  )
  
  # Histograma + densidad de nFeature_RNA
  p_features <- ggplot(df, aes(x = nFeature_RNA)) +
    geom_histogram(aes(y = ..density..),
                   bins = 60,
                   fill = "grey80",
                   color = "black") +
    geom_density(color = "red", linewidth = 1) +
    labs(
      title = paste0(sample_name, " - Distribución nFeature_RNA"),
      x = "Número de genes detectados por célula",
      y = "Densidad"
    ) +
    theme_bw()
  
  # Histograma + densidad de percent.mt
  p_mt <- ggplot(df, aes(x = percent.mt)) +
    geom_histogram(aes(y = ..density..),
                   bins = 60,
                   fill = "grey80",
                   color = "black") +
    geom_density(color = "blue", linewidth = 1) +
    labs(
      title = paste0(sample_name, " - Distribución percent.mt"),
      x = "% genes mitocondriales por célula",
      y = "Densidad"
    ) +
    theme_bw()
  
  # Guardar ambos en un único PDF (dos paneles)
  pdf_file <- file.path(out_dir, paste0(sample_name, "_QC_histograms.pdf"))
  pdf(pdf_file, width = 10, height = 5)
  print(p_features + p_mt)
  dev.off()
  
  message("✅ Guardado: ", pdf_file)
}

# Lista de muestras QC a procesar
samples_qc <- c(
  "HEALTHY_1_QC",
  "HEALTHY_2_QC",
  "HEALTHY_3_QC",
  "PCA_QC",
  "CRPC_1_QC",
  "CRPC_2_QC",
  "CRPC_3_QC"
)

# Bucle: cargar cada objeto QC y generar sus histogramas
for (fname in samples_qc) {
  sample_name <- sub("_QC$", "", fname)  # quitar el _QC para el título
  
  seu_qc <- readRDS(file.path(qc_dir, paste0(sample_name, "_QC.rds")))
  
  plot_qc_histograms(
    seu = seu_qc,
    sample_name = sample_name,
    out_dir = out_dir
  )
}