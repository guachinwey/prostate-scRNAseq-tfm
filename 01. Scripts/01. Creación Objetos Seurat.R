################################################################################
# 1. Creación de objetos Seurat a partir de datos crudos.
# - PCA (matriz 10x-like en txt.gz)
# - CRPC (matrices densas en csv.gz)
# - Healthy (matrices tabuladas en txt)

# Output:
# - Objetos Seurat individuales guardados en formato RDS.

################################################################################

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyverse)


################################################################################
# 1) Población PCA: lectura desde matriz cruda (txt.gz)
# Función para cargar la muestra PCA desde un fichero crudo.
# Formato esperado:
# - Fila 1: barcodes de células (separadas por tabulador).
# - Columna 1: nombres de genes.
# - Columnas 2..:: matriz de cuentas (genes x células).

# Pasos principales:
# - Leer barcodes de la primera línez.
# - Leer el resto como tabla numérica.
# - Forzar nombres de gen únicos.
# - Submuestrear células (downsample) para limitar RAM.
# - Crear objeto Seurat y calcular % de genes mitocondriales.
################################################################################


load_pca_from_raw <- function(file_path, sample_name = "PCA", downsample_cells = 5000) {
  
  # Lectura de barcodes (primera línea)
  con <- gzfile(file_path, open = "rt")
  first_line <- readLines(con, n = 1)
  close(con)
  cell_barcodes <- strsplit(first_line, "\t")[[1]]
  
  # Lectura de la matriz de cuentas (resto de líneas)
  raw_tab <- read.table(
    file = gzfile(file_path),
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    skip = 1) # saltar la fila de barcodes
  
  # Columna 1 = nombres de genes
  gene_names <- as.character(raw_tab[[1]])
  
  # Columasn 2..N = matriz de cuentas
  count_matrix <- raw_tab[, -1, drop = FALSE]
  
  # Asegurarse de que todas las columnas son numéricas
  count_matrix <- as.data.frame(
    lapply(count_matrix, function(x) {
      if (is.factor(x)) x <- as.character(x)
      as.numeric(x)}),
    check.names = FALSE)
  
  message("Dimensiones brutas antes de renombrar: ",
          nrow(count_matrix), " genes x ", ncol(count_matrix), " células")
  
  # Se asignan nombres de gen únicos (Seurat no admite duplicados)
  rownames(count_matrix) <- make.unique(gene_names)
 
  # Asignar barcodes originales como nombres de columna
  colnames(count_matrix) <- cell_barcodes
  
  # Submuestreo para proteger RAM:
  if (!is.infinite(downsample_cells) && ncol(count_matrix) > downsample_cells) {
    set.seed(123)
    keep_cols <- sample(colnames(count_matrix), downsample_cells)
    count_matrix <- count_matrix[, keep_cols, drop = FALSE]}
  
  # Renombrar las columnas para añadir prefijo de muestra (PCA_1, PCA_2...)
  colnames(count_matrix) <- paste0(sample_name, "_", seq_len(ncol(count_matrix)))
  
  # Convertir a matriz dispersa (dgCMatrix) para ahorrar memoria
  count_matrix <- as.matrix(count_matrix)
  mode(count_matrix) <- "numeric"
  count_matrix <- as(count_matrix, "dgCMatrix")
  
  # Creación del objeto Seurat:
  seurat_object <- CreateSeuratObject(
    counts = count_matrix,
    project = "scRNAseq")
  
  # Se añaden metadatos al objeto Seurat:
  seurat_object$sample_id <- sample_name
  seurat_object$condition <- "PCa"
  
  # Cálculo del porcentaje mitocondrial:
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(
    seurat_object,
    pattern = "^MT-")
  
  message("Objeto final: ",
          nrow(seurat_object), " genes x ", ncol(seurat_object), " células")
  
  return(seurat_object)}


# Ruta de la matriz de conteo PCA:
pca_file <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/PCA/GSE141445_RAW/GSM4203181_data.raw.matrix.txt.gz"

# Ruta para guardar el objeto Seurat final:
output_dir_pca <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_PCA/"

# Crear objeto PCA (con submuestreo a 5000 células)
PCA <- load_pca_from_raw(
  file_path        = pca_file,
  sample_name      = "PCA",
  downsample_cells = 5000)  

# Guardar el objeto Seurat:
saveRDS(
  PCA,
  file = file.path(output_dir_pca, "PCA_fixed.rds"))

# Comprobación del objeto:
PCA
head(rownames(PCA), 20)
head(PCA@meta.data)
summary(PercentageFeatureSet(PCA, pattern = "^MT-"))
dim(PCA)


################################################################################
# 2) Muestras CRPC: lectura desde matrices densas (csv.gz)
################################################################################


# Aumentar tamaño del buffer de lectura de 'vroom' para mejorar velocidad
# y evitar errores de conexión al usar readr::read_csv con ficheros grandes.
Sys.setenv("VROOM_CONNECTION_SIZE" = 5 * 131072 * 10) 

# Función genérica para crear un objeto Seurat a partir de los csv de CRPC.
# Formato esperado:
# - Columna ...1: índice sin usar (se elimina).
# - Columna CLUSTER: identidad de clúster original (se guarda en metadata).
# - Resto de columnas: matriz de cuentas (genes x células).
create_seurat_object <- function(file_path, sample_name) {
  
  # Leer tabla CSV
  data <- read_csv(file_path)
  
  # Extraer columna de clúster para guardar como metadato
  cluster_info <- data$CLUSTER
  
  # Eliminar columnas de índice (...1) y columna CLUSTER
  count_matrix <- data %>% select(-c(...1, CLUSTER))
  count_matrix <- as.data.frame(count_matrix)
  
  # Asignar nombres de fila (gene) como índices consecutivos con prefijo de muestra
  # Cada fila corresponde a un gen distinto
  rownames(count_matrix) <- paste0(sample_name, "_", seq_len(nrow(count_matrix)))
  
  # Transponer matriz (genes en filas y células en columnas)
  matr <- t(count_matrix)
  matr <- as.matrix(matr)
  matr <- as(matr, "dgCMatrix")
  
  #Crear objeto Seurat
  seurat_object <- CreateSeuratObject(counts = matr, project = "scRNAseq")
  
  # Metadatos
  seurat_object$Cluster <- cluster_info
  seurat_object$sample_id <- sample_name
  seurat_object$condition <- "CRPC"
  
  return(seurat_object)}

# Rutas a las muestras CRPC
file_paths <- list(
  CRPC_1 = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/CRPC/GSM6428952_1778_JZ_HMP_04_IGO_10726_2_dense.csv.gz",
  CRPC_2 = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/CRPC/GSM6428954_1845_HMP-08_IGO_10837_19_dense.csv.gz",
  CRPC_3 = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/CRPC/GSM6428955_1968_HMP11_1_IGO_11247_3_dense.csv.gz")

seurat_objects <- list()

for (sample in names(file_paths)) {
  seurat_objects[[sample]] <- create_seurat_object(file_paths[[sample]], sample)
  saveRDS(seurat_objects[[sample]], 
          file = file.path("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_CRPC", paste0(sample, ".rds")))}

# Comprobaciones de cada objeto CRPC
CRPC_1
head(rownames(CRPC_1), 20)
head(CRPC_1@meta.data)
summary(PercentageFeatureSet(CRPC_1, pattern = "^MT-"))
dim(CRPC_1)
CRPC_2
head(rownames(CRPC_2), 20)
head(CRPC_2@meta.data)
summary(PercentageFeatureSet(CRPC_2, pattern = "^MT-"))
dim(CRPC_2)
CRPC_3
head(rownames(CRPC_3), 20)
head(CRPC_3@meta.data)
summary(PercentageFeatureSet(CRPC_3, pattern = "^MT-"))
dim(CRPC_3)


################################################################################
# 3) Muestras HEALTHY: matrices tabuladas (txt)
# Este dataset ya viene en formato tabulado estándar:
# - Genes en filas.
# - Células en columnas.
# - Primera columna = nombres de genes.
# - Encabezado con barcodes de célula.
################################################################################


load_healthy_sample <- function(file_path, sample_name) {
  matr <- read.table(
    file = file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE)
  
  # Primera columna = nombre de gen:
  gene_names <- matr[[1]]
  count_matrix <- matr[, -1]
  
  # Nombres de gen únicos
  gene_names_unique <- make.unique(gene_names)
  rownames(count_matrix) <- gene_names_unique
  
  # Renombrar columnas añadiendo prefijo de muestra
  colnames(count_matrix) <- paste0(sample_name, "_", seq_len(ncol(count_matrix)))
  
  # Conversión a matriz dispersa
  count_matrix <- as.matrix(count_matrix)
  mode(count_matrix) <- "numeric"
  count_matrix <- as(count_matrix, "dgCMatrix")
  
  # Crear objeto Seurat
  seurat_object <- CreateSeuratObject(counts = count_matrix, project = "scRNAseq")
  seurat_object$sample_id <- sample_name
  seurat_object$condition <- "Healthy"
  return(seurat_object)}

# Lista de rutas a los archivos no comprimidos:
file_paths <- list(
  HEALTHY_1 = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/HEALTH/GSM4556600_MF002_human_cleanraw_counts_matrix.txt",
  HEALTHY_2 = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/HEALTH/GSM4556601_MM033_human_clean_raw_counts_matrix.txt",
  HEALTHY_3 = "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/HEALTH/GSM4556602_MM037_clean_raw_counts_matrix.txt"
)

seurat_objects <- list()
output_dir <- "C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/SEURAT_HEALTH"

for (sample in names(file_paths)) {
  seurat_objects[[sample]] <- load_healthy_sample(
    file_path   = file_paths[[sample]],
    sample_name = sample)
  
  saveRDS(
    seurat_objects[[sample]], file = file.path(
      output_dir, paste0(sample, ".rds")))}

# Comprobaciones básicas de cada objeto Healthy
HEALTHY_1
head(rownames(HEALTHY_1), 20)
head(HEALTHY_1@meta.data)
summary(PercentageFeatureSet(HEALTHY_1, pattern = "^MT-"))
dim(HEALTHY_1)
HEALTHY_2
head(rownames(HEALTHY_2), 20)
head(HEALTHY_2@meta.data)
summary(PercentageFeatureSet(HEALTHY_2, pattern = "^MT-"))
dim(HEALTHY_2)
HEALTHY_3
head(rownames(HEALTHY_3), 20)
head(HEALTHY_3@meta.data)
summary(PercentageFeatureSet(HEALTHY_3, pattern = "^MT-"))
dim(HEALTHY_3)