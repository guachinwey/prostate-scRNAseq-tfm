################################################################################
# 10. Selección de poblaciones comparables entre condiciones (CellType_Manual)
# - Carga del objeto anotado manualmente.
# - Conteo de células por combinación CellType–condición.
# - Identificación de tipos celulares presentes con N suficiente en las 3 condiciones.
################################################################################

library(Seurat)
library(dplyr)
library(tidyr)

# Cargar objeto anotado 
obj <- readRDS("C:/Users/silvi/Desktop/MSc Bioinformatica y Bioestadistica/202526-1/Trabajo final de máster/DATA/MERGED/merged_all_annotated_manual.rds")

# Verificación de columnas necesarias
stopifnot("CellType_Manual" %in% colnames(obj@meta.data))
stopifnot("condition" %in% colnames(obj@meta.data))

# Umbral mínimo de células por grupo 
min_cells <- 50  

# Tabla de conteos por CellType_Manual y condición
ct_counts <- obj@meta.data %>%
  select(CellType_Manual, condition) %>%
  count(CellType_Manual, condition, name = "n") %>%
  pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
  mutate(
    # Evaluar si hay suficientes células en las 3 condiciones
    ok_all3 = (Healthy >= min_cells) & (PCa >= min_cells) & (CRPC >= min_cells)
  ) %>%
  arrange(desc(ok_all3), desc(PCa + CRPC + Healthy))  # ordena por viabilidad y abundancia total

ct_counts

# Extraer lista de tipos celulares comparables entre condiciones (presencia y N ≥ min_cells)
ct_ok <- ct_counts %>% filter(ok_all3) %>% pull(CellType_Manual)
ct_ok