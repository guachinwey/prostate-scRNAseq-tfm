# Scripts de análisis.

Esta carpeta contiene los scripts en R utilizados para el análisis computacional de datos de scRNA-seq en cáncer de próstata, realizados en el marco del Trabajo Final de Máster.

Los scripts están numerados para indicar el orden lógico de ejecución del pipeline de análisis.

## Estructura del pipeline

- **01. Creación Objetos Seurat**  
  Importación de los datos originales y creación de los objetos Seurat.

- **02. Control de Calidad (QC)**  
  Filtrado de células y genes basado en métricas de calidad.

- **03. Integración, reducción dimensional y clustering**  
  Integración de muestras, reducción dimensional (PCA, UMAP) y clustering.

- **03.1. Decisión parámetros**  
  Evaluación y selección de parámetros para clustering y resolución.

- **04. TOP y genes característicos**  
  Identificación de genes característicos y marcadores por clúster.

- **05. Anotación automática de poblaciones celulares**  
  Anotación celular basada en referencias y métodos automáticos.

- **06. Estudio de poblaciones**  
  Análisis exploratorio de las poblaciones celulares identificadas.

- **07. Anotación manual de poblaciones celulares**  
  Curación manual y refinamiento de la anotación celular.

- **07.1. Luminal MALAT1**  
  Análisis específico de la población luminal MALAT1.

- **08. Estudio de poblaciones**  
  Análisis adicional de poblaciones celulares tras la anotación final.

- **09. DEGs por poblaciones**  
  Análisis de genes diferencialmente expresados entre poblaciones celulares.

- **10. Poblaciones a comparar**  
  Definición de comparaciones entre poblaciones celulares.

- **11. Volcano plots entre poblaciones**  
  Generación de volcano plots para las comparaciones de expresión diferencial.

- **12. Heatmaps de DEGs**  
  Visualización de genes diferencialmente expresados mediante heatmaps.

- **13. Tablas TOP DEGs (UP-DOWN)**  
  Generación de tablas con los genes diferencialmente expresados más relevantes.

  
*Para una descripción detallada del flujo de trabajo y de los parámetros utilizados, consulte el apartado [**Memoria**](../Memoria).*
