# Análisis computacional de datos de scRNA-seq en cáncer de próstata

Este repositorio contiene el código, los resultados y la documentación asociados al Trabajo Final de Máster (TFM) titulado:

**“Análisis computacional de datos de scRNA-seq en cáncer de próstata: integración, caracterización y comparación de poblaciones celulares en tejido sano, PCa y CRPC”**

El objetivo del trabajo es el análisis de datos de transcriptómica a nivel de célula única (scRNA-seq) para la identificación, anotación y comparación de poblaciones celulares en muestras de próstata sana, cáncer de próstata (PCa) y cáncer de próstata resistente a la castración (CRPC).

---

## Estructura del repositorio

El repositorio está organizado de la siguiente manera:

```
prostate-scRNAseq-tfm/
├── Scripts/
├── Resultados/
├── Memoria/
└── README.md
```

---

## 📁 Scripts

Esta carpeta contiene los scripts en R utilizados para el análisis de los datos de scRNA-seq.  
Los scripts están numerados para reflejar el orden lógico del flujo de trabajo (creación de objetos Seurat, control de calidad, integración, clustering, anotación celular y análisis de expresión diferencial).

*Para una descripción detallada del flujo de trabajo y de los parámetros utilizados, consulte el apartado [**Memoria**](./Memoria).*

---

## 📁 Resultados

Incluye las figuras, tablas y salidas finales del análisis, organizadas por bloques temáticos (control de calidad, clustering, genes característicos, anotación poblacional, estudio poblacional y DEGs).

Estos resultados corresponden a los análisis descritos y discutidos en la memoria del TFM.

---

## 📁 Memoria

Contiene el documento completo del Trabajo Final de Máster en formato PDF, donde se describen en detalle:
- los conjuntos de datos utilizados
- la metodología aplicada
- los parámetros empleados
- los resultados obtenidos y su interpretación

---

## 🧬 Disponibilidad de los datos

Los datos originales de scRNA-seq utilizados en este estudio no se incluyen en este repositorio debido a su gran tamaño y a las limitaciones de almacenamiento de GitHub.

Todos los conjuntos de datos empleados son de acceso público y pueden obtenerse a través de los repositorios originales (por ejemplo, GEO). Las referencias completas y los identificadores de acceso se encuentran debidamente documentados en la memoria del Trabajo Final de Máster.

Siguiendo el flujo de trabajo proporcionado en este repositorio, es posible reproducir los resultados descritos en el estudio en la medida en que lo permite la disponibilidad de los datos.

---

## 🔁 Reproducibilidad

Este repositorio proporciona los scripts de análisis, los resultados finales y la documentación metodológica necesaria con el objetivo de garantizar la transparencia y reproducibilidad del trabajo.

---

## 👤 Autoría

Trabajo realizado por **Silvia Arroitajauregui Avilés**  
Máster en Bioinformática y Bioestadística
