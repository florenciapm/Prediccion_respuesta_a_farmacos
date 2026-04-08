# Predicción de respuesta a Selumetinib mediante expresión génica

Modelo de machine learning para predecir la sensibilidad de líneas celulares tumorales
a Selumetinib (inhibidor de MEK1/2) a partir de datos de expresión génica (RNA-seq).

## Datos

Los datos no están incluidos en el repositorio debido a su tamaño. Descárgalos en:
- GDSC2 IC50: https://www.cancerrxgene.org/downloads
- RNA-seq TPM: https://cellmodelpassports.sanger.ac.uk

Archivos necesarios:
- `GDSC2_fitted_dose_response_27Oct23.xlsx`
- `rnaseq_merged_rsem_tpm_20260323.csv`
- `screened_compounds_rel_8.5.csv`

## Metodología

### 1. Preprocesado

#### 1.1 Selección del fármaco y variable target

El dataset GDSC2 contiene datos de respuesta para cientos de fármacos. Se seleccionó
**Selumetinib**, inhibidor de MEK1/2 en la vía RAS/MAPK, por su relevancia clínica en
oncología de precisión y su amplia cobertura de líneas celulares.

Durante la exploración se identificó que Selumetinib aparece con dos DRUG_IDs distintos
(1062 y 1736), correspondientes a dos lotes o versiones del compuesto incorporadas en
distintas fases del proyecto GDSC. Se seleccionó el **DRUG_ID 1062** por ser el original
y tener mayor cobertura (949 líneas celulares vs 717 del 1736). Cada DRUG_ID no contiene
duplicados internos, confirmado explícitamente.

Como variable target se utilizó **LN_IC50** (logaritmo natural del IC50) por tres razones:
- Distribución aproximadamente normal, adecuada para regresión
- Medida directa e interpretable de la potencia del fármaco
- Estándar en la literatura para este tipo de análisis con GDSC

Se descartó AUC por su distribución fuertemente sesgada hacia 1 (la mayoría de líneas
son resistentes) y Z_SCORE por ser más útil en comparaciones entre fármacos distintos.

#### 1.2 Reestructuración del RNA-seq

El archivo TPM original tiene un formato no estándar para ML:
- 3 filas de metadatos (SIDM IDs, nombres de modelos, fuente de datos)
- Filas = genes, columnas = líneas celulares

Se procesó para obtener un conjunto tidy, con el formato necesario para utilizar la libreria scikit-learn:
- Eliminación de filas de metadatos
- Transposición de la matriz: **filas = líneas celulares, columnas = genes**
- Asignación de SIDM IDs como índice

#### 1.3 Resolución de gen duplicado

Se detectó el gen **VARS1** duplicado en el dataset, con valores idénticos en ambas
columnas, siendo un error de copia en los datos originales. Se resolvió promediando
ambas columnas, resultando en una única columna para VARS1.

#### 1.4 Cruce de datasets

Se cruzaron ambos datasets por SIDM (Sanger Model ID), identificador único de cada
línea celular. De las 949 líneas con datos de Selumetinib y las 1898 con datos de
expresión génica, se obtuvieron **925 líneas celulares comunes**. Se verificó
explícitamente la alineación de índices antes de proceder al modelado.

#### 1.5 Tratamiento de valores ausentes

Se identificaron **4728 genes con valores ausentes** en el dataset de expresión,
distribuidos en tres grupos con patrones sistemáticos:

| Grupo | Nº genes | NaNs por gen | Causa |
|-------|----------|--------------|-------|
| Sin medición en ninguna línea | 247 | 925 | Gen no incluido en ningún panel |
| Sin medición en líneas Broad | 4459 | 536 | Gen no medido por el laboratorio Broad |
| Sin medición en líneas Sanger | 22 | 389 | Gen no medido por el laboratorio Sanger |

El patrón sistemático por laboratorio (Broad: 536 líneas, Sanger: 389 líneas) indica
que son diferencias de panel entre plataformas, no errores aleatorios. Por tanto,
imputar estos valores no sería apropiado ya que se estarían inventando datos para hasta
el 58% de las muestras en ciertos genes. Se optó por **eliminar los 4728 genes
afectados**, quedando 36416 genes con mediciones completas en todas las muestras.

Los genes eliminados se guardaron en `genes_eliminados.csv` para referencia posterior
y verificación de que no incluyen genes de relevancia biológica conocida para la vía
RAS/MAPK.

#### 1.6 Filtro de expresión mínima

Se eliminaron genes con expresión prácticamente nula, ya que no aportan información
predictiva. El criterio utilizado fue: **TPM > 0.8 en al menos el 10% de las muestras**
(≥93 líneas celulares). Este umbral es más robusto que filtrar por media, ya que un gen
con media alta podría tenerla por unos pocos outliers sin expresarse realmente en la
mayoría de muestras. Tras este filtro: **19223 genes**.

#### 1.7 Transformación logarítmica

Se aplicó la transformación **log2(TPM + 1)** estándar en análisis de RNA-seq. Los
valores TPM tienen una distribución muy sesgada hacia 0, con genes altamente expresados
dominando el espacio de features. La transformación logarítmica comprime estos rangos
y estabiliza la varianza, haciendo los datos más adecuados para ML. El +1 evita
log(0) para genes con expresión 0.

#### 1.8 Filtro de varianza

Se eliminaron genes con baja varianza entre líneas celulares, ya que si un gen tiene
expresión similar en todas las muestras no puede explicar diferencias en respuesta al
fármaco. Se aplicó un umbral de **Q50 sobre los datos log-transformados**, eliminando
el 50% de genes con menor varianza. Dataset final: **925 muestras × 9611 genes**.

**Resumen del pipeline de preprocesado:**

| Paso | Genes resultantes |
|------|-------------------|
| Dataset original | 41144 |
| Tras eliminar genes con NaNs | 36416 |
| Tras filtro de expresión (TPM > 0.8 en ≥10% muestras) | 19223 |
| Tras transformación log2(TPM+1) y filtro de varianza (Q50) | **9611** |

### 2. Modelos evaluados

Se evaluaron tres modelos de regresión para predecir LN_IC50 a partir del perfil de
expresión génica. Los datos se dividieron en train (80%, 740 muestras) y test (20%,
185 muestras). Adicionalmente se realizó validación cruzada de 5 folds para una
estimación más robusta del rendimiento, menos dependiente de la partición train/test.

#### 2.1 Random Forest

Se entrenó un Random Forest como modelo baseline por su robustez con datos de alta
dimensionalidad y su capacidad de capturar relaciones no lineales entre genes, 
características habituales en datos de expresión génica.

La validación cruzada (R²=0.368 ± 0.069) mostró variabilidad notable entre folds
(0.281 a 0.449), indicando cierta sensibilidad a la partición de los datos, esperable
con 925 muestras y alta dimensionalidad.

Se realizó tuning de hiperparámetros mediante RandomizedSearchCV (20 iteraciones, 5
folds), obteniendo como configuración óptima: `n_estimators=200`, `max_depth=20`,
`min_samples_leaf=2`, `max_features=0.3`. La mejora respecto al baseline fue marginal
(R²=0.400 vs 0.398), con un R² CV inferior al baseline (0.344 vs 0.368), lo que
indica que el modelo había alcanzado su techo predictivo con los datos disponibles.

#### 2.2 XGBoost

Se evaluó XGBoost como alternativa de gradient boosting. En el split test obtuvo
resultados similares al Random Forest (R²=0.394), sin embargo la validación cruzada
reveló un rendimiento inferior y menos estable (R²=0.307 ± 0.060), con folds que
llegaron a 0.232. Dado este comportamiento en CV, se descartó como modelo final.

#### 2.3 ElasticNet

Se evaluó ElasticNet como alternativa lineal, dado que los modelos lineales con
regularización L1+L2 son competitivos en datos de expresión génica cuando las
relaciones gen-respuesta son aproximadamente lineales. Requiere escalado previo de
features (StandardScaler). Los hiperparámetros óptimos fueron `alpha=1.0`,
`l1_ratio=0.1`.

A diferencia del Random Forest, ElasticNet mostró R² test y CV prácticamente idénticos
(0.372 vs 0.372), indicando mejor generalización y ausencia de sobreajuste, aunque con
rendimiento global inferior al Random Forest.

#### 2.4 Comparativa final

| Modelo | R² test | R² CV (5-fold) |
|--------|---------|----------------|
| Random Forest (baseline) | 0.398 | 0.368 ± 0.069 |
| XGBoost | 0.394 | 0.307 ± 0.060 |
| Random Forest (tuned) | 0.400 | 0.344 |
| ElasticNet | 0.372 | 0.372 |

El **Random Forest tuned** fue seleccionado como modelo final por su mayor R² en test.
La convergencia de todos los modelos en torno a R²≈0.37-0.40 sugiere que este es el
techo predictivo alcanzable con dichos datos de expresión génica como única fuente de información,
apuntando a la necesidad de incorporar datos adicionales, como mutaciones somáticas.

### 3. Interpretabilidad (SHAP)

El análisis SHAP sobre el modelo Random Forest (R² = 0.40) permitió identificar los
genes con mayor impacto en la predicción de sensibilidad a Selumetinib, inhibidor de
MEK1/2 en la vía RAS/MAPK.

El gen con mayor importancia fue **TMSB15A**, cuya relación con la respuesta a
Selumetinib no está bien establecida en la literatura, por lo que requeriría análisis
adicionales para esclarecer su papel en este contexto.

Entre los genes con interpretación biológica más clara destacan **DUSP6** y **SPRY2**,
reguladores de feedback negativo de la vía RAS/MAPK. Su expresión se asoció con
cambios en la sensibilidad al fármaco, lo que resulta coherente con su mecanismo de
acción: una mayor actividad basal de la vía MAPK puede aumentar la dependencia celular
de esta señalización y, por tanto, su vulnerabilidad a la inhibición farmacológica.
En la misma línea, **ETV4**, un factor de transcripción regulado por ERK, refuerza la
evidencia de activación de esta vía.

Asimismo, la presencia de **MITF**, regulador clave en la diferenciación de melanocitos,
es consistente con el análisis exploratorio previo, donde el melanoma era el tipo tumoral
con mayor sensibilidad al tratamiento.

En conjunto, los resultados sugieren que el modelo ha capturado señal biológica relevante
respecto al mecanismo de acción de Selumetinib. No obstante, el valor de R² = 0.40
indica que solo se explica parcialmente la variabilidad en la respuesta al fármaco, lo
que apunta a la necesidad de incorporar información adicional. En este sentido, la
inclusión de datos de mutaciones somáticas podría aportar información complementaria
clave para mejorar la capacidad predictiva del modelo.

Estos resultados refuerzan el valor de la expresión génica como fuente de información
predictiva, aunque evidencian que la respuesta a fármacos oncológicos es un fenómeno
multifactorial.

## Estructura del repositorio

- `prediccion_respuesta_farmacos.ipynb`         # Pipeline completo
- `genes_eliminados.csv`   # Genes eliminados durante el preprocesado
- `README.md`

## Requisitos

- pandas
- numpy
- matplotlib
- seaborn
- scikit-learn
- xgboost
- shap
- openpyxl
