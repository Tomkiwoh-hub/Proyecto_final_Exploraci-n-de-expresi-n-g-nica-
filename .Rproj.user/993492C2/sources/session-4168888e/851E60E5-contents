


install.packages("tidyverse")                                                   ## Comprovado
install.packages("here")
install.packages("fs")
install.packages("pheatmap")
install.packages("patchwork")
install.packages("knitr")
install.packages("rmarkdown")
install.packages("ggpubr")
install.packages("RColorBrewer")
install.packages("plotly") 
install.packages("rgbif")
install.packages("sf")
install.packages("rworldxtra")
install.packages("geodata")
install.packages("ggspatial")
install.packages("terra")
install.packages("tidyterra")
install.packages("paletteer")
install.packages("ggcorrplot")
install.packages("ggridges")
install.packages("magick")
install.packages("BiocManager")# Instalar el gestor de paquetes de Bioconductor
BiocManager::install("tximport")
BiocManager::install("DESeq2")  


# Paquetes de CRAN
library(tidyverse)    # Manipulación (dplyr, tidyr), lectura (readr), gráficos (ggplot2) 
library(here)         # Gestión de directorios
library(fs)           # Herramientas para trabajar con el sistema de archivos (crear directorios, gestionar archivos)
library(pheatmap)     # Genera heatmaps para visualizar expresión génica
library(patchwork)    # Combina múltiples gráficos de ggplot2 en una sola figura
library(knitr)        # Tablas dinámicas
library(rmarkdown)    # Genera informes dinámicos y formateados, útil para RMarkdown
library(ggpubr)       # Simplifica la creación de gráficos de publicación con ggplot2
library(RColorBrewer) # Proporciona paletas de colores optimizadas para visualizaciones
library(rgbif) #descargar datos de ocurrencias
library(sf) #manipulaci?n  de datos vectoriales
library(rworldxtra) #datos vectoriales de los paises del mundo
library(geodata) #datos geoespaciales complemenatarios
library(ggspatial)#auxiliar para visualizar datos espaciales
library(terra) #datos raster
library(tidyterra) #maniipulaci?n de raster
library(paletteer) #colores
library(ggcorrplot) #diagrama de correlaciones
library(ggridges) #gr?fico de ridges
library(plotly) #gr?ficos avanzados
library(magick) #para manejo de imagenes
# Paquetes de Bioconductor
library(tximport)     # Para importar Salmon
library(DESeq2)       #Realiza análisis de expresión diferencial y normalización de datos de RNA-seq  ## Comprovado


getwd() #Comprovando direccion de carpeta "PROYECTO"                            ##comprovado
#E:/PROGRAMACION/PROYECTO1/PROYECTO                             
#Generando las carpetas
dir_create("datos")        # Para archivos TSV
dir_create("scripts")      # Para scripts .R
dir_create("salidas")      # Para resultados
dir_create("salidas_data") # Para datos procesados                              ##comprovado


# Creación de la Lista de Muestras (samples 18 SRR) 
samples <- c("SRR27137980", "SRR27137981", "SRR27137982", "SRR27137983", "SRR27137984",
             "SRR27137985", "SRR27137986", "SRR27137987", "SRR27137988", "SRR27137989",
             "SRR27137990", "SRR27137991", "SRR27137992", "SRR27137993", "SRR27137994",
             "SRR27137995", "SRR27137996", "SRR27137997")
print(samples)
# Rutas a quant.sf (en E:/PROGRAMACION/PROYECTO1/PROYECTO/datos/SRR271379**_quant/quant.sf)  ##comprobado
files <- file.path(here("datos"), paste0(samples, "_quant"), "quant.sf") #Construcción de las Rutas a los Archivos (files)
names(files) <- samples #Asignación de Nombres a las Rutas
file.exists(files) #Verificación de Existencia de Archivos                     ##comprobado = todos TRUE, fueron creados correctamente
names(files)
str(files)


metadata <- tibble(
  Sample = c("SRR27137980", "SRR27137981", "SRR27137982", "SRR27137983", "SRR27137984",
             "SRR27137985", "SRR27137986", "SRR27137987", "SRR27137988", "SRR27137989",
             "SRR27137990", "SRR27137991", "SRR27137992", "SRR27137993", "SRR27137994",
             "SRR27137995", "SRR27137996", "SRR27137997"),
  Condition = c("RNAi_OPI", "RNAi_LPI", "RNAi_LPI", "RNAi_LPI", "EV_OPI",
                "EV_OPI", "EV_OPI", "EV_LPI", "Ox_OPI", "Ox_OPI",
                "Ox_OPI", "Ox_LPI", "Ox_LPI", "Ox_LPI", "RNAi_OPI",
                "RNAi_OPI", "EV_LPI", "EV_LPI"),
  Pi_Status = c("OPI", "LPI", "LPI", "LPI", "OPI",
                "OPI", "OPI", "LPI", "OPI", "OPI",
                "OPI", "LPI", "LPI", "LPI", "OPI",
                "OPI", "LPI", "LPI"),
  Treatment = c("RNAi", "RNAi", "RNAi", "RNAi", "EV",
                "EV", "EV", "EV", "Ox", "Ox",
                "Ox", "Ox", "Ox", "Ox", "RNAi",
                "RNAi", "EV", "EV"),
  Replicate = c(1, 3, 2, 1, 3,
                2, 1, 3, 3, 2,
                1, 3, 2, 1, 3,
                2, 2, 1)
)

# Guardar en la raíz
write_csv(metadata, here("samples.csv"))
# O en datos
# write_csv(metadata, here("datos", "samples.csv"))


#Revisar el estado de la tabla samples.csv (metadatos)
metadata <-  read_csv(here ("samples.csv"))
summary(metadata)
view(metadata)
head(metadata)
getwd()  # Verificar working directory "E:/PROGRAMACION/PROYECTO1/PROYECTO"


# Cargar metadatos
metadata_df <- as.data.frame(metadata) # Convertir a data frame tradicional
row.names(metadata_df) <- metadata_df$Sample
metadata_df$Sample <- NULL  # Remover columna Sample (ya está en row names)

# Verificar estructura final
print("Estructura de metadata_df:")
str(metadata_df)
head(metadata)
print("Primeras filas:")
head(metadata_df)
head(rownames(metadata_df))


# Lista de muestras y rutas (verificar que existan)
samples <- rownames(metadata_df)  # Usar row names directamente
files <- file.path(here("datos"), paste0(samples, "_quant"), "quant.sf")
names(files) <- samples

#VERIFICACIÓN CRÍTICA: ¿Existen todos los archivos?
file_check <- file.exists(files)
print(file_check)
print(paste("Archivos encontrados:", sum(file_check), "de", length(file_check)))


# Importar con tximport (nivel transcrito)
txi <- tximport(files, type = "salmon", txOut = TRUE)
print(paste("Dimensiones TPM:", dim(txi$abundance)))
print(paste("Dimensiones Counts:", dim(txi$counts)))
colnames(txi$abundance)[1:18]
all(colnames(txi$abundance) == samples) # Verificar que las columnas coincidan con samples


dir.create(here("salidas_data"), showWarnings = FALSE)
tpm_matrix <- as.data.frame(txi$abundance) %>% rownames_to_column("Transcript")
counts_matrix <- as.data.frame(txi$counts) %>% rownames_to_column("Transcript")


write_csv(tpm_matrix, here("salidas_data", "tpm_matrix.csv"))
write_csv(counts_matrix, here("salidas_data", "counts_matrix.csv"))

counts <-  read_csv(here("salidas_data", "counts_matrix.csv"))
view (counts)
summary(counts)
head(counts)

tpm <-  read_csv(here("salidas_data", "tpm_matrix.csv")) 
view (tpm)
summary(tpm)
head(tpm)

print(paste("TPM:", nrow(tpm_matrix), "transcripts ×", ncol(tpm_matrix)-1, "samples"))
print(paste("Counts:", nrow(counts_matrix), "transcripts ×", ncol(counts_matrix)-1, "samples"))

####################################################
txi <- tximport(files, type = "salmon", txOut = TRUE)
dim(txi$abundance)  # Debe mostrar ~[transcripts, 18]

# Verificar y guardar
head(colnames(txi$abundance))
head(rownames(metadata_df))

# Guardar matrices
dir.create(here("salidas_data"), showWarnings = FALSE)
write_csv(as.data.frame(txi$abundance) %>% rownames_to_column("Transcript"), 
          here("salidas_data", "tpm_matrix.csv"))
write_csv(as.data.frame(txi$counts) %>% rownames_to_column("Transcript"), 
          here("salidas_data", "counts_matrix.csv"))

################################### Análisis Exploratorio (Después de tximport)

# Crear DESeqDataSet (AQUÍ SÍ usa metadata_df)
dds <- DESeqDataSetFromTximport(txi, colData = metadata_df, design = ~ Treatment + Pi_Status)

# Normalización para PCA/heatmap
vsd <- vst(dds)

# PCA
pca_data <- plotPCA(vsd, intgroup = c("Treatment", "Pi_Status"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = Treatment, shape = Pi_Status)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Phaseolus vulgaris - PHR1 RNAi/Ox/EV bajo Pi") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom")
ggsave(here("salidas", "pca_phr1.png"), pca_plot, width = 10, height = 8)

print(pca_plot)

# Guardar PCA
dir.create(here("salidas"), showWarnings = FALSE)
ggsave(here("salidas", "pca_phr1.png"), pca_plot, width = 10, height = 8, dpi = 300)
print(pca_plot)


# Gráfico 3D: PCA interactivo con PC1, PC2 y PC3

# Instalamos y cargamos 'plotly' si no está presente

# Calculamos PCA manualmente para incluir PC3
pca <- prcomp(t(assay(vsd)), scale. = FALSE)
pca_data_3d <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3 = pca$x[,3],
  Treatment = colData(vsd)$Treatment,
  Pi_Status = colData(vsd)$Pi_Status
)
percentVar_3d <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:3])

pca_plot_3d <- plot_ly(
  pca_data_3d,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Treatment,
  symbol = ~Pi_Status,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 8)
) %>%
  layout(
    title = "PCA 3D: Phaseolus vulgaris - PHR1 RNAi/Ox/EV bajo Pi",
    scene = list(
      xaxis = list(title = paste0("PC1: ", percentVar_3d[1], "% variance")),
      yaxis = list(title = paste0("PC2: ", percentVar_3d[2], "% variance")),
      zaxis = list(title = paste0("PC3: ", percentVar_3d[3], "% variance"))
    ),
    legend = list(orientation = "h", y = -0.1)
  )

# Guardamos el gráfico 3D como HTML
htmlwidgets::saveWidget(pca_plot_3d, here::here("salidas", "pca_phr1_3d.html"))
cat("Gráfico 3D PCA guardado en:", here::here("salidas", "pca_phr1_3d.html"), "\n")


####################################################

# Heatmap de top 50 genes más variables
top_var_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[top_var_genes, ]
mat <- t(scale(t(mat)))  # Z-score por gen

pheatmap(mat, 
         annotation_col = metadata_df[, c("Treatment", "Pi_Status"), drop = FALSE],
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Top 50 Transcritos Más Variables",
         fontsize_row = 8,
         filename = here("salidas", "heatmap_top50.png"),
         width = 12, height = 10)

##################################

# 
mat_long <- as.data.frame(mat) %>%
  rownames_to_column("Transcript") %>%
  pivot_longer(cols = -Transcript, names_to = "Sample", values_to = "Z_score") %>%
  left_join(metadata_df %>% rownames_to_column("Sample"), by = "Sample")

# 
box_plot <- ggplot(mat_long, aes(x = Treatment, y = Z_score, fill = Pi_Status)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Expression Distribution of Top 50 Variable Transcripts",
       x = "Treatment", y = "Z-score Expression") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom")

# 
ggsave(here("salidas", "boxplot_top50.png"), box_plot, width = 10, height = 6, dpi = 300)
print(box_plot)


# Sección 4: Generacion de un Mapa

# ===========================================================================
## 4.1 : Seleccion de especies
# ===========================================================================


consulta_A <- name_backbone("Nyctibius grandis") # esta funcion obtiene el resultado que mejor coincida

#vamos a trabajar con una especie del g?nero Pseudoeurycea
#puedes consultar informaci?n en https://enciclovida.mx/especies/25951-pseudoeurycea

consulta_B <- name_backbone("Nyctibius jamaicensis")

consulta_B


consulta_A2 <- name_suggest("Nyctibius grandis")$data

consulta_A2


consulta_B2 <- name_suggest("Nyctibius jamaicensis")$data

consulta_B2




# Descarga de ocurrencias -------------------------------------------------

Pse_lep <- occ_search(scientificName = "Nyctibius grandis",
                      hasCoordinate = TRUE,
                      hasGeospatialIssue = FALSE
)$data

Pse_lep


names(Pse_lep)


#vamos a explorar algunas variables
unique(Pse_lep$country)

#revisar de qu? tipo de registros se trata
unique(Pse_lep$basisOfRecord)


Pse_lep <- Pse_lep %>% 
  filter(basisOfRecord %in% c(c("HUMAN_OBSERVATION")))

#Instituci?n que hizo el registro
unique(Pse_lep$institutionCode)

Pse_lep %>% 
  ggplot(aes(x= institutionCode, fill= institutionCode))+
  geom_bar()+
  coord_flip()+
  theme(legend.position = "none")

#vamos a filtrar los datos que no tienen instituci?n de registro
Pse_lep <- Pse_lep %>% 
  filter(!is.na(institutionCode))


#CRS (Coordinate Reference System)
Pse_lep_sp <- Pse_lep %>% 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs= 4326)

#los CRS tienen c?digos ?nicos que los identifican (ESPG)
#el m?s com?n es el 4326 que corresponde a WGS 84 - WGS84 - Sistema Geod?sico Mundial 1984,
#https://epsg.io/4326

Pse_lep_sp


#cargo una capa raster de altitud
alt <- worldclim_global(var="elev", res=5, path=tempdir())

#cargo archivos vectoriales
data(countriesHigh)
#desino un objeto que tiene datos del mundo de manera espacial
Mundo <- st_as_sf(countriesHigh) 

##repetimos para la especie 
Inc_val <- occ_search(scientificName = "Nyctibius jamaicensis",
                      hasCoordinate = TRUE,
                      hasGeospatialIssue = FALSE
)$data


Inc_val$country %>% unique


Inc_val <- Inc_val %>% 
  filter(basisOfRecord %in% c(c("HUMAN_OBSERVATION"))) %>% 
  filter(!is.na(institutionCode))

Inc_val_sp <- Inc_val %>% 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs= 4326)


# ===========================================================================
## 4.2 : Mapas para especies
# ===========================================================================


#Mapa simple para la especie
ggplot()+
  geom_spatraster(data= alt)+
  geom_sf(data= Pse_lep_sp, aes(col = species), col="red4")+
  coord_sf(xlim = c(-107, -60), ylim = c(27, -10))+
  scale_fill_paletteer_c("grDevices::terrain.colors",
                         limits = c(0, 5000),
                         na.value = "transparent")

#repetir los pasos para la otra especie


ggplot()+
  geom_spatraster(data= alt)+
  geom_sf(data= Inc_val_sp, aes(col = species), col="blue4")+
  coord_sf(xlim = c(-115, -77), ylim = c(32, 10))+
  scale_fill_paletteer_c("grDevices::terrain.colors",
                         limits = c(0, 5000),
                         na.value = "transparent") +
  annotation_north_arrow(location = "bl",
                       which_north="true",
                       pad_x = unit(0.2, "in"),
                       pad_y = unit(0.7, "in"),
                       style=north_arrow_fancy_orienteering(fill = c("white", "grey60")))+
  annotation_scale(location = "bl",
                   bar_cols = c("grey60", "white"), 
                   text_family = "ArcherPro Book")+
  labs(fill= "Altitud (msnm)", col= "Especies")

#juntamos los dos mapas 
map_occ <-ggplot()+
  geom_spatraster(data= alt, alpha= 0.5)+
  geom_sf(data= Pse_lep_sp, aes(col = species))+
  geom_sf(data= Inc_val_sp, aes(col = species))+
  geom_sf(data= Mundo, fill= NA, linewidth=1)+
  coord_sf(xlim = c(-115, -60), ylim = c(32, -10))+
  scale_color_manual(values = c("#8d62fc", "#fc8d62")) +
  scale_fill_paletteer_c("grDevices::terrain.colors",
                         limits = c(0, 5000),
                         na.value = "transparent")+
  annotation_north_arrow(location = "bl",
                         which_north="true",
                         pad_x = unit(0.2, "in"),
                         pad_y = unit(0.7, "in"),
                         style=north_arrow_fancy_orienteering(fill = c("white", "grey60")))+
  annotation_scale(location = "bl",
                   bar_cols = c("grey60", "white"), 
                   text_family = "ArcherPro Book")+
  labs(fill= "Altitud (msnm)", col= "Especies")

map_occ







