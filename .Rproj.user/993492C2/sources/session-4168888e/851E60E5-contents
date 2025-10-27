


install.packages("tidyverse")                                                   ## Comprovado
install.packages("here")
install.packages("fs")
install.packages("pheatmap")
install.packages("patchwork")
install.packages("knitr")
install.packages("rmarkdown")
install.packages("ggpubr")
install.packages("RColorBrewer")
install.packages("BiocManager")# Instalar el gestor de paquetes de Bioconductor
BiocManager::install("tximport")
BiocManager::install("DESeq2")                                                  ## Comprovado
#install.packages("ggthemes")
#BiocManager::install("ComplexHeatmap")


library(tidyverse)    # Manipulación (dplyr, tidyr), lectura (readr), gráficos (ggplot2)  ## Comprovado 
library(here)         # Gestión de directorios
library(fs) 
library(pheatmap)     # Heatmaps
library(patchwork)    # Combinar gráficos
library(knitr)        # Tablas dinámicas
library(rmarkdown)    # Informes reproducibles
library(ggpubr)       # Gráficos con estadísticas
library(RColorBrewer) # Paletas de colores
library(tximport)      # Para importar Salmon
library(DESeq2)        # Para                                ## Comprovado
#library(ggpubr)
#library(ComplexHeatmap)


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

# PASO 2: Verificar y guardar
head(colnames(txi$abundance))
head(rownames(metadata_df))

# PASO 3: Guardar matrices
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

# --------------------------------------
# Gráfico 3D: PCA interactivo con PC1, PC2 y PC3
# --------------------------------------
# Instalamos y cargamos 'plotly' si no está presente
if (!require("plotly", character.only = TRUE)) {
  install.packages("plotly", dependencies = TRUE)
  library(plotly, character.only = TRUE)
}

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


##################################################################

txi <- tximport(files, type = "salmon", txOut = TRUE)
dim(txi$abundance)
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
library(tidyverse)
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








