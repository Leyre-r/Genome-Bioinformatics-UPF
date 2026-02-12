library(readxl)
library(ggplot2)
library(factoextra)
library(janitor)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(lubridate)
library(ggpubr)
library(patchwork)
library(pheatmap)
library(Rtsne)
library(RColorBrewer)

#df_metadata
df_metadata <- read_excel("data/raw_data/Table-S1.Additional-Demographical-and-Baseline-Characteristics-of-COVID-19-Patients-and-Control-Groups.xlsx", sheet = 2)
df_metadata <- clean_names(df_metadata)

df_metadata <- df_metadata %>%
  mutate(
    group = case_when(
      group_d == 0 ~ "Healthy Control",
      group_d == 1 ~ "Non-COVID-19",
      group_d == 2 ~ "Non-severe",
      group_d == 3 ~ "Severe",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group))  

#Cleaning variables used in df_metadata
df_metadata <- df_metadata %>%
  mutate(
    # sex coding
    sex = case_when(
      sex_g == 1 ~ "Male",
      sex_g == 0 ~ "Female",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Male","Female")),
    
    # numeric age
    age = parse_number(as.character(age_year)),
    
    # numeric BMI 
    bmi = na_if(as.character(bmi_h), "/") %>%
      str_replace_all(",", ".") %>%
      str_replace_all("[^0-9.\\-]", "") %>%
      na_if("") %>%
      parse_number(),
    
    # dates + time from onset to admission
    onset_date = as.Date(onset_date_f),
    admission_date = as.Date(admission_date),
    onset_to_admission = as.numeric(admission_date - onset_date + 1),
    
    # progression date
    prog_num = suppressWarnings(as.numeric(na_if(as.character(date_of_progression_to_severe_state), "/"))),
    prog_date = as.Date(prog_num, origin = "1899-12-30"),
    
    # time from admission to severe (only meaningful for Severe)
    admission_to_severe = as.numeric(prog_date - admission_date + 1)
  )

# df_proteom
df_proteom <- read_excel("data/raw_data/Table_S2_COVID.xlsx", sheet = 2, skip = 1) |> 
  rename("Proteins" = 1, "Gene Symbol" = 2) |>  
  clean_names()

#Removing NAs from df_proteom
data_proteins <- df_proteom[, 3:72]
es_valido <- (data_proteins != "NA")
proporcion_rellenos <- rowMeans(es_valido)
df_proteomics_filtrado <- df_proteom[proporcion_rellenos >= 0.85, ]

#df_metadata updated with patients kept after removing NAs from df_proteom
id_proteom <- str_to_upper(colnames(df_proteomics_filtrado[3:72]))
df_f_metadata <- df_metadata[df_metadata$ms_id_b %in% id_proteom,]
df_f2_metadata <- df_f_metadata[match(id_proteom,df_f_metadata$ms_id_b),]

#Adding patients' severity-level to df_proteomics
df_f2_proteomics <- df_proteomics_filtrado %>%
  mutate(across(3:72, as.numeric))
matrix_proteomics <- t(df_f2_proteomics[,3:72])
colnames(matrix_proteomics) <-  df_proteomics_filtrado$proteins
rownames(matrix_proteomics) <- str_to_upper(rownames(matrix_proteomics))
df_f2_proteomics <- as.data.frame(matrix_proteomics)
df_f2_proteomics <-  df_f2_proteomics %>% mutate(group = df_f2_metadata$group)

#Imputing of missing data 
df_f3_proteomics <- df_f2_proteomics %>%
  group_by(group) %>%
  mutate(across(where(is.numeric), ~ replace_na(., mean(., na.rm = TRUE)))) %>%
  ungroup()
df_f3_proteomics <- mutate(df_f3_proteomics, id_proteomics = rownames(matrix_proteomics))
df_f3_proteomics <- df_f3_proteomics[complete.cases(df_f3_proteomics),]


#Clustering based on patients' serum proteomic profiles

# K-means clustering
matrix_proteomics_scaled <- scale(df_f3_proteomics[,1:381], center = T, scale = T) 
rownames(matrix_proteomics_scaled) <-  df_f3_proteomics$id_proteomics

# Number of clusters evaluation
fviz_nbclust(matrix_proteomics_scaled, 
             kmeans, 
             method = "silhouette") +
  labs(subtitle = "Silhouette method",
       x = "Number of clusters (k)") +
  theme_minimal()

k4 <- kmeans(matrix_proteomics_scaled,centers = 4, iter.max = 10, nstart = 25)

# Hierarchical clustering
d <- dist(matrix_proteomics_scaled, method = "euclidean")
result_hclust <- hclust(d, method = "ward.D2") 

plot(result_hclust, cex = 0.6, main = "Hierarchical Clustering of Patients") 
rect.hclust(result_hclust, k = 4, border = "red")

groups <- cutree(result_hclust, k = 4)

# Comparison of clustering results
tabla_hcluster <- table(groups, df_f3_proteomics$group)
tabla_hcluster

tabla_kcluster <- table(groups = k4$cluster, df_f3_proteomics$group)
tabla_kcluster

compara_h_k <- table(Jerarquico = groups, KMeans = k4$cluster)
compara_h_k


#Dimensionality reduction

# Principal Component Analysis (PCA)
pca_results <- prcomp(matrix_proteomics_scaled, scale = TRUE)

#Percent of variance explained
percent <- 100*pca_results$sdev^2/sum(pca_results$sdev^2)
perc_data <- data.frame(Percent=percent, PC=1:length(percent))
ggplot(perc_data, aes(x=PC, y=Percent)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=round(percent, 2)), size=4, vjust=-.5) +
  ylim(0, 80)

#Data frame with PCA data
pca_data <- data.frame(pca_results$x, Patient_group=df_f3_proteomics$group)
#Plot coloring samples by patient group
ggplot(pca_data, aes(x=PC1, y=PC2, color = Patient_group)) + 
  geom_point() +
  theme_minimal()

#TOP 20 proteins for PC1
loadings <- as.data.frame(pca_results$rotation)
top20_pc1 <- head(loadings[order(abs(loadings$PC1), decreasing = TRUE), ], 20)
print(rownames(top20_pc1))

#TOP 20 proteins for PC2
top20_pc2 <- head(loadings[order(abs(loadings$PC2), decreasing = TRUE), ], 20)
print(rownames(top20_pc2))

# Common proteins in Top 20 PC1 and PC2
proteinas_comunes <- intersect(rownames(top20_pc1), rownames(top20_pc2))
print(proteinas_comunes) 


# t-SNE
set.seed(42)

tsne_result1 <- Rtsne(matrix_proteomics_scaled, dims = 2, perplexity=21, verbose=TRUE, max_iter = 500)
tsne_result2 <- Rtsne(matrix_proteomics_scaled, dims = 2, perplexity=10, verbose=TRUE, max_iter = 500)

# Creamos un data frame para ggplot
tsne_plot_data1 <- data.frame(X = tsne_result1$Y[,1], 
                              Y = tsne_result1$Y[,2], 
                              Group = df_f3_proteomics$group)

ggplot(tsne_plot_data1, aes(x = X, y = Y, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(
    data = tsne_plot_data1 %>% filter(Group %in% c("Healthy Control", "Severe", "Non-severe")),
    geom = "polygon", alpha = 0.1, level = 0.95) +
  theme_minimal() +
  labs(title = "t-SNE of Proteomics COVID-19",
       subtitle = "Perplexity 21")


#  UMAP
set.seed(42)
x <- df_f3_proteomics[,1:381] 
y <- df_f3_proteomics$group

umap_result <- uwot::umap( 
  x, n_neighbors = 15, 
  min_dist = 0.1, n_components = 2 )

umap_df <- data.frame(
  UMAP1 = umap_result[,1],
  UMAP2 = umap_result[,2],
  groups = y )

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = groups, fill=groups)) +
  geom_point(size = 2) +
  stat_ellipse(
    data = umap_df %>% filter(groups %in% c("Healthy Control", "Severe", "Non-severe")),
    geom = "polygon", alpha = 0.1, level = 0.95) +
  theme_minimal() +
  labs(title = "UMAP on COVID-19 Proteomic data")


# Heatmap visualization of the differentially expressed proteins between Severe and Healthy patients
df_metabolite <- read_excel("data/raw_data/Table_S6_Differentially-Expressed-Proteins-and-Metabolites.xlsx", sheet = 4)
df_metabolite <- clean_names(df_metabolite)

# Differentialy Expressed proteins from Table S6
proteins_DE <- df_metabolite %>% 
  filter(p_value_adjust < 0.05) %>% 
  filter(abs(fd) > 1) %>% 
  pull(proteins)

# Filtering proteins not present in the proteomic dataset 
index <- proteins_DE %in% colnames(df_f3_proteomics[1:381])
proteins_DE <- proteins_DE[index]

# Extracting proteomic data
heatmap_data <- df_f3_proteomics %>%
  select(c(proteins_DE, group, id_proteomics)) %>%
  filter(group %in% c("Severe", "Healthy Control"))
heatmap_data <-  heatmap_data %>% left_join(df_f2_metadata %>% select(ms_id_b, age), by = c("id_proteomics" = "ms_id_b")) 
heatmap_data <-  heatmap_data %>% left_join(df_f2_metadata %>% select(ms_id_b, sex), by = c("id_proteomics" = "ms_id_b")) 

heatmap_matrix <- as.matrix(heatmap_data[1:12])
rownames(heatmap_matrix) <- heatmap_data$id_proteomics

# Scale the data
heatmap_matrix <- scale(heatmap_matrix)

# Creating clinical annotations
annotation_df <- data.frame(
  Group = heatmap_data$group,
  Gender = heatmap_data$sex,
  Age = heatmap_data$age
)
rownames(annotation_df) <- heatmap_data$id_proteomics

# Sort patients by group
order_samples <-order(annotation_df$Group)

heatmap_matrix <- heatmap_matrix[order_samples, ]
annotation_df <- annotation_df[order_samples, ]

annotation_colors <- list(
  Group = c("Healthy Control" = "red","Severe" = "purple"),
  Gender = c("Male" = "orange","Female" = "lightblue"))
my_color_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

pheatmap(
  t(heatmap_matrix), 
  annotation_col = annotation_df,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,           
  show_rownames = TRUE,            
  fontsize_row = 7,              
  scale = "row",                   
  clustering_method = "ward.D2",
  cluster_cols = FALSE,
  color = my_color_palette,        
  border_color = NA,
  
  main = "Differentially Expressed Proteins: Severe Vs Healthy",
  treeheight_row = FALSE,            
  treeheight_col = FALSE,           
  cellwidth = 5,                  
  cellheight = 6                 
)

#Searching associations between protein expression and the annotations
#Common proteins: PC1-heatmap
top50_pc1 <- head(loadings[order(abs(loadings$PC1), decreasing = TRUE), ], 50)
index_pc1 <- rownames(top50_pc1) %in% colnames(heatmap_data[1:12])
print(rownames(top50_pc1)[index_pc1]) 
#Common proteins: PC2-heatmap
top50_pc2 <- head(loadings[order(abs(loadings$PC2), decreasing = TRUE), ], 50)
index_pc2 <- rownames(top50_pc2) %in% colnames(heatmap_data[1:12])
print(rownames(top50_pc2)[index_pc2])
