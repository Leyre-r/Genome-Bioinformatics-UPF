library(ggplot2)
library(dplyr)
library(patchwork)

### CHIMPANzEE ###
valores_dN_ds_chimp <- read.delim(file="dN-dS_oct_2014_ortologs_chimpance.txt", header = TRUE)

#Adding dN/dS as new column nueva column and remove transcripts
valores_dN_ds_chimp <- data.frame(
  Ensembl.Gene.ID = c(valores_dN_ds_chimp$Ensembl.Gene.ID),
  dN = c(valores_dN_ds_chimp$dN),
  dS = c(valores_dN_ds_chimp$dS),
  dN_dS_chimp = c(valores_dN_ds_chimp$dN/valores_dN_ds_chimp$dS)
) 

#New data.frame with finite dN/dS genes 
reducido_dN_dS_chimp <- valores_dN_ds_chimp[is.finite(valores_dN_ds_chimp$dN_dS_chimp), ]
reducido_dN_dS_chimp <- unique(reducido_dN_dS_chimp)

#Olig2 gene dN/dS
dN_dS_olig2_chimp <- reducido_dN_dS_chimp[reducido_dN_dS_chimp$Ensembl.Gene.ID == "ENSG00000205927",4]

## Adding TF IDs
Human_TF <- read.delim("HumanTFs_DBD.txt")
FT_dN_dS_chimp <- inner_join(Human_TF, reducido_dN_dS_chimp, by = c("Ensembl.ID" = "Ensembl.Gene.ID"))

## Compaaring with other human bHLH genes
bHLH_genes <- Human_TF[Human_TF$DBD == "bHLH",]
bHLH_dN_dS_chimp <- inner_join(bHLH_genes, reducido_dN_dS_chimp, by = c("Ensembl.ID" = "Ensembl.Gene.ID"))

## Mean, median and percentile
chimp_dN_dS_statistics <- data.frame(Mean = c(Genomic = mean(reducido_dN_dS_chimp$dN_dS_chimp), 
                                              TF = mean(FT_dN_dS_chimp$dN_dS_chimp), 
                                              bHLH = mean(bHLH_dN_dS_chimp$dN_dS_chimp)),
                                     Median = c(Genomic = median(reducido_dN_dS_chimp$dN_dS_chimp),
                                                TF = median(FT_dN_dS_chimp$dN_dS_chimp),
                                                bHLH = median(bHLH_dN_dS_chimp$dN_dS_chimp)),
                                     Percentile = c(mean(reducido_dN_dS_chimp$dN_dS_chimp <= dN_dS_olig2_chimp),
                                                    mean(FT_dN_dS_chimp$dN_dS_chimp <= dN_dS_olig2_chimp),
                                                    mean(bHLH_dN_dS_chimp$dN_dS_chimp <= dN_dS_olig2_chimp)),
                                     Species_compared = c("chimpanzee", "chimpanzee", "chimpanzee"))

#Building Violin Plots
genomic_vplot_chimp <- ggplot(data = reducido_dN_dS_chimp, aes(y= dN_dS_chimp, x = "")) +
  geom_violin(fill ="orange", color = "black", alpha=0.5) +
  geom_boxplot(width = 0.25, color = "black", alpha = 0.3, outlier.shape = NA) +
  theme_classic() +
  labs(title = "dN/dS Human Vs Chimpanzee - Genome",
       y = "dN/dS",
       x = "Genomic Comparison (Human Vs Chimpanzee)") + 
  coord_cartesian(ylim=c(0,3)) +
  geom_hline(yintercept = dN_dS_olig2_chimp, color = "red")

TF_vplot_chimp <- ggplot(data = FT_dN_dS_chimp, aes(y= dN_dS_chimp, x = "")) +
  geom_violin(fill = "skyblue", color = "black", alpha=0.5) +
  geom_boxplot(width = 0.25, color = "black", alpha = 0.3, outlier.shape = NA) +
  theme_classic() +
  labs(title = "dN/dS Human Vs Chimpanzee - TF",
       y = "dN/dS",
       x = "TF Comparison (Human Vs Chimpanzee)") + 
  coord_cartesian(ylim=c(0,3)) +
  geom_hline(yintercept = dN_dS_olig2_chimp, color = "red")

bHLH_vplot_chimp <- ggplot(data = bHLH_dN_dS_chimp, aes(y= dN_dS_chimp, x = "")) +
  geom_violin(fill = "lightgreen", color = "black", alpha=0.5) +
  geom_boxplot(width = 0.25, color = "black", alpha = 0.3, outlier.shape = NA) +
  theme_classic() +
  labs(title = "dN/dS Human Vs Chimpanzee - bHLH",
       y = "dN/dS",
       x = "bHLH Comparison (Human Vs Chimpanzee)") + 
  coord_cartesian(ylim=c(0,3)) +
  geom_hline(yintercept = dN_dS_olig2_chimp, color = "red")

## Violin PLOTS
genomic_vplot_chimp
TF_vplot_chimp
bHLH_vplot_chimp
genomic_vplot_chimp | TF_vplot_chimp | bHLH_vplot_chimp



### MOUSE ###
valores_dN_ds_mouse <- read.delim(file="dN-dS_oct_2014_ortologs_mouse.txt", header = TRUE)

#Adding dN/dS as new column and removing quito transcripts
valores_dN_ds_mouse <- data.frame(
  Ensembl.Gene.ID = c(valores_dN_ds_mouse$Ensembl.Gene.ID),
  dN = c(valores_dN_ds_mouse$dN),
  dS = c(valores_dN_ds_mouse$dS),
  dN_dS_mouse = c(valores_dN_ds_mouse$dN/valores_dN_ds_mouse$dS)
) 

#New data.frame with finite dN/dS value genes
reducido_dN_dS_mouse <- valores_dN_ds_mouse[is.finite(valores_dN_ds_mouse$dN_dS_mouse), ]
reducido_dN_dS_mouse <- unique(reducido_dN_dS_mouse)

#Olig2 gene dN/dS
dN_dS_olig2_mouse <- reducido_dN_dS_mouse[reducido_dN_dS_mouse$Ensembl.Gene.ID == "ENSG00000205927",4]

##Adding TF IDs
FT_dN_dS_mouse <- inner_join(Human_TF, reducido_dN_dS_mouse, by = c("Ensembl.ID" = "Ensembl.Gene.ID"))

##Comparing with other bHLH 
bHLH_dN_dS_mouse <- inner_join(bHLH_genes, reducido_dN_dS_mouse, by = c("Ensembl.ID" = "Ensembl.Gene.ID"))

## Mean, Median and Percentile
mouse_dN_dS_statistics <- data.frame(Mean = c(Genomic = mean(reducido_dN_dS_mouse$dN_dS_mouse), 
                                        TF = mean(FT_dN_dS_mouse$dN_dS_mouse), 
                                        bHLH = mean(bHLH_dN_dS_mouse$dN_dS_mouse)),
                               Median = c(Genomic = median(reducido_dN_dS_mouse$dN_dS_mouse),
                                          TF = median(FT_dN_dS_mouse$dN_dS_mouse),
                                          bHLH = median(bHLH_dN_dS_mouse$dN_dS_mouse)),
                               Percentile = c(mean(reducido_dN_dS_mouse$dN_dS_mouse <= dN_dS_olig2_mouse),
                                              mean(FT_dN_dS_mouse$dN_dS_mouse <= dN_dS_olig2_mouse),
                                              mean(bHLH_dN_dS_mouse$dN_dS_mouse <= dN_dS_olig2_mouse)),
                               Species_compared = c("mouse", "mouse", "mouse"))


#Building Violin Plots
genomic_vplot_mouse <- ggplot(data = reducido_dN_dS_mouse, aes(y= dN_dS_mouse, x = "")) +
  geom_violin(fill = "orange", color = "black", alpha=0.5) +
  geom_boxplot(width = 0.25, color = "black", alpha = 0.3, outlier.shape = NA) +
  theme_classic() +
  labs(title = "dN/dS Human Vs Mouse - Genome",
       y = "dN/dS",
       x = "Genomic Comparison (Human Vs Mouse)") + 
  coord_cartesian(ylim=c(0,3)) +
  geom_hline(yintercept = dN_dS_olig2_mouse, color = "red")

TF_vplot_mouse <- ggplot(data = FT_dN_dS_mouse, aes(y= dN_dS_mouse, x = "")) +
  geom_violin(fill = "skyblue", color = "black", alpha=0.5) +
  geom_boxplot(width = 0.25, color = "black", alpha = 0.3, outlier.shape = NA) +
  theme_classic() +
  labs(title = "dN/dS Human Vs Mouse - TF",
       y = "dN/dS",
       x = "TF Comparison (Human Vs Mouse)") + 
  coord_cartesian(ylim=c(0,3)) +
  geom_hline(yintercept = dN_dS_olig2_mouse, color = "red")

bHLH_vplot_mouse <- ggplot(data = bHLH_dN_dS_mouse, aes(y= dN_dS_mouse, x = "")) +
  geom_violin(fill = "lightgreen", color = "black", alpha=0.5) +
  geom_boxplot(width = 0.25, color = "black", alpha = 0.3, outlier.shape = NA) +
  theme_classic() +
  labs(title = "dN/dS Human Vs Mouse - bHLH",
       y = "dN/dS",
       x = "bHLH Comparison (Human Vs Mouse)") + 
  coord_cartesian(ylim=c(0,3)) +
  geom_hline(yintercept = dN_dS_olig2_mouse, color = "red")

## Violin PLOTS
genomic_vplot_mouse
TF_vplot_mouse
bHLH_vplot_mouse
genomic_vplot_mouse | TF_vplot_mouse | bHLH_vplot_mouse
