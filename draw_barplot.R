library(ggplot2)
library(Seurat)
library(cowplot)
zeb.integrated <- readRDS("zebrafish_endothelial_4stages_merged.rds")
ss <- summary(zeb.integrated@meta.data$samples)
samples <- zeb.integrated@meta.data$samples
samples <- as.factor(samples)
ss <- summary(samples)
df <- data.frame("embryo"=names(ss), "cells"=ss)
pdf("number_of_cells_per_sample.pdf")
p <- ggplot(data=df, aes(y=cells,x=embryo))+geom_bar(stat="identity",fill="steelblue")
print(p)
dev.off()
