library(Seurat)
obj=readRDS("../zebrafish_HighGlucosevsControl_merged.rds")
clusters=obj$seurat_clusters
clusters=data.frame(clusters)
write.table(clusters,file="cell_cluster_ids.tsv",sep="\t",row.names=T,col.names=F,quote=F)
