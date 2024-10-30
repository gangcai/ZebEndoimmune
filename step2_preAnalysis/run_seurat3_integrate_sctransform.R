library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
options(future.globals.maxSize = 5000 * 1024^2) # 5G memory
raw.data.list=list()
samples=c("Control","HighGlucose")
raw.data.merged=""
matrix.dir="/home/db/private/CollaboratorLab/LiuDongLab/HighGlucose_vs_Control/RunCellRanger/"
i=0
for(sample in samples){
        dir=paste0(matrix.dir,sample,"/outs/filtered_feature_bc_matrix/")
	raw.data.list[[sample]]=Read10X(data.dir = dir)
	colnames=colnames(raw.data.list[[sample]])
	colnames(raw.data.list[[sample]])=paste0(sample,"_",colnames)
	i=i+1
	if(i == 1){
		raw.data.merged=raw.data.list[[sample]]
	}else{
		raw.data.merged=cbind(raw.data.merged,raw.data.list[[sample]])
	}
}

all.cells=colnames(raw.data.merged)
all.samples=sapply(all.cells,function(x){strsplit(as.character(x),"_")[[1]][1]})
all.samples=as.character(all.samples)
metadata = data.frame("samples"=all.samples,"cells"=all.cells,row.names=all.cells)
zeb <- CreateSeuratObject(raw.data.merged, meta.data = metadata, project = "zeb3k", min.cells = 3, min.features = 200)
zeb <- PercentageFeatureSet(zeb, pattern = "^mt-", col.name = "percent.mt")

p0 <- VlnPlot(zeb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, group.by = "samples")
pdf("zebrafish_HighGlucosevsControl_sctransform_QC_before.pdf")
print(p0)
dev.off()

#filter
zeb <- subset(zeb, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA < 50000)

zeb.list <- SplitObject(zeb, split.by = "samples")
for (i in 1:length(zeb.list)) {
    #zeb.list[[i]] <- NormalizeData(zeb.list[[i]], verbose = FALSE)
    #zeb.list[[i]] <- FindVariableFeatures(zeb.list[[i]], selection.method = "vst", 
    #    nfeatures = 2000, verbose = FALSE)
    zeb.list[[i]] <- SCTransform(zeb.list[[i]], verbose = FALSE,
				 vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
}

zeb.features <- SelectIntegrationFeatures(object.list = zeb.list, nfeatures = 3000)
zeb.list <- PrepSCTIntegration(object.list = zeb.list, anchor.features = zeb.features, 
				        verbose = FALSE)
zeb.anchors <- FindIntegrationAnchors(object.list = zeb.list, normalization.method = "SCT", 
					       anchor.features = zeb.features, verbose = FALSE)
zeb.integrated <- IntegrateData(anchorset = zeb.anchors, normalization.method = "SCT", 
				         verbose = FALSE)


#reference.list <- zeb.list[samples]
#zeb.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
#zeb.integrated <- IntegrateData(anchorset = zeb.anchors, dims = 1:30)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(zeb.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
#zeb.integrated <- ScaleData(zeb.integrated, verbose = FALSE)
zeb.integrated <- RunPCA(zeb.integrated, npcs = 30, verbose = FALSE)
zeb.integrated <- RunUMAP(zeb.integrated, reduction = "pca", dims = 1:30)
zeb.integrated <- RunTSNE(zeb.integrated, dims = 1:30, verbose = FALSE)
zeb.integrated <- FindNeighbors(zeb.integrated, dims = 1:30, verbose = FALSE)
zeb.integrated <- FindClusters(zeb.integrated,resolution = 0.15,  verbose = FALSE)

saveRDS(zeb.integrated, file = "zebrafish_HighGlucosevsControl_merged.rds")

p1 <- DimPlot(zeb.integrated, reduction = "umap", group.by = "samples")
p2 <- DimPlot(zeb.integrated, reduction = "tsne", group.by = "samples")
p3 <- DimPlot(zeb.integrated, label = TRUE, reduction="umap") + NoLegend()
p4 <- DimPlot(zeb.integrated, label = TRUE, reduction="tsne") + NoLegend()

#merged clusters with label
p5 <- DimPlot(zeb.integrated, split.by="samples",reduction="umap",ncol=2)
p6 <- DimPlot(zeb.integrated, split.by="samples",reduction="tsne",ncol=2)
#merged clusters without label
p7 <- DimPlot(zeb.integrated, split.by="samples",reduction="umap",label = TRUE ,ncol=2) + NoLegend()
p8 <- DimPlot(zeb.integrated, split.by="samples",reduction="tsne",label = TRUE, ncol=2) + NoLegend()


#QC plot
p9 <- VlnPlot(zeb.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, group.by = "samples")

#p2 <- DimPlot(zeb.integrated, reduction = "umap")
#p2 <- DimPlot(zeb.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
#    repel = TRUE) + NoLegend()
p <- plot_grid(p1, p2,p3,p4,ncol=2)
pdf("zebrafish_HighGlucosevsControl_merged_sctransform_umap_tsne.pdf")
print(p)
dev.off()

pdf("zebrafish_HighGlucosevsControl_split_sctransform_umap.pdf",width=12,height=6)
print(p5)
dev.off()

pdf("zebrafish_HighGlucosevsControl_split_sctransform_tsne.pdf",width=12,height=6)
print(p6)
dev.off()

pdf("zebrafish_HighGlucosevsControl_split_sctransform_umap_labeled.pdf",width=12,height=6)
print(p7)
dev.off()

pdf("zebrafish_HighGlucosevsControl_split_sctransform_tsne_labeled.pdf",width=12,height=6)
print(p8)
dev.off()

pdf("zebrafish_HighGlucosevsControl_sctransform_QC.pdf")
print(p9)
dev.off()


zeb.markers <- FindAllMarkers(zeb.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(zeb.markers,file="zebrafish_HighGlucosevsControl_merged_markers.tsv",sep="\t",quote=F,row.names=F)
