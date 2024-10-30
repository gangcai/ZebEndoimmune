library(Seurat)
library(cowplot)
#library("RColorBrewer")
library(dplyr)
library(ggplot2)
#mycolor=colorRampPalette(c("tomato2","snow1","skyblue2"))(10)
mycolor=colorRampPalette(c("skyblue2","tomato2"))(10)
obj = readRDS("../zebrafish_HighGlucosevsControl_merged.rds")
obj <- NormalizeData(obj,verbose = FALSE, normalization.method = "LogNormalize",
                                          assay="RNA",scale.factor = 1e6) #log1p(RPM)

obj <- ScaleData(obj,do.scale = TRUE , do.center = TRUE,
                                  assay="RNA")
clusters=levels(obj)
DefaultAssay(obj) <- "RNA"

marker.data=read.table("../../EC_Markers.txt",header=T,sep="\t")
marker.n=nrow(marker.data)
for(i in c(1:marker.n)){
	gene=marker.data[i,"GeneName"]
	cell=marker.data[i,"CellType"]
        try({	
		p = FeaturePlot(obj, features = gene, cols = mycolor,slot = "scale.data",ncol=1)
		pdf(paste0(cell,"_",gene,"_","zebTestis_FeaturePlot.pdf"),height=5,width=5)
		print(p)
		dev.off()

		p = VlnPlot(obj,assay = "RNA", features = gene,slot = "scale.data", ncol=1,pt.size=0.5)
		pdf(paste0(cell,"_",gene,"_","zebTestis_VlnPlot.pdf"),height=5,width=5)
		print(p)
		dev.off()
	})
}
