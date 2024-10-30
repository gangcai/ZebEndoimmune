library(DESeq2)
library(Seurat)
library(ggplot2)
obj=readRDS("../zebrafish_HighGlucosevsControl_merged_addNorm.rds")
DefaultAssay(obj)="RNA"
#obj <- NormalizeData(obj,verbose = FALSE, normalization.method = "LogNormalize",
#				                 assay="RNA",scale.factor = 1e6) #log1p(RPM)

#obj <- ScaleData(obj,do.scale = TRUE , do.center = TRUE,
#			               assay="RNA", vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

obj$group=paste(Idents(obj),obj$samples,sep="_")
obj$celltype=Idents(obj)
Idents(obj)="group"
clusters=unique(obj$celltype)
samples=unique(obj$samples)
results=""
i=0
#log2FC: positive means higher in HighGlucose.
genes2=c()
for(cid in clusters){
   print(cid)
   g1=paste(cid,"HighGlucose",sep="_")
   g2=paste(cid,"Control",sep="_")
   high.glucose.response=FindMarkers(obj,ident.1=g1,ident.2=g2,
				     slot = "data",verbose=F,
				     assay= "RNA",  test.use = "wilcox",
				     pseudocount.use = 0.1)
   cluster=rep(cid,nrow(high.glucose.response))
   genes2=c(genes2,rownames(high.glucose.response))
   hgr=cbind(cluster,high.glucose.response)
   i=i+1
   if(i==1){
     results=hgr
   }else{
     results=rbind(results,hgr) #if the rownames appeared twice, a number will be added to the rownames, which makes the rownames not the same as the gene name
   }
}

#genes=rownames(results)
#results2=cbind(genes,results)
results3=cbind(genes2,results)
write.table(results3,file="HighGlucose_vs_Control_DEGs_WilCox.tsv",sep="\t",quote=F,row.names=F,col.names=T)


