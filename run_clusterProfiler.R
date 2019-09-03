library("clusterProfiler")
args = commandArgs(trailingOnly=TRUE)
cluster.choose=args[1]
cluster.choose=as.character(cluster.choose)

print(cluster.choose)

data=read.table("../zebrafish_endothelial_4stages_merged_markers.tsv",sep="\t",header=T)
#filter=data$cluster == cluster.choose & data$p_val_adj < 0.1
#genes=as.character(data[filter,"gene"])
genes=subset(data,cluster == cluster.choose & p_val_adj < 0.1, select = gene)
genes=as.character(genes$gene)

#check keyType
#library(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
#mouse: org.Mm.eg.db
eg = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Dr.eg.db")
gene_id=as.character(eg$ENTREZID)


#GO
#ggo <- groupGO(gene     = gene_id,
#               OrgDb    = org.Dr.eg.db,
#               ont      = "CC",
#               level    = 3,
#               readable = TRUE)

#pdf(paste0("GOEnrichment_Barplot_cluster_",cluster,".pdf"))
#barplot(ggo, drop=TRUE, showCategory=15)
#dev.off()

for(ont_type in c("CC","MF","BP")){
	ego <- enrichGO(gene          = gene_id,
			OrgDb         = org.Dr.eg.db,
			ont           = ont_type,
			pAdjustMethod = "BH",
			pvalueCutoff  = 0.05,
		        readable      = TRUE)

	pdf(paste0("GOEnrichment_",ont_type,"_DotPlot_cluster_",cluster.choose,".pdf"),width=14)
	d=dotplot(ego)
	print(d)
	dev.off()
}

#KEGG
kk <- enrichKEGG(gene         = gene_id,
                 organism     = 'dre',
                 pvalueCutoff = 0.05)
pdf(paste0("KEGGEnrichment_DotPlot_cluster_",cluster.choose,".pdf"))
dotplot(kk,showCategory=20)
dev.off()

