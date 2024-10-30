library(ggplot2)
library(Seurat)
library(cowplot)
zeb.integrated <- readRDS("../zebrafish_HighGlucosevsControl_merged.rds")
clusters <- zeb.integrated$seurat_clusters
samples <- zeb.integrated$orig.ident
n1 <- sum(names(clusters) == names(samples))
n2 <- length(clusters)
#n1 == n2
cells <- names(clusters)
cell.df <- data.frame("cell"=cells, "sample"=samples, "cluster"=clusters)

result <- ""
i=0
for(sample in unique(samples)){
     n1 <- sum(samples == sample)
     for(cluster in unique(clusters)){
     	n2 <- sum(samples == sample & clusters == cluster)
        per <- 100*round(n2/n1,4)
	info <- c(sample,cluster,n1,n2,per)
	i <- i+1
	if (i == 1){
		result <- info
	}else{
		result <- rbind(result,info)
	}
     }
}

#times <- c(1,2,4,8)
#times.samples <- c("EC1dpf","EC2dpf","EC4dpf","EC8dpf")
colnames(result)=c("sample","cluster","total_cells","cluster_cells","percentage")
#experiment <- sub("EC","",result[,"sample"])
#experiment <- sub("dpf","",experiment)
experiment <- result[,"sample"]
result=cbind(result,experiment)
result=as.data.frame(result)

plot.list <- list()
i=0
clusters.u <- as.numeric(as.character(unique(clusters)))
clusters.u <- clusters.u[order(clusters.u)]
for(c in clusters.u){
	i=i+1
	#pdf(paste0(cluster,"_","cell_proportion_changes.pdf"))
	df <- subset(result,cluster==c)
        plot.list[[i]] <- ggplot(data=df, aes(x=experiment,y=as.numeric(as.character(percentage)), group = 1)) +
	    geom_line()+
	    ylab("percentage")+
	    ggtitle(c)+
	    theme( plot.title = element_text(color="darkorange", size=10,hjust = 0.5))+
	    geom_point()
	#print(p)
	#dev.off()
}


p <- plot_grid(plot.list[[1]],plot.list[[2]],plot.list[[3]],
	  plot.list[[4]],plot.list[[5]],ncol=3)

pdf("cluster_cell_proportion_changes.pdf",height=10)
print(p)
dev.off()

write.table(result,file="cluster_cell_proportions.tsv",sep="\t",row.names=F,quote=F)


p <- VlnPlot(zeb.integrated, features = c("marco","coro1a","mfap4","lyz","mpx"),pt.size= 0.02,ncol=2)
pdf("immune_marker_VlnPlot.pdf",height=10)
print(p)
dev.off()
p <- FeaturePlot(zeb.integrated, features = c("marco","coro1a","mfap4","lyz","mpx"),ncol=2)
pdf("immune_marker_FeaturePlot.pdf",height=10)
print(p)
dev.off()

#check kdrl gene
p <- FeaturePlot(zeb.integrated, features = c("kdrl"))
pdf("kdrl_FeaturePlot.pdf")
print(p)
dev.off()
