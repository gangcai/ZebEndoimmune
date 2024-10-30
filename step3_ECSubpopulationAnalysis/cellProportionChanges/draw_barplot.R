library(ggplot2)
library(plyr)
library(RColorBrewer)
data=read.table("cluster_cell_proportions.tsv",header=T,sep="\t")
#df_sorted <- arrange(data, experiment, percentage) 
mycolors = c(brewer.pal(name="Dark2", n = 3), brewer.pal(name="Paired", n = 10))
#p=ggplot(data=df_cumsum, aes(x=experiment, y=percentage, fill=cluster)) +
#  geom_bar(stat="identity")+
#  geom_text(aes(y=label_ypos, label=percentage), vjust=1.6, 
#            color="white", size=3.5)+
#  scale_fill_brewer(palette="Paired")+
#  theme_minimal()
data$cluster=as.character(data$cluster)
data$cluster=factor(data$cluster,levels=c("0","1","2","3","4","5"))
p=ggplot(data=data, aes(x=experiment, y=percentage, fill=cluster)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=mycolors)+
  theme_minimal()

pdf("cluster_cellProportion_comparison.pdf",height=6,width=3)
print(p)

dev.off()
