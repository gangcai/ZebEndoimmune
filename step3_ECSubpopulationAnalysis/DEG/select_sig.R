library(ggplot2)
data=read.table("HighGlucose_vs_Control_DEGs_WilCox.tsv",header=T,sep="\t")
filter=data$p_val_adj < 0.05 & abs(data$avg_log2FC)>1
data.f=data[filter,]

change.f=data.f$avg_log2FC > 0
change.type=sapply(change.f,function(x){
 if(x){return("UP")}else(return("DOWN"))
})

data.f.added=cbind(data.f,change.type)

write.table(data.f.added,file="HighGlucose_vs_Control_DEGs_WilCox_Sig.tsv",col.names=T,row.names=F,sep="\t",quote=F)

library(dplyr)
set.seed(1)
dat <- data.frame(ID = sample(letters,100,rep=TRUE))
data.f$cluster=factor(as.numeric(data.f$cluster))
sum.d <- data.f.added %>% 
    group_by(cluster, change.type) %>%
    summarise(no_rows = length(cluster))
write.table(sum.d,file="HighGlucose_vs_Control_DEGs_WilCox_Sig_Summary.tsv",
	    col.names=T,row.names=F,sep="\t",quote=F)
pdf("DEG_summary.pdf",height=5,width=6)
p=ggplot(data=sum.d, aes(x=as.character(cluster), y=no_rows, fill=change.type)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=no_rows), vjust=-0.3, color="black",
            position = position_dodge(0.9), size=3.5)+
   ylab("# of DEGs")+
   xlab("")+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()
print(p)
dev.off()
