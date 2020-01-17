#----------
setwd("F:/science/science/2018_Xi_Huang_CLIC1/Datasets/GSE85217")
data<-read.delim("GSE85217.txt",header=TRUE)
library(psych)
library(pheatmap)
which(colnames(data) %in% c("CLIC1"))
ncol(data)
CLIC1_correlated_genes<-corr.test(data[3185],data[2:20197],adjust="fdr")
r<-CLIC1_correlated_genes$r
p<-CLIC1_correlated_genes$p
CLIC1MBsta<-rbind(r,p)
CLIC1MBsta<-t(CLIC1MBsta)
colnames(CLIC1MBsta)<-c('rvalue','pvalue')
CLIC1MBsta <- as.data.frame(CLIC1MBsta)
CLIC1MBsta$Gene <- rownames(CLIC1MBsta)
ionchannelgenelist<-read.delim("Index_for_ion_channel.txt",header=T)
ionchannelgenesta<-merge(CLIC1MBsta,ionchannelgenelist,by="Gene",all=FALSE)
ionchannelgenesig<-subset(ionchannelgenesta,pvalue < 0.05)
ionchannelgenesig<-subset(ionchannelgenesig,rvalue > 0.3 | rvalue<(-0.3))
rownames(data) <- data$Sample
data <- data[,-1]
datat<-t(data)
Sig_ion_matrix <- datat[which(rownames(datat) %in% ionchannelgenesig$Gene),]
zscore<-t(scale(t(Sig_ion_matrix)))
Cairo(file="heatmap_CLIC1_correlated_ion_channel_genes.png",type="png",units="in",bg="white",width=6,height=8,pointsize=16,dpi=300)
pheatmap(zscore,scale="column",cluster_cols=T,color=colorRampPalette(c("Green","black","red"))(50))
dev.off()
write.table(zscore,"zscore.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
###The zscore.txt file was used for heatmap visualization by gitools.###