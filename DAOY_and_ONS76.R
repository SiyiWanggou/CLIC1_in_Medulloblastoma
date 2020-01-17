#Combined RNAseq analysis for Michelle's CLIC1 project.
#This script includes RNAseq data from DAOY and ONS76.
#------Vocalno plot of DE genes-------
#ONS76 cell line
setwd("O:/Michelle_CLIC1/Combined_DAOY_and_ONS76/shCLIC1_vs_shScr_ONS76/Statistics")
#shCLIC1_2_vs_shScr
shCLIC1_2_vs_shScr <- read.delim("DESeq2_normalized_counts_and_statistics_shCLIC1_2_vs_control.txt",header=TRUE,row.names = 1)
shCLIC1_2_vs_shScr <- na.omit(shCLIC1_2_vs_shScr)
library(ggplot2)
library(ggthemes)
library(Cairo)
shCLIC1_2_vs_shScr$threshold=as.factor(ifelse(shCLIC1_2_vs_shScr$padj<0.1 & abs(shCLIC1_2_vs_shScr$log2FoldChange) >= 1.0,
                                              ifelse(shCLIC1_2_vs_shScr$log2FoldChange > 1.0,
                                                     'Up Regulated in shCLIC1-2','Down Regulated in shCLIC1-2'),'None'))
Cairo(file="Vocalno_plot_of_shCLIC1-2_vs_shScr.png",type="png",units="in",bg="white",width=8,height=5,pointsize=16,dpi=300)
ggplot(data=shCLIC1_2_vs_shScr, aes(x=log2FoldChange, y = -log10(padj), colour=threshold,fill=threshold)) +
  scale_color_manual(values=c("Green", "black","Red"))+
  geom_point(alpha=0.4, size=1.6) +
  xlim(c(-4, 4)) +
  theme_bw(base_size = 16, base_family = "Times") +
  geom_vline(xintercept=c(-1.0,1.0),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.1),lty=4,col="grey",lwd=0.6)+
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=16),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=16),
        axis.text.y = element_text(face="bold",  color="black", size=16),
        axis.title.x = element_text(face="bold", color="black", size=16),
        axis.title.y = element_text(face="bold",color="black", size=16))+
  labs(x="Log2 (fold change) ",y="-Log10 (Adjust p-value)",title="shCLIC1-2 vs shScr",size=16)
dev.off()
#shCLIC1_3_vs_shScr
shCLIC1_3_vs_shScr <- read.delim("DESeq2_normalized_counts_and_statistics_shCLIC1_3_vs_control.txt",header=TRUE,row.names = 1)
shCLIC1_3_vs_shScr <- na.omit(shCLIC1_3_vs_shScr)
library(ggplot2)
library(ggthemes)
library(Cairo)
shCLIC1_3_vs_shScr$threshold=as.factor(ifelse(shCLIC1_3_vs_shScr$padj<0.1 & abs(shCLIC1_3_vs_shScr$log2FoldChange) >= 1.0,
                                              ifelse(shCLIC1_3_vs_shScr$log2FoldChange > 1.0,
                                                     'Up Regulated in shCLIC1-3','Down Regulated in shCLIC1-3'),'None'))
Cairo(file="Vocalno_plot_of_shCLIC1-3_vs_shScr.png",type="png",units="in",bg="white",width=8,height=5,pointsize=16,dpi=300)
ggplot(data=shCLIC1_3_vs_shScr, aes(x=log2FoldChange, y = -log10(padj), colour=threshold,fill=threshold)) +
  scale_color_manual(values=c("Green", "black","Red"))+
  geom_point(alpha=0.4, size=1.6) +
  xlim(c(-4, 4)) +
  theme_bw(base_size = 16, base_family = "Times") +
  geom_vline(xintercept=c(-1.0,1.0),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.1),lty=4,col="grey",lwd=0.6)+
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=16),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=16),
        axis.text.y = element_text(face="bold",  color="black", size=16),
        axis.title.x = element_text(face="bold", color="black", size=16),
        axis.title.y = element_text(face="bold",color="black", size=16))+
  labs(x="Log2 (fold change) ",y="-Log10 (Adjust p-value)",title="shCLIC1-3 vs shScr",size=16)
dev.off()

#DAOY cell line
setwd("O:/Michelle_CLIC1/Combined_DAOY_and_ONS76/shCLIC1_vs_shScr_DAOY/Statistics")
#shCLIC1_vs_shScr
shCLIC1_vs_shScr <- read.delim("Gene_DESeq2_normalized_counts_and_statistics.txt",header=TRUE,row.names = 1)
shCLIC1_vs_shScr <- na.omit(shCLIC1_vs_shScr)
library(ggplot2)
library(ggthemes)
library(Cairo)
shCLIC1_vs_shScr$threshold=as.factor(ifelse(shCLIC1_vs_shScr$padj<0.1 & abs(shCLIC1_vs_shScr$log2FoldChange) >= 1.0,
                                            ifelse(shCLIC1_vs_shScr$log2FoldChange > 1.0,
                                                   'Up Regulated in shCLIC1','Down Regulated in shCLIC1'),'None'))
Cairo(file="Vocalno_plot_of_shCLIC1_vs_shScr.png",type="png",units="in",bg="white",width=8,height=5,pointsize=16,dpi=300)
ggplot(data=shCLIC1_vs_shScr, aes(x=log2FoldChange, y = -log10(padj), colour=threshold,fill=threshold)) +
  scale_color_manual(values=c("Green", "black","Red"))+
  geom_point(alpha=0.4, size=1.6) +
  xlim(c(-4, 4)) +
  theme_bw(base_size = 16, base_family = "Times") +
  geom_vline(xintercept=c(-1.0,1.0),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.1),lty=4,col="grey",lwd=0.6)+
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=16),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=16),
        axis.text.y = element_text(face="bold",  color="black", size=16),
        axis.title.x = element_text(face="bold", color="black", size=16),
        axis.title.y = element_text(face="bold",color="black", size=16))+
  labs(x="Log2 (fold change) ",y="-Log10 (Adjust p-value)",title="shCLIC1 vs shScr",size=16)
dev.off()
#------Done------

#------Analysis focusing on p38 cascade pathway------
#The Gene expression profile from both cell lines are sent to GSEA software for enrichment analysis.
#We used enrichplot here to plot individualized pathways.
library(clusterProfiler)
library(enrichplot)
library(plyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(GSEABase)
library(Cairo)
#ONS76 cell lines
setwd("O:/Michelle_CLIC1/Combined_DAOY_and_ONS76/shCLIC1_vs_shScr_ONS76/Statistics")
shCLIC1_2_vs_shScr <- read.delim("DESeq2_merged_shCLIC1_2_vs_control.txt",header = TRUE)
shCLIC1_3_vs_shScr <- read.delim("DESeq2_merged_shCLIC1_3_vs_control.txt",header = TRUE)
#For shCLIC1_2_vs_shScr
Rank<-data.frame(shCLIC1_2_vs_shScr$geneNames,shCLIC1_2_vs_shScr$log2FoldChange)
Rank <- na.omit(Rank)
colnames(Rank)<-c("SYMBOL","SCORE")
write.table(Rank,"Rank.txt",col.names = T,sep='\t',quote=FALSE,row.names=FALSE)
Rank<-read.delim("Rank.txt",header=T)
Rank_id<-bitr(Rank$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
Rank_id_score<-merge(Rank,Rank_id,by="SYMBOL",all=F)
Rank_id_score_sorted<-Rank_id_score[order(Rank_id_score$SCORE,decreasing = T),]
id_score<-Rank_id_score_sorted$SCORE
names(id_score)<-Rank_id_score_sorted$SYMBOL
gmtfileGO<-system.file("extdata","c5.bp.v7.0.symbols.gmt",package="clusterProfiler")
c5 <- read.gmt(gmtfileGO)
kk<-GSEA(id_score, TERM2GENE = c5, pvalueCutoff = 1,verbose=TRUE)
#Plot p38 MAPK pathway
Cairo(file="GO_p38_cascade_ONS76_shCLIC1_2_vs_shScr.png",type="png",units="in",bg="white",width=7.5,height=6,pointsize=16,dpi=300)
gseaplot2(kk, title="GO_P38MAPK_CASCADE", geneSetID = "GO_P38MAPK_CASCADE", color = "Red",ES_geom='line',base_size=15, rel_heights = c(1.5,0.5,1),subplots=1:3) 
dev.off()

#For shCLIC1_3_vs_shScr
Rank<-data.frame(shCLIC1_3_vs_shScr$geneNames,shCLIC1_3_vs_shScr$log2FoldChange)
Rank <- na.omit(Rank)
colnames(Rank)<-c("SYMBOL","SCORE")
write.table(Rank,"Rank.txt",col.names = T,sep='\t',quote=FALSE,row.names=FALSE)
Rank<-read.delim("Rank.txt",header=T)
Rank_id<-bitr(Rank$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
Rank_id_score<-merge(Rank,Rank_id,by="SYMBOL",all=F)
Rank_id_score_sorted<-Rank_id_score[order(Rank_id_score$SCORE,decreasing = T),]
id_score<-Rank_id_score_sorted$SCORE
names(id_score)<-Rank_id_score_sorted$SYMBOL
gmtfileGO<-system.file("extdata","c5.bp.v7.0.symbols.gmt",package="clusterProfiler")
c5 <- read.gmt(gmtfileGO)
kk<-GSEA(id_score, TERM2GENE = c5, pvalueCutoff = 1,verbose=TRUE)
#Plot p38 MAPK pathway
Cairo(file="GO_p38_cascade_ONS76_shCLIC1_3_vs_shScr.png",type="png",units="in",bg="white",width=7.5,height=6,pointsize=16,dpi=300)
gseaplot2(kk, title="GO_P38MAPK_CASCADE", geneSetID = "GO_P38MAPK_CASCADE",
          color = "Red",ES_geom='line',base_size=15, rel_heights = c(1.5,0.5,1),subplots=1:3) 
dev.off()

#DAYO cell lines
setwd("O:/Michelle_CLIC1/Combined_DAOY_and_ONS76/shCLIC1_vs_shScr_DAOY/Statistics")
shCLIC1_vs_shScr <- read.delim("Gene_DESeq2_normalized_counts_and_statistics.txt",header = TRUE)
Rank<-data.frame(shCLIC1_vs_shScr$geneNames,shCLIC1_vs_shScr$log2FoldChange)
Rank <- na.omit(Rank)
colnames(Rank)<-c("SYMBOL","SCORE")
write.table(Rank,"Rank.txt",col.names = T,sep='\t',quote=FALSE,row.names=FALSE)
Rank<-read.delim("Rank.txt",header=T)
Rank_id<-bitr(Rank$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
Rank_id_score<-merge(Rank,Rank_id,by="SYMBOL",all=F)
Rank_id_score_sorted<-Rank_id_score[order(Rank_id_score$SCORE,decreasing = T),]
id_score<-Rank_id_score_sorted$SCORE
names(id_score)<-Rank_id_score_sorted$SYMBOL
gmtfileGO<-system.file("extdata","c5.bp.v7.0.symbols.gmt",package="clusterProfiler")
c5 <- read.gmt(gmtfileGO)
kk<-GSEA(id_score, TERM2GENE = c5, pvalueCutoff = 1,verbose=TRUE)
#Plot GSEA results
Cairo(file="GO_p38_cascade_DAOY_shCLIC1_vs_shScr.png",type="png",units="in",bg="white",width=7.5,height=6,pointsize=16,dpi=300)
gseaplot2(kk, title="GO_P38MAPK_CASCADE", geneSetID = "GO_P38MAPK_CASCADE", color = "Red",ES_geom='line',base_size=15, rel_heights = c(1.5,0.5,1),subplots=1:3) 
dev.off()

