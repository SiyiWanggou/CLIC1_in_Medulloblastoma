#Combined RNAseq analysis for Michelle's CLIC1 project.
#This script includes RNAseq data from ONS76 and Vandy MB 11.
#------Charpt 1. Vocalno plot-------
#ONS76 cell line
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76/shCLIC1_vs_shScr_ons76/Statistics")
#shCLIC1_2_vs_shScr
shCLIC1_2_vs_shScr <- read.delim("DESeq2_normalized_counts_and_statistics_shClicl_2_vs_control.txt",header=TRUE,row.names = 1)
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
shCLIC1_3_vs_shScr <- read.delim("DESeq2_normalized_counts_and_statistics_shClicl_3_vs_control.txt",header=TRUE,row.names = 1)
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

#Vandy-MB-11 cell line
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76/shCLIC1_vs_shScr_vandy11/Statistics")
#shCLIC1_vs_shScr
shCLIC1_vs_shScr <- read.delim("Gene_DESeq2_normalized_counts_and_statistics.txt",header=TRUE,row.names = 1)
shCLIC1_vs_shScr <- na.omit(shCLIC1_vs_shScr)
library(ggplot2)
library(ggthemes)
library(Cairo)
shCLIC1_vs_shScr$threshold=as.factor(ifelse(shCLIC1_vs_shScr$padj<0.1 & abs(shCLIC1_vs_shScr$log2FoldChange) >= 1.0,
                                              ifelse(shCLIC1_vs_shScr$log2FoldChange > 1.0,
                                                     'Up Regulated in shCLIC1-1','Down Regulated in shCLIC1-1'),'None'))
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
#------Charpter 1 Done------

#------Charpter 2 Heatmap of DE genes------
#Identify CLIC1 target genes from ONS76 cell line.
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76/shCLIC1_vs_shScr_ons76/Statistics")
Ons76_shCLIC1_2_vs_shScr <- read.delim("Sig_DESeq2_merged_shClicl_2_vs_control.txt",header=TRUE,row.names = "geneNames")
Ons76_shCLIC1_3_vs_shScr <- read.delim("Sig_DESeq2_merged_shClicl_3_vs_control.txt",header=TRUE,row.names = "geneNames")
Ons76_shCLIC1_2_vs_shScr$geneNames <- row.names(Ons76_shCLIC1_2_vs_shScr)
Ons76_shCLIC1_3_vs_shScr$geneNames <- row.names(Ons76_shCLIC1_3_vs_shScr)
#Up or down regulated genes in shCLIC1_2 vs shScr
Up_Ons76_shCLIC1_2_vs_shScr <- subset(Ons76_shCLIC1_2_vs_shScr,log2FoldChange > 0)
Down_Ons76_shCLIC1_2_vs_shScr <- subset(Ons76_shCLIC1_2_vs_shScr,log2FoldChange < 0)
#Up or down regulated genes in shCLIC1_3 vs shScr
Up_Ons76_shCLIC1_3_vs_shScr <- subset(Ons76_shCLIC1_3_vs_shScr,log2FoldChange > 0)
Down_Ons76_shCLIC1_3_vs_shScr <- subset(Ons76_shCLIC1_3_vs_shScr,log2FoldChange < 0)
#Overlap the Up or down regulated genes in both shCLIC1_2 and shCLIC1_3
Up_Ons76_genes <- merge(Up_Ons76_shCLIC1_2_vs_shScr,Up_Ons76_shCLIC1_3_vs_shScr,by="geneNames",all=FALSE)
Down_Ons76_genes <- merge(Down_Ons76_shCLIC1_2_vs_shScr,Down_Ons76_shCLIC1_3_vs_shScr,by="geneNames",all=FALSE)
write.table(Up_Ons76_genes,"Up_in_Ons76.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
write.table(Down_Ons76_genes,"Down_in_Ons76.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)

#Identify CLIC1 target genes from VandyMB-11 cell line.
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76/shCLIC1_vs_shScr_vandy11/Statistics")
Vandy11_shCLIC1_vs_shScr <- read.delim("Results_significant_DE_Gene_level.txt",header=TRUE,row.names = 1)
#Up or down regulated genes in shCLIC1 vs shScr
Up_Vandy11_shCLIC1_vs_shScr <- subset(Vandy11_shCLIC1_vs_shScr,log2FoldChange > 0)
Down_Vandy11_shCLIC1_vs_shScr <- subset(Vandy11_shCLIC1_vs_shScr,log2FoldChange < 0)

#Generate CLIC1 target genes from both ONS76 and Vandy-MB-11 cell lines.
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76")
#option 1------Ons76_shCLIC1_2, Ons76_shCLIC1_3 and Vandy11_shCLIC1
Up_in_option1 <- merge(Up_Ons76_genes,Up_Vandy11_shCLIC1_vs_shScr,by="geneNames",all=FALSE)
Down_in_option1 <- merge(Down_Ons76_genes,Down_Vandy11_shCLIC1_vs_shScr,by="geneNames",all=FALSE)
#option 2------Ons76_shCLIC1_2 and Vandy11_shCLIC1
Up_in_option2 <- merge(Up_Ons76_shCLIC1_2_vs_shScr,Up_Vandy11_shCLIC1_vs_shScr,by="geneNames",all=FALSE)
Down_in_option2 <- merge(Down_Ons76_shCLIC1_2_vs_shScr,Down_Vandy11_shCLIC1_vs_shScr,by="geneNames",all=FALSE)
#option 3------Ons76_shCLIC1_3 and Vandy11_shCLIC1
Up_in_option3 <- merge(Up_Ons76_shCLIC1_3_vs_shScr,Up_Vandy11_shCLIC1_vs_shScr,by="geneNames",all=FALSE)
Down_in_option3 <- merge(Down_Ons76_shCLIC1_3_vs_shScr,Down_Vandy11_shCLIC1_vs_shScr,by="geneNames",all=FALSE)
#Save the Genelists
write.table(Up_in_option1,"Up_in_option1.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
write.table(Up_in_option2,"Up_in_option2.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
write.table(Up_in_option3,"Up_in_option3.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
write.table(Down_in_option1,"Down_in_option1.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
write.table(Down_in_option2,"Down_in_option2.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
write.table(Down_in_option3,"Down_in_option3.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
#Finally,we choose option1 as the CLIC1 target genes in both Vandy-MB-11 and Ons76.
#Plot Heatmap of Target genes from each design.
#Generate Matrix of DE Genes in Each Design.
#ONS76 shCLIC1-2 vs shScr
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76/shCLIC1_vs_shScr_ons76/Statistics")
Ons76_shCLIC1_2_vs_shScr <- read.delim("Sig_DESeq2_merged_shClicl_2_vs_control.txt",header=TRUE,row.names = "geneNames")
Matrix_Ons76_shCLIC1_2_vs_shScr <- Ons76_shCLIC1_2_vs_shScr[,2:7]
Ons76_shCLIC1_3_vs_shScr <- read.delim("Sig_DESeq2_merged_shClicl_3_vs_control.txt",header=TRUE,row.names = "geneNames")
Matrix_Ons76_shCLIC1_3_vs_shScr <- Ons76_shCLIC1_3_vs_shScr[,2:7]
#Vandy11-MB shCLIC1 vs shScr
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76/shCLIC1_vs_shScr_vandy11/Statistics")
Vandy11_shCLIC1_vs_shScr <- read.delim("Results_significant_DE_Gene_level.txt",header=TRUE,row.names = 1)
Matrix_Vandy11_shCLIC1_vs_shScr <- Vandy11_shCLIC1_vs_shScr[,1:7]
Matrix_Vandy11_shCLIC1_vs_shScr <- Matrix_Vandy11_shCLIC1_vs_shScr[!duplicated(Matrix_Vandy11_shCLIC1_vs_shScr[c("geneNames")]),]
row.names(Matrix_Vandy11_shCLIC1_vs_shScr) <- Matrix_Vandy11_shCLIC1_vs_shScr$geneNames
Matrix_Vandy11_shCLIC1_vs_shScr <- Matrix_Vandy11_shCLIC1_vs_shScr[,-1]

#Generate Target Gene Matrix for each design.
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76")
Up_target <- Up_in_option1$geneNames
Down_target <- Down_in_option1$geneNames
Target <- c(Up_target,Down_target)
#Matrix for Ons76_shCLIC1_2_vs_shScr
Matrix_1<- Matrix_Ons76_shCLIC1_2_vs_shScr[which(row.names(Matrix_Ons76_shCLIC1_2_vs_shScr) %in% Target),]
#Matrix for Ons76_shCLIC1_3_vs_shScr
Matrix_2<- Matrix_Ons76_shCLIC1_3_vs_shScr[which(row.names(Matrix_Ons76_shCLIC1_3_vs_shScr) %in% Target),]
#Matrix for Vandy-MB-11_shCLIC1_vs_shScr
Matrix_3<- Matrix_Vandy11_shCLIC1_vs_shScr[which(row.names(Matrix_Vandy11_shCLIC1_vs_shScr) %in% Target),]

#Merge for whole heatmap of Vandy_11
Matrix_Vandy_11 <- Matrix_3
Matrix_Vandy_11 <- t(scale(t(Matrix_Vandy_11)))
Cairo(file="Heatmap_plot_of_matrix_Vandy-11.png",type="png",units="in",bg="white",width=8,height=5,pointsize=16,dpi=300)
pheatmap(Matrix_Vandy_11)
dev.off()
#Merge for whole heatmap of Ons76
Matrix_1<- Matrix_Ons76_shCLIC1_2_vs_shScr[which(row.names(Matrix_Ons76_shCLIC1_2_vs_shScr) %in% Target),]
Matrix_2<- Matrix_Ons76_shCLIC1_3_vs_shScr[which(row.names(Matrix_Ons76_shCLIC1_3_vs_shScr) %in% Target),]
Matrix_ONS_76 <- cbind(Matrix_1,Matrix_2)
Matrix_ONS_76 <- Matrix_ONS_76[,-4:-6]
Matrix_ONS_76 <- t(scale(t(Matrix_ONS_76)))
Cairo(file="Heatmap_plot_of_Matrix_ONS_76.png",type="png",units="in",bg="white",width=8,height=5,pointsize=16,dpi=300)
pheatmap(Matrix_ONS_76)
dev.off()
#Save the matrix and plot heatmap in gitools.
write.table(Matrix_Vandy_11,"Matrix_Vandy_11.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
write.table(Matrix_ONS_76,"Matrix_ONS_76.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
#------Charpter 2 Done------

#------Charpter 3-------
#Enrichment analysis of DE genes across ONS76.
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76")
#Enrichment of Up-regulated genes in ONS76
up_genes <- read.delim("Up_in_Ons76.txt",header=TRUE,row.names = 1)
down_genes <- read.delim("Down_in_Ons76.txt",header=TRUE,row.names = 1)
library(clusterProfiler)
library(GSEABase)
gmtfileGO<-system.file("extdata","c5.bp.v7.0.symbols.gmt",package="clusterProfiler")
c5 <- read.gmt(gmtfileGO)
c5 <- enricher(up_genes$geneNames,TERM2GENE = c5,qvalueCutoff = 0.5)
write.table(c5,"enrichment_up_genes.txt",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
c5 <- read.gmt(gmtfileGO)
c5 <- enricher(down_genes$geneNames,TERM2GENE = c5)
write.table(c5,"enrichment_down_genes.txt",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")

#------Charpter 4.Total RNA readcounts-------
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76/shCLIC1_vs_shScr_ons76/Statistics/")
shCLIC1_2_vs_shScr <- read.delim("DESeq2_merged_shClicl_2_vs_control.txt",header = TRUE)
shCLIC1_3_vs_shScr <- read.delim("DESeq2_merged_shClicl_3_vs_control.txt",header = TRUE)
#Generate Gene Expression Readcounts shCLIC1_2 vs shScr
Matrix_shCLIC1_2_vs_shScr <- shCLIC1_2_vs_shScr[,2:8]
row.names(Matrix_shCLIC1_2_vs_shScr) <- Matrix_shCLIC1_2_vs_shScr$geneNames
write.table(Matrix_shCLIC1_2_vs_shScr,"Matrix_shCLIC1_2_vs_shScr.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
#Generate Gene Expression Readcounts shCLIC1_3 vs shScr
Matrix_shCLIC1_3_vs_shScr <- shCLIC1_3_vs_shScr[,2:8]
row.names(Matrix_shCLIC1_3_vs_shScr) <- Matrix_shCLIC1_3_vs_shScr$geneNames
write.table(Matrix_shCLIC1_3_vs_shScr,"Matrix_shCLIC1_3_vs_shScr.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
#load gene index of housekeeping genes
Gene_index <- read.delim("index.txt",header = TRUE)
#Generate housekeeping gene expression on both designs
Housekeeping_genes_in_design_1 <- merge(Matrix_shCLIC1_2_vs_shScr,Gene_index,by="geneNames",all = FALSE)
Housekeeping_genes_in_design_2 <- merge(Matrix_shCLIC1_3_vs_shScr,Gene_index,by="geneNames",all = FALSE)
#Save the housekeeping gene expression.
write.table(Housekeeping_genes_in_design_1,"Housekeeping_genes_in_design_1.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)
write.table(Housekeeping_genes_in_design_2,"Housekeeping_genes_in_design_2.txt",row.names = TRUE,col.names = TRUE,sep='\t',quote=FALSE)

#------Charpter 5. Analysis focusing on p38 cascade pathway------
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
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76/shCLIC1_vs_shScr_ons76/Statistics/")
shCLIC1_2_vs_shScr <- read.delim("DESeq2_merged_shClicl_2_vs_control.txt",header = TRUE)
shCLIC1_3_vs_shScr <- read.delim("DESeq2_merged_shClicl_3_vs_control.txt",header = TRUE)
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
write.table(kk,"gsea_output.txt", quote = F, row.names = F,sep='\t',col.names = TRUE)
#Plot GSEA results
Cairo(file="GO_p38_cascade_ONS76_shCLIC1_2_vs_shScr.png",type="png",units="in",bg="white",width=7.5,height=6,pointsize=16,dpi=300)
gseaplot2(kk, title="GO_P38MAPK_CASCADE", geneSetID = "GO_P38MAPK_CASCADE", color = "Red",ES_geom='line',base_size=15, rel_heights = c(1.5,0.5,1),subplots=1:3) 
dev.off()
Cairo(file="GO_POSITIVE_REGULATION_OF_P38MAPK_CASCADE_ONS76_shCLIC1_2_vs_shScr.png",type="png",units="in",bg="white",width=7.5,height=6,pointsize=16,dpi=300)
gseaplot2(kk, title="GO_POSITIVE_REGULATION_OF_P38MAPK_CASCADE", geneSetID = "GO_POSITIVE_REGULATION_OF_P38MAPK_CASCADE", color = "Red",ES_geom='line',base_size=15, rel_heights = c(1.5,0.5,1),subplots=1:3) 
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
write.table(kk,"gsea_output.txt", quote = F, row.names = F,sep='\t',col.names = TRUE)
#Plot GSEA results
Cairo(file="GO_p38_cascade_ONS76_shCLIC1_3_vs_shScr.png",type="png",units="in",bg="white",width=7.5,height=6,pointsize=16,dpi=300)
gseaplot2(kk, title="GO_P38MAPK_CASCADE", geneSetID = "GO_P38MAPK_CASCADE",
          color = "Red",ES_geom='line',base_size=15, rel_heights = c(1.5,0.5,1),subplots=1:3) 
dev.off()
Cairo(file="GO_POSITIVE_REGULATION_OF_P38MAPK_CASCADE_ONS76_shCLIC1_3_vs_shScr.png",type="png",units="in",bg="white",width=7.5,height=6,pointsize=16,dpi=300)
gseaplot2(kk, title="GO_POSITIVE_REGULATION_OF_P38MAPK_CASCADE", geneSetID = "GO_POSITIVE_REGULATION_OF_P38MAPK_CASCADE", 
          color = "Red",ES_geom='line',base_size=15, rel_heights = c(1.5,0.5,1),subplots=1:3,pvalue_table = FALSE) 
dev.off()

#vandyMB-11 cell lines
setwd("O:/Michelle_CLIC1/Combined_MB11_and_ONS76/shCLIC1_vs_shScr_vandy11/Statistics")
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
write.table(kk,"gsea_output.txt", quote = F, row.names = F,sep='\t',col.names = TRUE)
#Plot GSEA results
Cairo(file="GO_p38_cascade_Vandy11_shCLIC1_vs_shScr.png",type="png",units="in",bg="white",width=7.5,height=6,pointsize=16,dpi=300)
gseaplot2(kk, title="GO_P38MAPK_CASCADE", geneSetID = "GO_P38MAPK_CASCADE", color = "Red",ES_geom='line',base_size=15, rel_heights = c(1.5,0.5,1),subplots=1:3) 
dev.off()
Cairo(file="GO_POSITIVE_REGULATION_OF_P38MAPK_CASCADE_Vandy11_shCLIC1_vs_shScr.png",type="png",units="in",bg="white",width=7.5,height=6,pointsize=16,dpi=300)
gseaplot2(kk, title="GO_POSITIVE_REGULATION_OF_P38MAPK_CASCADE", geneSetID = "GO_POSITIVE_REGULATION_OF_P38MAPK_CASCADE", color = "Red",ES_geom='line',base_size=15, rel_heights = c(1.5,0.5,1),subplots=1:3) 
dev.off()


