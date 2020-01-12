setwd("Z:/Michelle_Clic1/Y4BT8P1/CHE11495/Ballgawn_and_DEseq2")
library(ballgown)
pheno_data<-read.delim("phenotype.txt")
bg<-ballgown(dataDir="shCLIC1_vs_shCtrl",samplePattern="",pData=pheno_data)
names<-data.frame(geneNames=ballgown::geneNames(bg),geneIDs=ballgown::geneIDs(bg))
setwd("Z:/Michelle_Clic1/Y4BT8P1/CHE11495/Stringtie_DESeq2")
#Compare DE genes between shCLIC1_2 vs Control
coldata<-read.delim("phenotype_shClic1_2_vs_control.txt",header=T,row.names = 1)
cts_all<-read.csv("gene_count_matrix.csv",header=T,row.names="gene_id")
cts_shClic1_2_vs_control <- cts_all[,c(1:3,7:9)]
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts_shClic1_2_vs_control))
all(rownames(coldata) == colnames(cts_shClic1_2_vs_control))
library(DESeq2)
dds_shClic1_2_vs_control<-DESeqDataSetFromMatrix(countData=cts_shClic1_2_vs_control,colData=coldata,design= ~ phenotype)
dds_shClic1_2_vs_control<-DESeq(dds_shClic1_2_vs_control)
res_shClic1_2_vs_control<-results(dds_shClic1_2_vs_control)
normalized_counts<-counts(dds_shClic1_2_vs_control,normalized=TRUE)
merged_file<-data.frame(normalized_counts,res_shClic1_2_vs_control)
write.table(merged_file,"DESeq2_normalized_counts_and_statistics_shClicl_2_vs_control.txt",row.names = T,col.names = T,sep='\t',quote=FALSE)
results<-read.delim("DESeq2_normalized_counts_and_statistics_shClicl_2_vs_control.txt",header=T)
results$geneIDs<-rownames(results)
c<-merge(names,results,by="geneIDs",all=FALSE)
c<-c[!duplicated(c[c("geneNames")]),]
write.table(c,"DESeq2_merged_shClicl_2_vs_control.txt",row.names = FALSE,col.names = T,sep='\t',quote=FALSE)

#Compare DE genes between shCLIC1_3 vs Control
coldata<-read.delim("phenotype_shClic1_3_vs_control.txt",header=T,row.names = 1)
cts_all<-read.csv("gene_count_matrix.csv",header=T,row.names="gene_id")
cts_shClic1_3_vs_control <- cts_all[,c(4:6,7:9)]
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts_shClic1_3_vs_control))
all(rownames(coldata) == colnames(cts_shClic1_3_vs_control))
library(DESeq2)
dds_shClic1_3_vs_control<-DESeqDataSetFromMatrix(countData=cts_shClic1_3_vs_control,colData=coldata,design= ~ phenotype)
dds_shClic1_3_vs_control<-DESeq(dds_shClic1_3_vs_control)
res_shClic1_3_vs_control<-results(dds_shClic1_3_vs_control)
normalized_counts<-counts(dds_shClic1_3_vs_control,normalized=TRUE)
merged_file<-data.frame(normalized_counts,res_shClic1_3_vs_control)
write.table(merged_file,"DESeq2_normalized_counts_and_statistics_shClicl_3_vs_control.txt",row.names = T,col.names = T,sep='\t',quote=FALSE)
results<-read.delim("DESeq2_normalized_counts_and_statistics_shClicl_3_vs_control.txt",header=T)
results$geneIDs<-rownames(results)
c<-merge(names,results,by="geneIDs",all=FALSE)
c<-c[!duplicated(c[c("geneNames")]),]
write.table(c,"DESeq2_merged_shClicl_3_vs_control.txt",row.names = FALSE,col.names = T,sep='\t',quote=FALSE)

#vennplot to find high confident genes
Des_1 <- read.delim("DESeq2_merged_shClicl_2_vs_control.txt",header = T)
Des_2 <- read.delim("DESeq2_merged_shClicl_3_vs_control.txt",header = T)
Sig_Des_1 <- subset(Des_1,padj < 0.1)
Sig_Des_1 <- subset(Sig_Des_1, log2FoldChange > 1 | log2FoldChange < (-1))
Sig_Des_2 <- subset(Des_2,padj < 0.1)
Sig_Des_2 <- subset(Sig_Des_2, log2FoldChange > 1 | log2FoldChange < (-1))
write.table(Sig_Des_1,"Sig_DESeq2_merged_shClicl_2_vs_control.txt",row.names = FALSE,col.names = T,sep='\t',quote=FALSE)
write.table(Sig_Des_2,"Sig_DESeq2_merged_shClicl_3_vs_control.txt",row.names = FALSE,col.names = T,sep='\t',quote=FALSE)

#Identify high confident CLIC1 target genes
Up_Sig_Des_1 <- subset(Sig_Des_1, log2FoldChange > 1)
Down_Sig_Des_1 <- subset(Sig_Des_1, log2FoldChange < (-1))
Up_Sig_Des_2 <- subset(Sig_Des_2, log2FoldChange > 1)
Down_Sig_Des_2 <- subset(Sig_Des_2, log2FoldChange < (-1))
Up_Genes_final <- merge(Up_Sig_Des_1,Up_Sig_Des_2, by="geneNames",all=FALSE)
Down_Genes_final <- merge(Down_Sig_Des_1,Down_Sig_Des_2, by="geneNames",all=FALSE)
High_Confident_Genes <- rbind(Up_Genes_final,Down_Genes_final)
write.table(High_Confident_Genes,"High_confident_CLIC1_target_genes_in_ons76.txt",row.names = FALSE,col.names = T,sep='\t',quote=FALSE)


