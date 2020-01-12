#-----shCLIC1_vs_shScr------
setwd("O:/CLIC1_RNAseq/Ballgawn")
library(ballgown)
pheno_data<-read.csv("shCLIC1_vs_shScr_phenotype.csv")
bg<-ballgown(dataDir="shCLIC1_vs_shScr",samplePattern="CHE8695_",pData=pheno_data)
names<-data.frame(geneNames=ballgown::geneNames(bg),geneIDs=ballgown::geneIDs(bg))
names_unique<-names[!duplicated(names[c("geneIDs")]),]
setwd("O:/CLIC1_RNAseq/Ballgown_to_DEseq2/shCLIC1_vs_shScr")
write.table(names_unique,"geneNames.txt",row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
coldata<-read.delim("shCLIC1_vs_shScr_phenotype.txt",header=T,row.names = "Sample")
cts<-read.delim("gene_count_matrix.txt",header=T,row.names="gene_id")
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

#Analysis at gene level
library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData=cts,colData=coldata,design= ~ Condition)
dds<-dds[rowSums(counts(dds))>=10,]
dds<-DESeq(dds)
res<-results(dds)
normalized_counts<-counts(dds,normalized=TRUE)
merged_file<-data.frame(normalized_counts,res)
write.table(res,"DESeq2_results_statistics.txt",row.names = T,col.names = T,sep='\t',quote=FALSE)
write.table(normalized_counts,"DESeq2_normalized_counts.txt",row.names = T,col.names = T,sep='\t',quote=FALSE)
write.table(merged_file,"DESeq2_normalized_counts_and_statistics.txt",row.names = T,col.names = T,sep='\t',quote=FALSE)
results<-read.delim("DESeq2_normalized_counts_and_statistics.txt",header=T,row.names = 1)
results$geneIDs<-row.names(results)
names<-read.delim("geneNames.txt",header=T)
c<-merge(names,results,by="geneIDs",all=FALSE)
write.table(c,"Gene_DESeq2_normalized_counts_and_statistics.txt",row.names = FALSE,col.names = T,sep='\t',quote=FALSE)
c_sig<-subset(c,padj<0.05)
c_sig<-subset(c_sig,log2FoldChange > 1 | log2FoldChange < -1)
write.table(c_sig,"Results_significant_DE_Gene_level.txt",row.names = FALSE,col.names = T,sep='\t',quote=FALSE)

#Analysis at transcripts level
setwd("O:/CLIC1_RNAseq/Ballgawn")
library(ballgown)
pheno_data<-read.csv("shCLIC1_vs_shScr_phenotype.csv")
bg<-ballgown(dataDir="shCLIC1_vs_shScr",samplePattern="CHE8695_",pData=pheno_data)
bg_table_transcripts=texpr(bg,'all')

setwd("O:/CLIC1_RNAseq/Ballgown_to_DEseq2/shCLIC1_vs_shScr")
coldata<-read.delim("shCLIC1_vs_shScr_phenotype.txt",header=T,row.names = "Sample")
cts<-read.delim("transcript_count_matrix.txt",header=T,row.names="transcript_id")
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData=cts,colData=coldata,design= ~ Condition)
dds<-dds[rowSums(counts(dds))>=10,]
dds<-DESeq(dds)
res<-results(dds)
normalized_counts<-counts(dds,normalized=TRUE)
merged_file<-data.frame(normalized_counts,res)
merged_file$t_name<-row.names(merged_file)
combined_results<-merge(merged_file,bg_table_transcripts,by="t_name",all=FALSE)
combined_results_sig<-subset(combined_results,padj<0.05)
combined_results_sig<-subset(combined_results_sig,log2FoldChange > 1 | log2FoldChange < -1)

write.table(combined_results_sig,"Results_significant_DE_transcripts_level.txt",col.names = TRUE,row.names = TRUE,sep="\t",quote=FALSE)
write.table(combined_results,"Results_transcripts_all.txt",col.names = TRUE,row.names = TRUE,sep="\t",quote=FALSE)

#Visualization of Vocalno plot
setwd("O:/CLIC1_RNAseq/Ballgown_to_DEseq2/Visualization_of_shCLIC1_vs_shScr/")
a<-read.delim("Gene_DESeq2_normalized_counts_and_statistics.txt",header=T)
a<-na.omit(a)
library(ggplot2)
library(ggthemes)
library(Cairo)
a$threshold=as.factor(ifelse(a$padj<0.1 & abs(a$log2FoldChange) >= 1.0,
                             ifelse(a$log2FoldChange > 1.0,
                                    'Up Regulated in shCLIC1','Down Regulated in shCLIC1'),'None'))
Cairo(file="Vocalno_plot_of_shCLIC1_vs_shScr.png",type="png",units="in",bg="white",width=8,height=6,pointsize=16,dpi=300)
ggplot(data=a, aes(x=log2FoldChange, y = -log10(padj), colour=threshold,fill=threshold)) +
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
  labs(x="Log2 (fold change) ",y="-Log10 (Adjust p-value)",title="shCLIC1_vs_shScr",size=16)
dev.off()