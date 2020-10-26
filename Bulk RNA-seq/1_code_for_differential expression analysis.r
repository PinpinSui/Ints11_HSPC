#Step 1 Data normalizing and DEGs Identifying
library(DESeq2)
x<-as.matrix(read.table("raw_count.txt",row.names=1,header=TRUE))
x_filter<-c()
gene_name<-rownames(x)
gene_name_filter<-c()
for( i in 1:nrow(x)){
	if(sum(x[i,]==0)<5){
		x_filter<-rbind(x_filter,x[i,])
		gene_name_filter<-c(gene_name_filter,gene_name[i])
		}
}
x_filter<-as.matrix(x_filter)
rownames(x_filter)<-gene_name_filter
DE_file<-x_filter
#countData <- all
colData <- data.frame(row.names=colnames(DE_file),condition=factor(c("KO","KO","KO","WT","WT","WT"),levels=c("KO","WT")))
dds <- DESeqDataSetFromMatrix(countData = DE_file,colData = colData,design = ~ condition)
dds <- DESeq(dds)


res <- results(dds,contrast=c("condition","KO","WT"))
write.table(res,file="DEseq2_differential_genes_WT_KO.txt",sep="\t",na="NA",quote=TRUE)
pdf("DESeq2_MA_WT_KO.pdf")
plotMA(res)
dev.off()
pdf("DESeq2_MA_WT_KO_new_ylim.pdf")
plotMA(res,ylim=c(-6,8))
dev.off()
pdf("DESeq2_hist_WT_KO.pdf")
hist(res$pvalue, breaks=100, col="skyblue", border="slateblue", main="",xlab="P-value")
dev.off()
resOrdered <- res[order(res$padj),]
head(resOrdered)


dds<-estimateSizeFactors(dds)
dds_normalized<-counts(dds,normalized=TRUE)
write.table(dds_normalized,"Normalize_count_DEGseq2.txt",sep="\t",quote=FALSE)


vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
############PCAplot
pdf("PCA_plot.pdf")
plotPCA(vsd, intgroup=c("condition"))
dev.off()



sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library(pheatmap)
pdf("pheatmap_distance_between_samples.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,cellwidth=30,cellheight=30)
dev.off()




#Step 2 Functional enrichment analysis for DEGs
library(clusterProfiler)
library(org.Mm.eg.db)
#GO
data_A_up<-as.matrix(read.table("differential_expression_gene_up.txt",head=FALSE,sep="\t"))
data_A_down<-as.matrix(read.table("differential_expression_gene_down.txt",head=FALSE,sep="\t"))

go_A_up <- enrichGO(data_A_up[,1], OrgDb = org.Mm.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2,keyType = 'SYMBOL')
go_A_down <- enrichGO(data_A_down[,1], OrgDb = org.Mm.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2,keyType = 'SYMBOL')

pdf("cluster_profile_up_GO_annotation_barplot.pdf",width=10)
barplot(go_A_up,showCategory=30,drop=T)
dev.off()
pdf("cluster_profile_up_GO_annotation_dotplot.pdf",width=10)
dotplot(go_A_up,showCategory=30)
dev.off()

pdf("cluster_profile_down_GO_annotation_barplot.pdf",width=10)
barplot(go_A_down,showCategory=30,drop=T)
dev.off()
pdf("cluster_profile_down_GO_annotation_dotplot.pdf",width=10)
dotplot(go_A_down,showCategory=30)
dev.off()



#KEGG
data_A_up<-as.matrix(read.table("differential_expression_gene_up.txt",head=FALSE,sep="\t"))
data_A_down<-as.matrix(read.table("differential_expression_gene_down.txt",head=FALSE,sep="\t"))

genelist_A_up <- bitr(data_A_up[,1], fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
genelist_A_down <- bitr(data_A_down[,1], fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")


kegg_A_up<- enrichKEGG(genelist_A_up[,2], organism = 'mmu', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                     minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
kegg_A_down<- enrichKEGG(genelist_A_down[,2], organism = 'mmu', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                     minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)

pdf("cluster_profile_up_KEGG_annotation_barplot.pdf",width=10)					 
barplot(kegg_A_up, showCategory=30)
dev.off()
pdf("cluster_profile_down_KEGG_annotation_barplot.pdf",width=10)					 
barplot(kegg_A_down, showCategory=30)
dev.off()


pdf("cluster_profile_up_KEGG_annotation_dotplot.pdf",width=10)					 
dotplot(kegg_A_up, showCategory=30)
dev.off()
pdf("cluster_profile_down_KEGG_annotation_dotplot.pdf",width=10)					 
dotplot(kegg_A_down, showCategory=30)
dev.off()





