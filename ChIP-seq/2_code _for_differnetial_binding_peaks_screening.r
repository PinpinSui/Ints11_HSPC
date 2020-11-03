library(DiffBind)
dbObj <- dba(sampleSheet="/path/to/peakinformation_samples.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
#Quality Control
pdf("WT_VS_KO_samples_PCA.pdf")
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
dev.off()
pdf("WT_VS_KO_samples_clustering.pdf")
plot(dbObj)
dev.off()

###########Differential analysis
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
#  summary of results
dba.show(dbObj, bContrasts=T)
#  overlapping peaks identified by the two different tools (DESeq2 and edgeR)
pdf("WT_VS_KO_samples_plot_overlap.pdf")
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dev.off()
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
###################Results saving
# EdgeR
out <- as.data.frame(comp1.edgeR)
write.table(out, file="WT_VS_KO_samples_edgeR.txt", sep="\t", quote=F, col.names = NA)
# DESeq2
out <- as.data.frame(comp1.deseq)
write.table(out, file="WT_VS_KO_samples_deseq2.txt", sep="\t", quote=F, col.names = NA)



######################Significant differential results saving with bed format
# Create bed files for each keeping only significant peaks (p < 0.05)
# EdgeR
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed,file="WT_VS_KO_samples_edgeR_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="WT_VS_KO_samples_deseq2_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

######################annotation
library(ChIPpeakAnno)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(clusterProfiler)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
files = "WT_VS_KO_samples_deseq2_sig.bed"
peakAnno <- annotatePeak(files,tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Mm.eg.db")
final_annotation<-as.data.frame(peakAnno)
write.table(final_annotation,"WT_VS_KO_samples_deseq2_sig.txt",quote=FALSE,sep="'\t")
pdf("WT_VS_KO_Diff_peaks_deseq2_annotations_samples.pdf")
plotAnnoPie(peakAnno)
dev.off()

files = "WT_VS_KO_samples_edgeR_sig.bed"
peakAnno <- annotatePeak(files,tssRegion=c(-3000, 3000),TxDb=txdb,annoDb="org.Mm.eg.db")
final_annotation<-as.data.frame(peakAnno)
write.table(final_annotation,"WT_VS_KO_samples_EDGER_sig.txt",quote=FALSE,sep="'\t")
pdf("WT_VS_KO_Diff_peaks_EDGER_annotations_samples.pdf")
plotAnnoPie(peakAnno)
dev.off()
