library(Seurat)
library(cowplot)
library(magrittr)
library(dplyr)
library(sctransform)

wt.data <- Read10X(data.dir = "/path/to/WT/filtered_feature_bc_matrix/")
ko.data <- Read10X(data.dir = "/path/to/KO/filtered_feature_bc_matrix/")


############### Set up WT object
wt <- CreateSeuratObject(counts = wt.data, project = "WT_INTS11", min.cells = 3)
#Doublets marking and removing
doublets<-read.table("/path/to/WT_doublet.txt",head=TRUE,sep=",")
wt$doublet<-doublets[,3]
wt <- subset(x = wt, subset = doublet=="False")
wt[["percent.mt"]] <- PercentageFeatureSet(object = wt, pattern = "^mt-")
wt[["percent.rp"]] <- PercentageFeatureSet(object = wt, pattern = "^Rp[sl]")

pdf("WT_VlnPlot_gene_UMI_mito.pdf")
VlnPlot(object = wt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
dev.off()

pdf("WT_visualize feature-feature relationships.pdf")
plot1 <- FeatureScatter(object = wt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = wt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

2.Cell filtering
wt <- subset(x = wt, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
wt$Genotype <- "WT"

##############3.sctransform for WT
wt <- SCTransform(wt, vars.to.regress = "percent.mt", verbose = TRUE)




############### Set up KO object
ko <- CreateSeuratObject(counts = ko.data, project = "KO_INTS11", min.cells = 3)
doublets<-read.table("/path/to/KO_doublet.txt",head=TRUE,sep=",")
ko$doublet<-doublets[,3]
ko <- subset(x = ko, subset = doublet=="False")
ko[["percent.mt"]] <- PercentageFeatureSet(object = ko, pattern = "^mt-")
ko[["percent.rp"]] <- PercentageFeatureSet(object = ko, pattern = "^Rp[sl]")
pdf("KO_VlnPlot_gene_UMI_mito.pdf")
VlnPlot(object = ko, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
dev.off()

pdf("KO_visualize feature-feature relationships.pdf")
plot1 <- FeatureScatter(object = ko, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = ko, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

ko <- subset(x = ko, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
ko$Genotype <- "KO"

##############3.sctransform for KO
ko <- SCTransform(ko, vars.to.regress = "percent.mt", verbose = TRUE)
###############################################
#############################integration
INTS11.list<-list(wt,ko)
INTS11.features <- SelectIntegrationFeatures(object.list = INTS11.list, nfeatures = 3000)
INTS11.list <- PrepSCTIntegration(object.list = INTS11.list, anchor.features = INTS11.features, verbose = TRUE)

INTS11.anchors <- FindIntegrationAnchors(object.list = INTS11.list, normalization.method = "SCT", anchor.features = INTS11.features, verbose = TRUE)
INTS11.combined <- IntegrateData(anchorset = INTS11.anchors, normalization.method = "SCT", verbose = TRUE)

INTS11.combined <- RunPCA(INTS11.combined, verbose = TRUE)
pdf("figure00_VizDimLoadings.pdf")
VizDimLoadings(object = INTS11.combined, dims = 1:4, reduction = "pca")
dev.off()
pdf("figure01_Dimplot.pdf")
DimPlot(object = INTS11.combined, reduction = "pca")
dev.off()
INTS11.combined <- JackStraw(object = INTS11.combined,num.replicate = 100,dims=50,assay="integrated")
INTS11.combined <- ScoreJackStraw(object = INTS11.combined, dims = 1:50)
pdf("figure02_PCA_summray.pdf")
JackStrawPlot(object =INTS11.combined , dims = 1:50)
dev.off()
pdf("figure03_PCA_elbow.pdf")
ElbowPlot(object = INTS11.combined,ndims=50)
dev.off()

pdf("figure04_DimHeatmap_dims_1.pdf")
DimHeatmap(object = INTS11.combined, dims = 1, cells = 500, balanced = TRUE)
dev.off()
pdf("figure05_DimHeatmap_dims_1_20.pdf")
DimHeatmap(object = INTS11.combined, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()
pdf("figure05_DimHeatmap_dims_21_40.pdf")
DimHeatmap(object = INTS11.combined, dims = 21:40, cells = 500, balanced = TRUE)
dev.off()
pdf("figure05_DimHeatmap_dims_41_50.pdf")
DimHeatmap(object = INTS11.combined, dims = 41:50, cells = 500, balanced = TRUE)
dev.off()


# t-SNE and Clustering
INTS11.combined <- RunUMAP(object = INTS11.combined,, dims = 1:16)
INTS11.combined <- FindNeighbors(object = INTS11.combined, reduction = "pca", dims = 1:16)
INTS11.combined <- FindClusters(INTS11.combined, verbose=TRUE,resolution = 1.2)
DimPlot(INTS11.combined, label = TRUE)

#INTS11.combined <- FindClusters(INTS11.combined, resolution = 1.5)
# Visualization
p1 <- DimPlot(object = INTS11.combined, reduction = "umap", group.by = "Genotype")
p2 <- DimPlot(object = INTS11.combined, reduction = "umap", label = TRUE)
pdf("figure1_overall_landscape.pdf",height=8,width=12)
plot_grid(p1, p2)
dev.off()
#To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster
pdf("figure2_visualize the two conditions side-by-side.pdf",height=8,width=12)
DimPlot(object = INTS11.combined, reduction = "umap", split.by = "Genotype", label = TRUE)
dev.off()
#overall landscape
pdf("figure3_overall landscape of all cell clustering.pdf")
DimPlot(object = INTS11.combined, label = TRUE)
dev.off()
pdf("figure3_overall landscape of all cell clustering_new_size.pdf",height=8,width=8)
DimPlot(object = INTS11.combined, label = TRUE)
dev.off()
#############save temp result
DefaultAssay(object = INTS11.combined) <- "RNA"
INTS11.combined <- NormalizeData(INTS11.combined, verbose = TRUE)
#################markers selection
INTS11.combined.markers <- FindAllMarkers(object = INTS11.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,assay="RNA",slot="data")
write.table(INTS11.combined.markers,"table1_INTS11.combined.markers.txt",sep="\t",quote=FALSE)
write.table(INTS11.combined@meta.data,"table2_all_data.txt",sep="\t",quote=FALSE)

####################################################annotation############################################################
INTS11.combined <- RenameIdents(INTS11.combined, `0` = "preNeu", `1` = "preNeu", `2` = "HSC", `3` = "proNeu", 
`4` = "MonoP", `5` = "proNeu", `6` = "Ery/MK", `7` = "preNeu",  `8` = "preNeu",
`9` = "proNeu", `10` = "preNeu", `11` = "preNeu", `12` = "proNeu", `13` = "preNeu",
`14` = "preNeu", `15` = "immNeu", `16` = "immNeu", `17` = "Bas/Mast", `18` = "LMPP",
`19` = "immNeu", `20` = "Eosi", `21` = "DC", `22` = "Macro", `23` = "NK", `24` = "B cell", `25` = "Bas/Mast")
#To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster
p1 <- DimPlot(object = INTS11.combined, reduction = "umap", group.by = "Genotype")
p2 <- DimPlot(object = INTS11.combined, reduction = "umap", label = TRUE,label.size=5)
pdf("After_with_B_figure1_overall_landscape.pdf",height=8,width=12)
plot_grid(p1, p2)
dev.off()
pdf("After_with_B_figure2_visualize the two conditions side-by-side.pdf",height=8,width=12)
DimPlot(object = INTS11.combined, split.by = "Genotype", label = TRUE,label.size=5)
dev.off()
#overall landscape
pdf("After_with_B_figure3_overall landscape of all cell clustering.pdf")
DimPlot(object = INTS11.combined, label = TRUE,label.size=5)
dev.off()
pdf("After_with_B_figure3_overall landscape of all cell clustering_new_size.pdf",height=8,width=8)
DimPlot(object = INTS11.combined, label = TRUE,label.size=5)
dev.off()

INTS11.combined <- subset(INTS11.combined, idents = c("HSC","LMPP","MonoP","Ery/MK","proNeu","preNeu","immNeu","Macro","Bas/Mast","Eosi","DC","NK"))
p1 <- DimPlot(object = INTS11.combined, reduction = "umap", group.by = "Genotype")
p2 <- DimPlot(object = INTS11.combined, reduction = "umap", label = TRUE,label.size=5)
pdf("After_remove_B_figure1_overall_landscape.pdf",height=8,width=12)
plot_grid(p1, p2)
dev.off()
pdf("After_remove_B_figure2_visualize the two conditions side-by-side.pdf",height=8,width=12)
DimPlot(object = INTS11.combined, split.by = "Genotype", label = TRUE,label.size=5)
dev.off()
#overall landscape
pdf("After_remove_B_figure3_overall landscape of all cell clustering.pdf")
DimPlot(object = INTS11.combined, label = TRUE,label.size=5)
dev.off()
pdf("After_remove_B_figure3_overall landscape of all cell clustering_new_size.pdf",height=8,width=8)
DimPlot(object = INTS11.combined, label = TRUE,label.size=5)
dev.off()


saveRDS(INTS11.combined, file = "0_initiated_dim_16_res_1.2_RDS_after_annotation.rds")
#INTS11.combined.markers <- FindAllMarkers(object = INTS11.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,assay="RNA",slot="data")
#write.table(INTS11.combined.markers,"after_table1_INTS11.combined.markers.txt",sep="\t",quote=FALSE)
write.table(INTS11.combined@meta.data,"after_table2_all_data.txt",sep="\t",quote=FALSE)
#write.table(Idents(INTS11.combined),"after_table3_ident.txt",sep="\t",quote=FALSE)




################Reclustering for HSC population
HSC_subset <- subset(INTS11.combined, idents = "HSC")
DefaultAssay(HSC_subset)<-"integrated"
HSC_subset <- RunPCA(HSC_subset, verbose = TRUE)
pdf("figure00_VizDimLoadings.pdf")
VizDimLoadings(object = HSC_subset, dims = 1:4, reduction = "pca")
dev.off()
pdf("figure01_Dimplot.pdf")
DimPlot(object = HSC_subset, reduction = "pca")
dev.off()
HSC_subset <- JackStraw(object = HSC_subset,num.replicate = 100,dims=50,assay="integrated")
HSC_subset <- ScoreJackStraw(object = HSC_subset, dims = 1:50)
pdf("figure02_PCA_summray.pdf")
JackStrawPlot(object =HSC_subset , dims = 1:50)
dev.off()
pdf("figure03_PCA_elbow.pdf")
ElbowPlot(object = HSC_subset,ndims=50)
dev.off()

pdf("figure04_DimHeatmap_dims_1.pdf")
DimHeatmap(object = HSC_subset, dims = 1, cells = 500, balanced = TRUE)
dev.off()
pdf("figure05_DimHeatmap_dims_1_20.pdf")
DimHeatmap(object = HSC_subset, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()
pdf("figure05_DimHeatmap_dims_21_40.pdf")
DimHeatmap(object = HSC_subset, dims = 21:40, cells = 500, balanced = TRUE)
dev.off()
pdf("figure05_DimHeatmap_dims_41_50.pdf")
DimHeatmap(object = HSC_subset, dims = 41:50, cells = 500, balanced = TRUE)
dev.off()


HSC_subset <- RunUMAP(object = HSC_subset, dims = 1:10)
HSC_subset <- FindNeighbors(object = HSC_subset, reduction = "pca", dims = 1:10)
HSC_subset <- FindClusters(HSC_subset, verbose=TRUE,resolution = 0.5)
DimPlot(HSC_subset, label = TRUE)
p1 <- DimPlot(object = HSC_subset, reduction = "umap", group.by = "Genotype")
p2 <- DimPlot(object = HSC_subset, reduction = "umap", label = TRUE)
pdf("figure1_overall_landscape.pdf",height=8,width=12)
plot_grid(p1, p2)
dev.off()
#To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster
pdf("figure2_visualize the two conditions side-by-side.pdf",height=8,width=12)
DimPlot(object = HSC_subset, reduction = "umap", split.by = "Genotype", label = TRUE)
dev.off()
#overall landscape
pdf("figure3_overall landscape of all cell clustering.pdf")
DimPlot(object = HSC_subset, label = TRUE)
dev.off()
HSC_subset.markers <- FindAllMarkers(object = HSC_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,assay="RNA",slot="data")
write.table(HSC_subset.markers,"HSC_subset_table1_INTS11.combined.markers.txt",sep="\t",quote=FALSE)
write.table(HSC_subset@meta.data,"HSC_subset_table2_all_data.txt",sep="\t",quote=FALSE)


HSC_subset <- RenameIdents(HSC_subset, `0` = "ST-HSC", `1` = "LT-HSC", `2` = "MPP3", `3` = "MPP2")
#To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster
p1 <- DimPlot(object = HSC_subset, reduction = "umap", group.by = "Genotype")
p2 <- DimPlot(object = HSC_subset, reduction = "umap", label = TRUE,label.size=5)
pdf("After_figure1_overall_landscape.pdf",height=8,width=12)
plot_grid(p1, p2)
dev.off()
pdf("After_figure2_visualize the two conditions side-by-side.pdf",height=8,width=12)
DimPlot(object = HSC_subset, split.by = "Genotype", label = TRUE,label.size=5)
dev.off()
#overall landscape
pdf("After_figure3_overall landscape of all cell clustering.pdf")
DimPlot(object = HSC_subset, label = TRUE,label.size=5)
dev.off()
pdf("After_figure3_overall landscape of all cell clustering_new_size.pdf",height=8,width=8)
DimPlot(object = HSC_subset, label = TRUE,label.size=5)
dev.off()

HSC_subset.markers <- FindAllMarkers(object = HSC_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,assay="RNA",slot="data")
write.table(HSC_subset.markers,"after_table1_HSC_subset.markers.txt",sep="\t",quote=FALSE)
write.table(HSC_subset@meta.data,"after_table2_all_data.txt",sep="\t",quote=FALSE)
write.table(Idents(HSC_subset),"after_table3_ident.txt",sep="\t",quote=FALSE)


INTS11.combined$umap@cell.embeddings[,1]<- -(INTS11.combined$umap@cell.embeddings[,1])
INTS11.combined$umap@cell.embeddings[,2]<- -(INTS11.combined$umap@cell.embeddings[,2])

#To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster
p1 <- DimPlot(object = INTS11.combined, reduction = "umap", group.by = "Genotype")
p2 <- DimPlot(object = INTS11.combined, reduction = "umap", label = TRUE,label.size=5)
pdf("Reorder_After_figure1_overall_landscape.pdf",height=8,width=12)
plot_grid(p1, p2)
dev.off()
pdf("Reorder_After_figure2_visualize the two conditions side-by-side.pdf",height=8,width=12)
DimPlot(object = INTS11.combined, split.by = "Genotype", label = TRUE,label.size=5)
dev.off()
#overall landscape
pdf("Reorder_After_figure3_overall landscape of all cell clustering.pdf")
DimPlot(object = INTS11.combined, label = TRUE,label.size=5)
dev.off()
pdf("Reorder_After_figure3_overall landscape of all cell clustering_new_size.pdf",height=8,width=8)
DimPlot(object = INTS11.combined, label = TRUE,label.size=5)
dev.off()



HSC_subset$umap@cell.embeddings[,1]<- -(HSC_subset$umap@cell.embeddings[,1])
HSC_subset$umap@cell.embeddings[,2]<- -(HSC_subset$umap@cell.embeddings[,2])

#To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster
p1 <- DimPlot(object = HSC_subset, reduction = "umap", group.by = "Genotype")
p2 <- DimPlot(object = HSC_subset, reduction = "umap", label = TRUE,label.size=5)
pdf("Reorder_After_figure1_overall_landscape.pdf",height=8,width=12)
plot_grid(p1, p2)
dev.off()
pdf("Reorder_After_figure2_visualize the two conditions side-by-side.pdf",height=8,width=12)
DimPlot(object = HSC_subset, split.by = "Genotype", label = TRUE,label.size=5)
dev.off()
#overall landscape
pdf("Reorder_After_figure3_overall landscape of all cell clustering.pdf")
DimPlot(object = HSC_subset, label = TRUE,label.size=5)
dev.off()
pdf("Reorder_After_figure3_overall landscape of all cell clustering_new_size.pdf",height=8,width=8)
DimPlot(object = HSC_subset, label = TRUE,label.size=5)
dev.off()



######################GSEA
INTS11.combined$celltype.Genotype <- paste(Idents(INTS11.combined), INTS11.combined$Genotype, sep = "_")
INTS11.combined$celltype <- Idents(INTS11.combined)
Idents(INTS11.combined) <- INTS11.combined$celltype.Genotype


HSC_subset$celltype.Genotype <- paste(Idents(HSC_subset), HSC_subset$Genotype, sep = "_")
HSC_subset$celltype <- Idents(HSC_subset)
Idents(HSC_subset) <- HSC_subset$celltype.Genotype


#################################################################Generating datasets for GSEA
#################################overall
all_cell_name_overall<-rownames(INTS11.combined@meta.data)
expression_matrix_overall<-as.matrix(INTS11.combined[["RNA"]]@data)
tracking_id_overall<-rownames(expression_matrix_overall)
Description_overall<-rep("na",length(tracking_id_overall))

position_HSC_KO<-which(INTS11.combined@meta.data[,12]=="HSC_KO")
position_HSC_WT<-which(INTS11.combined@meta.data[,12]=="HSC_WT")
position_B_KO<-which(INTS11.combined@meta.data[,12]=="B cell_KO")
position_B_WT<-which(INTS11.combined@meta.data[,12]=="B cell_WT")
position_NK_KO<-which(INTS11.combined@meta.data[,12]=="NK_KO")
position_NK_WT<-which(INTS11.combined@meta.data[,12]=="NK_WT")
position_LMPP_KO<-which(INTS11.combined@meta.data[,12]=="LMPP_KO")
position_LMPP_WT<-which(INTS11.combined@meta.data[,12]=="LMPP_WT")
position_DC_KO<-which(INTS11.combined@meta.data[,12]=="DC_KO")
position_DC_WT<-which(INTS11.combined@meta.data[,12]=="DC_WT")
position_Mono_KO<-which(INTS11.combined@meta.data[,12]=="Mono_KO")
position_Mono_WT<-which(INTS11.combined@meta.data[,12]=="Mono_WT")
position_Macro_KO<-which(INTS11.combined@meta.data[,12]=="Macro_KO")
position_Macro_WT<-which(INTS11.combined@meta.data[,12]=="Macro_WT")
position_Bas_Mast_KO<-which(INTS11.combined@meta.data[,12]=="Bas/Mast_KO")
position_Bas_Mast_WT<-which(INTS11.combined@meta.data[,12]=="Bas/Mast_WT")
position_Ery_Mega_KO<-which(INTS11.combined@meta.data[,12]=="Ery/Mega_KO")
position_Ery_Mega_WT<-which(INTS11.combined@meta.data[,12]=="Ery/Mega_WT")
position_proNeu_KO<-which(INTS11.combined@meta.data[,12]=="proNeu_KO")
position_proNeu_WT<-which(INTS11.combined@meta.data[,12]=="proNeu_WT")
position_preNeu_KO<-which(INTS11.combined@meta.data[,12]=="preNeu_KO")
position_preNeu_WT<-which(INTS11.combined@meta.data[,12]=="preNeu_WT")
position_immNeu_KO<-which(INTS11.combined@meta.data[,12]=="immNeu_KO")
position_immNeu_WT<-which(INTS11.combined@meta.data[,12]=="immNeu_WT")
position_Eosi_KO<-which(INTS11.combined@meta.data[,12]=="Eosi_KO")
position_Eosi_WT<-which(INTS11.combined@meta.data[,12]=="Eosi_WT")

HSC_KO_cell<-all_cell_name_overall[position_HSC_KO]
HSC_WT_cell<-all_cell_name_overall[position_HSC_WT]
HSC_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,HSC_KO_cell],expression_matrix_overall[,HSC_WT_cell])
write.table(HSC_expression,"./GSEA/GSEA_matrix_HSC_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_HSC<-t(as.matrix(c(rep("KO",length(HSC_KO_cell)),rep("WT",length(HSC_WT_cell)))))
write.table(label_HSC,"./GSEA/GSEA_label_HSC_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


B_KO_cell<-all_cell_name_overall[position_B_KO]
B_WT_cell<-all_cell_name_overall[position_B_WT]
B_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,B_KO_cell],expression_matrix_overall[,B_WT_cell])
write.table(B_expression,"./GSEA/GSEA_matrix_B_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_B<-t(as.matrix(c(rep("KO",length(B_KO_cell)),rep("WT",length(B_WT_cell)))))
write.table(label_B,"./GSEA/GSEA_label_B_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

NK_KO_cell<-all_cell_name_overall[position_NK_KO]
NK_WT_cell<-all_cell_name_overall[position_NK_WT]
NK_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,NK_KO_cell],expression_matrix_overall[,NK_WT_cell])
write.table(NK_expression,"./GSEA/GSEA_matrix_NK_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_NK<-t(as.matrix(c(rep("KO",length(NK_KO_cell)),rep("WT",length(NK_WT_cell)))))
write.table(label_NK,"./GSEA/GSEA_label_NK_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

DC_KO_cell<-all_cell_name_overall[position_DC_KO]
DC_WT_cell<-all_cell_name_overall[position_DC_WT]
DC_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,DC_KO_cell],expression_matrix_overall[,DC_WT_cell])
write.table(DC_expression,"./GSEA/GSEA_matrix_DC_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_DC<-t(as.matrix(c(rep("KO",length(DC_KO_cell)),rep("WT",length(DC_WT_cell)))))
write.table(label_DC,"./GSEA/GSEA_label_DC_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

Macro_KO_cell<-all_cell_name_overall[position_Macro_KO]
Macro_WT_cell<-all_cell_name_overall[position_Macro_WT]
Macro_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,Macro_KO_cell],expression_matrix_overall[,Macro_WT_cell])
write.table(Macro_expression,"./GSEA/GSEA_matrix_Macro_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_Macro<-t(as.matrix(c(rep("KO",length(Macro_KO_cell)),rep("WT",length(Macro_WT_cell)))))
write.table(label_Macro,"./GSEA/GSEA_label_Macro_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

Bas_Mast_KO_cell<-all_cell_name_overall[position_Bas_Mast_KO]
Bas_Mast_WT_cell<-all_cell_name_overall[position_Bas_Mast_WT]
Bas_Mast_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,Bas_Mast_KO_cell],expression_matrix_overall[,Bas_Mast_WT_cell])
write.table(Bas_Mast_expression,"./GSEA/GSEA_matrix_Bas_Mast_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_Bas_Mast<-t(as.matrix(c(rep("KO",length(Bas_Mast_KO_cell)),rep("WT",length(Bas_Mast_WT_cell)))))
write.table(label_Bas_Mast,"./GSEA/GSEA_label_Bas_Mast_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

LMPP_KO_cell<-all_cell_name_overall[position_LMPP_KO]
LMPP_WT_cell<-all_cell_name_overall[position_LMPP_WT]
LMPP_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,LMPP_KO_cell],expression_matrix_overall[,LMPP_WT_cell])
write.table(LMPP_expression,"./GSEA/GSEA_matrix_LMPP_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_LMPP<-t(as.matrix(c(rep("KO",length(LMPP_KO_cell)),rep("WT",length(LMPP_WT_cell)))))
write.table(label_LMPP,"./GSEA/GSEA_label_LMPP_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

Mono_KO_cell<-all_cell_name_overall[position_Mono_KO]
Mono_WT_cell<-all_cell_name_overall[position_Mono_WT]
Mono_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,Mono_KO_cell],expression_matrix_overall[,Mono_WT_cell])
write.table(Mono_expression,"./GSEA/GSEA_matrix_Mono_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_Mono<-t(as.matrix(c(rep("KO",length(Mono_KO_cell)),rep("WT",length(Mono_WT_cell)))))
write.table(label_Mono,"./GSEA/GSEA_label_Mono_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

Ery_Mega_KO_cell<-all_cell_name_overall[position_Ery_Mega_KO]
Ery_Mega_WT_cell<-all_cell_name_overall[position_Ery_Mega_WT]
Ery_Mega_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,Ery_Mega_KO_cell],expression_matrix_overall[,Ery_Mega_WT_cell])
write.table(Ery_Mega_expression,"./GSEA/GSEA_matrix_Ery_Mega_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_Ery_Mega<-t(as.matrix(c(rep("KO",length(Ery_Mega_KO_cell)),rep("WT",length(Ery_Mega_WT_cell)))))
write.table(label_Ery_Mega,"./GSEA/GSEA_label_Ery_Mega_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

Eosi_KO_cell<-all_cell_name_overall[position_Eosi_KO]
Eosi_WT_cell<-all_cell_name_overall[position_Eosi_WT]
Eosi_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,Eosi_KO_cell],expression_matrix_overall[,Eosi_WT_cell])
write.table(Eosi_expression,"./GSEA/GSEA_matrix_Eosi_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_Eosi<-t(as.matrix(c(rep("KO",length(Eosi_KO_cell)),rep("WT",length(Eosi_WT_cell)))))
write.table(label_Eosi,"./GSEA/GSEA_label_Eosi_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

immNeu_KO_cell<-all_cell_name_overall[position_immNeu_KO]
immNeu_WT_cell<-all_cell_name_overall[position_immNeu_WT]
immNeu_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,immNeu_KO_cell],expression_matrix_overall[,immNeu_WT_cell])
write.table(immNeu_expression,"./GSEA/GSEA_matrix_immNeu_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_immNeu<-t(as.matrix(c(rep("KO",length(immNeu_KO_cell)),rep("WT",length(immNeu_WT_cell)))))
write.table(label_immNeu,"./GSEA/GSEA_label_immNeu_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

preNeu_KO_cell<-all_cell_name_overall[position_preNeu_KO]
preNeu_WT_cell<-all_cell_name_overall[position_preNeu_WT]
preNeu_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,preNeu_KO_cell],expression_matrix_overall[,preNeu_WT_cell])
write.table(preNeu_expression,"./GSEA/GSEA_matrix_preNeu_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_preNeu<-t(as.matrix(c(rep("KO",length(preNeu_KO_cell)),rep("WT",length(preNeu_WT_cell)))))
write.table(label_preNeu,"./GSEA/GSEA_label_preNeu_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

proNeu_KO_cell<-all_cell_name_overall[position_proNeu_KO]
proNeu_WT_cell<-all_cell_name_overall[position_proNeu_WT]
proNeu_expression<-cbind(tracking_id_overall,Description_overall,expression_matrix_overall[,proNeu_KO_cell],expression_matrix_overall[,proNeu_WT_cell])
write.table(proNeu_expression,"./GSEA/GSEA_matrix_proNeu_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_proNeu<-t(as.matrix(c(rep("KO",length(proNeu_KO_cell)),rep("WT",length(proNeu_WT_cell)))))
write.table(label_proNeu,"./GSEA/GSEA_label_proNeu_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


#################################HSC_subset
all_cell_name_HSC_subset<-rownames(HSC_subset@meta.data)
expression_matrix_HSC_subset<-as.matrix(HSC_subset[["RNA"]]@data)
tracking_id_HSC_subset<-rownames(expression_matrix_HSC_subset)
Description_HSC_subset<-rep("na",length(tracking_id_HSC_subset))

position_LT_HSC_KO<-which(HSC_subset@meta.data[,13]=="LT-HSC_KO")
position_LT_HSC_WT<-which(HSC_subset@meta.data[,13]=="LT-HSC_WT")
position_ST_HSC_KO<-which(HSC_subset@meta.data[,13]=="ST-HSC_KO")
position_ST_HSC_WT<-which(HSC_subset@meta.data[,13]=="ST-HSC_WT")
position_MPP2_KO<-which(HSC_subset@meta.data[,13]=="MPP2_KO")
position_MPP2_WT<-which(HSC_subset@meta.data[,13]=="MPP2_WT")
position_MPP3_KO<-which(HSC_subset@meta.data[,13]=="MPP3_KO")
position_MPP3_WT<-which(HSC_subset@meta.data[,13]=="MPP3_WT")

LT_HSC_KO_cell<-all_cell_name_HSC_subset[position_LT_HSC_KO]
LT_HSC_WT_cell<-all_cell_name_HSC_subset[position_LT_HSC_WT]
LT_HSC_expression<-cbind(tracking_id_HSC_subset,Description_HSC_subset,expression_matrix_HSC_subset[,LT_HSC_KO_cell],expression_matrix_HSC_subset[,LT_HSC_WT_cell])
write.table(LT_HSC_expression,"./GSEA/GSEA_matrix_LT_HSC_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_LT_HSC<-t(as.matrix(c(rep("KO",length(LT_HSC_KO_cell)),rep("WT",length(LT_HSC_WT_cell)))))
write.table(label_LT_HSC,"./GSEA/GSEA_label_LT_HSC_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

ST_HSC_KO_cell<-all_cell_name_HSC_subset[position_ST_HSC_KO]
ST_HSC_WT_cell<-all_cell_name_HSC_subset[position_ST_HSC_WT]
ST_HSC_expression<-cbind(tracking_id_HSC_subset,Description_HSC_subset,expression_matrix_HSC_subset[,ST_HSC_KO_cell],expression_matrix_HSC_subset[,ST_HSC_WT_cell])
write.table(ST_HSC_expression,"./GSEA/GSEA_matrix_ST_HSC_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_ST_HSC<-t(as.matrix(c(rep("KO",length(ST_HSC_KO_cell)),rep("WT",length(ST_HSC_WT_cell)))))
write.table(label_ST_HSC,"./GSEA/GSEA_label_ST_HSC_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

MPP2_KO_cell<-all_cell_name_HSC_subset[position_MPP2_KO]
MPP2_WT_cell<-all_cell_name_HSC_subset[position_MPP2_WT]
MPP2_expression<-cbind(tracking_id_HSC_subset,Description_HSC_subset,expression_matrix_HSC_subset[,MPP2_KO_cell],expression_matrix_HSC_subset[,MPP2_WT_cell])
write.table(MPP2_expression,"./GSEA/GSEA_matrix_MPP2_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_MPP2<-t(as.matrix(c(rep("KO",length(MPP2_KO_cell)),rep("WT",length(MPP2_WT_cell)))))
write.table(label_MPP2,"./GSEA/GSEA_label_MPP2_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

MPP3_KO_cell<-all_cell_name_HSC_subset[position_MPP3_KO]
MPP3_WT_cell<-all_cell_name_HSC_subset[position_MPP3_WT]
MPP3_expression<-cbind(tracking_id_HSC_subset,Description_HSC_subset,expression_matrix_HSC_subset[,MPP3_KO_cell],expression_matrix_HSC_subset[,MPP3_WT_cell])
write.table(MPP3_expression,"./GSEA/GSEA_matrix_MPP3_expression_KO_WT.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
label_MPP3<-t(as.matrix(c(rep("KO",length(MPP3_KO_cell)),rep("WT",length(MPP3_WT_cell)))))
write.table(label_MPP3,"./GSEA/GSEA_label_MPP3_expression_KO_WT.cls",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


#####################Score calculating & visualization
#################merge 
INTS11_no_HSC<- subset(INTS11.combined, idents = c("MonoP","Ery/MK","proNeu","preNeu","immNeu","LMPP","Eosi","DC","Macro","NK","Bas/Mast"))
merged<-merge(x=INTS11_no_HSC,y=HSC_subset)
merged$celltype.Genotype <- paste(Idents(merged), merged$Genotype, sep = "_")
merged$celltype <- Idents(merged)
Idents(merged) <- merged$celltype.Genotype
Idents(merged) <- "celltype"
merged_HSPC<-subset(merged, idents = c("LT-HSC","ST-HSC","MPP2","MPP3","LMPP","Ery/MK","proNeu","MonoP"))
Idents(merged_HSPC) <- "celltype"
merged@active.ident<-factor(as.vector(merged@active.ident),levels=rev(c("NK","DC","Eosi","Bas/Mast","Macro","MonoP","immNeu","preNeu","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC")))
merged_HSPC@active.ident<-factor(merged_HSPC@active.ident,levels=rev(c("MonoP","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC")))
merged_HSPC@meta.data$Genotype<-factor(merged_HSPC@meta.data$Genotype, levels=c("WT","KO"))

#############Cdkn1a
pdf("cdkn1a_among_populations.pdf",width=6,height=5)
VlnPlot(merged_HSPC, features = c("Cdkn1a"), split.by = "Genotype",pt.size = 0.000, combine = FALSE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE)
dev.off()

####################addscore
##monocyte
monocyte<-list(as.matrix(read.table("monocyte_marker.txt"))[,1])
merged<-AddModuleScore(merged, features=monocyte, name="monocyte_score",assay="RNA")
merged@meta.data$Genotype<-factor(merged@meta.data$Genotype, levels=c("WT","KO"))
Idents(merged) <- "celltype"
merged@active.ident<-factor(merged@active.ident,levels=rev(c("NK","DC","Eosi","Bas/Mast","Macro","MonoP","immNeu","preNeu","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC")))
pdf("monocyte_in_different_population_split_WT_first_HSPC.pdf")
VlnPlot(merged, features = c("monocyte_score1"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("LMPP","MPP3","MPP2","ST-HSC","LT-HSC"))
dev.off()


proliferation<-list(as.matrix(read.table("proliferation.txt"))[,1])
merged<-AddModuleScore(merged, features=proliferation, name="proliferation_score",assay="RNA")
merged@meta.data$Genotype<-factor(merged@meta.data$Genotype, levels=c("WT","KO"))
Idents(merged) <- "celltype"
merged@active.ident<-factor(merged@active.ident,levels=rev(c("NK","DC","Eosi","Bas/Mast","Macro","MonoP","immNeu","preNeu","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC")))
pdf("proliferation_in_different_population_split_WT_first_HSPC.pdf")
VlnPlot(merged, features = c("proliferation_score1"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("LMPP","MPP3","MPP2","ST-HSC","LT-HSC"))
dev.off()


#############################cell cycle
s.genes<-as.matrix(read.table("cellcycle_s_genes.txt",head=FALSE))[,1]
g2m.genes<-as.matrix(read.table("cellcycle_g2m_genes.txt",head=FALSE))[,1]
merged <- CellCycleScoring(merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE,assay="RNA",slot="data")


Ery_MK<-list(as.matrix(read.table("Ery_MK_genelist_from_ZP_reduced.txt"))[,1])
merged<-AddModuleScore(merged, features=Ery_MK, name="Ery_MK_score",assay="RNA")
merged@meta.data$Genotype<-factor(merged@meta.data$Genotype, levels=c("WT","KO"))
Idents(merged) <- "celltype"
merged@active.ident<-factor(merged@active.ident,levels=rev(c("NK","DC","Eosi","Bas/Mast","Macro","MonoP","immNeu","preNeu","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC")))
pdf("0_Ery_MK_in_different_population_split_WT_first_HSPC.pdf")
VlnPlot(merged, features = c("Ery_MK_score1"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("LMPP","MPP3","MPP2","ST-HSC","LT-HSC"))
dev.off()



stemness<-list(as.matrix(read.table("stemness_gene.txt"))[,1])
merged<-AddModuleScore(merged, features=stemness, name="stemness_score",assay="RNA")
merged@meta.data$Genotype<-factor(merged@meta.data$Genotype, levels=c("WT","KO"))
Idents(merged) <- "celltype"
merged@active.ident<-factor(merged@active.ident,levels=rev(c("NK","DC","Eosi","Bas/Mast","Macro","MonoP","immNeu","preNeu","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC")))
pdf("0_stemness_in_different_population_split_WT_first_HSPC.pdf")
VlnPlot(merged, features = c("stemness_score1"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("LMPP","MPP3","MPP2","ST-HSC","LT-HSC"))
dev.off()



apoptosis<-list(as.matrix(read.table("apoptosis_gene.txt"))[,1])
merged<-AddModuleScore(merged, features=apoptosis, name="apoptosis_score",assay="RNA")
merged@meta.data$Genotype<-factor(merged@meta.data$Genotype, levels=c("WT","KO"))
Idents(merged) <- "celltype"
merged@active.ident<-factor(merged@active.ident,levels=rev(c("NK","DC","Eosi","Bas/Mast","Macro","MonoP","immNeu","preNeu","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC")))
pdf("apoptosis_in_different_population_split_WT_first_HSPC.pdf")
VlnPlot(merged, features = c("apoptosis_score1"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("LMPP","MPP3","MPP2","ST-HSC","LT-HSC"))
dev.off()
write.table(merged@meta.data,"merged.table.txt",quote=FALSE,sep="\t")



pdf("proliferation_in_different_population_split_WT_first_candidate_populations.pdf")
VlnPlot(merged, features = c("proliferation_score1"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("MonoP","proNeu","Ery/MK","LMPP","MPP3","ST-HSC","LT-HSC"))
dev.off()


pdf("stemnes_in_different_population_split_WT_first_candidate_populations.pdf")
VlnPlot(merged, features = c("stemness_score1"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("MonoP","proNeu","Ery/MK","LMPP","MPP3","ST-HSC","LT-HSC"))
dev.off()

pdf("apoptosis_in_different_population_split_WT_first_candidate_populations.pdf")
VlnPlot(merged, features = c("apoptosis_score1"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("MonoP","proNeu","Ery/MK","LMPP","MPP3","ST-HSC","LT-HSC"))
dev.off()


###########candidate genes
pdf("candidate_genes_in_different_population_split_WT_first_candidate_populations.pdf",width=8,height=8)
VlnPlot(merged, features = c("Gata2","Id2","Perp","Plk2","Trib1"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,ncol=2,
idents=c("MonoP","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC"))
dev.off()

pdf("candidate_genes_in_different_population_split_WT_first_candidate_populations_Gata2.pdf",width=6,height=4)
VlnPlot(merged, features = c("Gata2"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("MonoP","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC"))
dev.off()



pdf("candidate_genes_in_different_population_split_WT_first_candidate_populations_Id2.pdf",width=6,height=4)
VlnPlot(merged, features = c("Id2"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("MonoP","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC"))
dev.off()

pdf("candidate_genes_in_different_population_split_WT_first_candidate_populations_Perp.pdf",width=6,height=4)
VlnPlot(merged, features = c("Perp"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("MonoP","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC"))
dev.off()

pdf("candidate_genes_in_different_population_split_WT_first_candidate_populations_Plk2.pdf",width=6,height=4)
VlnPlot(merged, features = c("Plk2"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("MonoP","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC"))
dev.off()

pdf("candidate_genes_in_different_population_split_WT_first_candidate_populations_Trib1.pdf",width=6,height=4)
VlnPlot(merged, features = c("Trib1"), split.by = "Genotype",pt.size = 0.000, combine = TRUE,
assay="RNA",slot="data",cols=c("firebrick3","gray"),sort=FALSE,
idents=c("MonoP","proNeu","Ery/MK","LMPP","MPP3","MPP2","ST-HSC","LT-HSC"))
dev.off()

#end


