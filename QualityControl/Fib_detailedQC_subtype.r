suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(ggthemr)
  library(patchwork)
  library(Matrix)
  library(stats)
  
  library(Seurat)
})




#############################################
################# preparation ###############
#############################################
set.seed(1)
wkdir = "/home/mztao/scRNA_analysis/data_analysis/HPC/results/"
setwd(wkdir)
source("/home/mztao/scRNA_analysis/data_analysis/HPC/scripts/common_functions.r")
SampleNames <- c("T0_CD45", "T0_tumor", "T1_1", "T1_2", "T2_1", "T3_1", "T3_2", "T4_1", "T4_2", "NT4", "T5_1", "T5_2", "NT5", "T6_1", "T6_2", "T7")

## load data
current.object <- readRDS(paste0(wkdir,"/step1_MergeMainTypes/AllSamples.RDS"))
dim(current.object);table(current.object$maintypes)







#############################################
# primary analysis on Endo cells 
#############################################
dir.create("./step2_6_Fib/")
setwd("./step2_6_Fib/")

## load data
sub.object <- subset(current.object, subset = maintypes == "Fibroblast");dim(sub.object)
sub.name <- "AllSample_Fib"
rm(current.object)
table(sub.object$orig.ident);table(sub.object$sample.patients)

## basic steps
# exclusing some genes expressed less than min cells
sub.object <- ObjectUpdate(sub.object, min.cells = 10);dim(sub.object)
## seurat pipeline
sub.object <- NormalizeData(sub.object, normalization.method = "LogNormalize", scale.factor = 10000)
sub.object <- FindVariableFeatures(sub.object, selection.method = "vst", nfeatures = max(round(dim(sub.object)[1] * 0.1), 2000))
sub.object <- ScaleData(sub.object, features = rownames(sub.object))

# checking cell cycle effect
dir.create("./CellCycle-Effect/")
sub.object <- CellCycleScoring(sub.object, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE);table(sub.object$Phase)
sub.object <- RunPCA(sub.object, features = c(cc.genes$s.genes, cc.genes$g2m.genes))
ggsave(paste0("./CellCycle-Effect/", sub.name, "_BeforeCC.jpg"), DimPlot(sub.object), dpi = 300, width = 10, height = 8)
#################################################################
sub.object <- RunPCA(sub.object, features = VariableFeatures(object = sub.object))
ggsave(paste0("./", sub.name, "_ElbowPlot.jpg"),ElbowPlot(sub.object,ndims = 40), dpi = 300)
sub.object <- RunUMAP(sub.object, dims = 1:20)
sub.object <- RunTSNE(sub.object, dims = 1:20)
sub.object <- FindNeighbors(sub.object, dims = 1:20)
sub.object <- FindClusters(sub.object, resolution = 0.8) # 分别度可以大一些方便做QC


p1 <- DimPlot(sub.object, reduction = "umap", label = TRUE)
p2 <- DimPlot(sub.object, reduction = "tsne", label = TRUE)
ggsave(paste0("./", sub.name, "_Reduction_raw.jpg"), p1|p2, dpi = 300, width = 12,height = 5)

table(sub.object$orig.ident);table(sub.object$seurat_clusters)


saveRDS(sub.object, file = paste0("./", sub.name, ".RDS"))









###########################################
# advanced QC in all Fib cells
###########################################
##### basic ones
# 1. plot
filepath = paste("./", sub.name, "_BasicStats.jpg", sep='')
ggsave(filepath, FeatureScatter(sub.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"), dpi = 300)

# 2. disssociation effect calculating for further use
DissGene.stats <- sub.object$percent.diss
OverRatio_DissGene <- length(DissGene.stats[DissGene.stats > 5.75 ]) / length(DissGene.stats)
pdf(file=paste("./", sub.name, "_DissliHist1.pdf", sep=''))
hist(DissGene.stats, breaks = 200, xlim = c(0,15))
abline(v = 5.75, col = "red", lwd =2)
text(10,500,paste0(as.character(round(OverRatio_DissGene*100,2)),"%"))
dev.off()


filepath = paste0("./", sub.name, "_Umap_DissEffectionScores.jpg")
ggsave(filepath, FeaturePlot(sub.object,reduction = "umap", features = "percent.diss"), dpi = 300,width = 10,height = 8)

subset_diss <- rep("NoEffected",dim(sub.object)[2])
subset_diss[sub.object$percent.diss > 5.75] <- "Effected"
sub.object$diss_info <- subset_diss
filepath = paste0("./", sub.name, "_Umap_DissEffection.jpg")
ggsave(filepath, DimPlot(sub.object,reduction = "umap", group.by = 'diss_info', cols = c("red","gray")), dpi = 300,width = 10,height = 8)

# 3. proliferation effect checking
ProliGene.stats <- sub.object$percent.proli
pdf(file=paste("./", sub.name, "_ProliHist.pdf", sep=''))
hist(ProliGene.stats, breaks = 200, xlim = c(0,0.5))
dev.off()

filepath = paste("./", sub.name, "_ProliEffection_Umap1.jpg", sep='')
ggsave(filepath,FeaturePlot(sub.object, reduction = "umap", features = c("MKI67","NUSAP1","PLK1","CDC20","CDK1")), dpi = 300, width = 10,height = 8)
filepath = paste("./", sub.name, "_ProliEffection_Umap2.jpg", sep='')
ggsave(filepath,FeaturePlot(sub.object, reduction = "umap", features = c("CDKN3","CENPA","BIRC5","PCNA","CCNA2")), dpi = 300, width = 10,height = 8)

# 4. red cells
filepath = paste("./", sub.name, "_HB_Umap.jpg", sep='')
ggsave(filepath,FeaturePlot(sub.object, reduction = "umap", features = c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ"), ncol = 3), dpi = 300, width = 10,height = 8)

# 5.1 Srublets
dir.create("./Srublets/")
Out10x(sub.object, "./Srublets/", sub.name)

# 在terminal运行程序Srublets.py，用python包Srubelts计算。python Scrublets.py > ./log
sub.scores <- as.data.frame(read.table(paste0("./Srublets/",sub.name,"_doublet_scores.txt")))
rownames(sub.scores) <- colnames(sub.object) ; colnames(sub.scores) <- "scores"
sub.object$sub.srublets.scores <- sub.scores

# plotting for visualization
pdf(file=paste("./", sub.name, "_sub_srublets_scores.pdf", sep=''))
hist(sub.object$sub.srublets.scores, breaks = 200, xlim = c(0,1))
dev.off()

threshold <- 0.2
pdf(file=paste("./", sub.name, "_sub_srublets_scores_threshold.pdf", sep=''))
hist(sub.object$sub.srublets.scores, breaks = 200, xlim = c(0,1))
abline(v = threshold, col = "red", lwd =2)
text(0.6,600,paste0(as.character(round(length(which(sub.object$sub.srublets.scores > threshold)) / dim(sub.object)[2]*100, 2)),"%"))
dev.off()

filepath = paste("./", sub.name, "_Umap_sub_srublets_scores.jpg", sep='')
ggsave(filepath,FeaturePlot(sub.object, reduction = "umap",features = "sub.srublets.scores"), dpi = 300, width = 10, height = 8)

subset_Sdoublets <- rep("NoEffected",dim(sub.object)[2])
subset_Sdoublets[sub.object$sub.srublets.scores > threshold] <- "Effected"
sub.object$Sub_HighSrubeltsEffect <- as.factor(subset_Sdoublets)
filepath = paste0("./", sub.name, "_Umap_sub_srublets_scores_group.jpg")
ggsave(filepath,DimPlot(sub.object,reduction = "umap", group.by = 'Sub_HighSrubeltsEffect', cols = c("red","gray")), dpi = 300,width = 10,height = 8)

saveRDS(sub.object, file = paste0("./", sub.name, ".RDS"))

# 5.2 DoubletDecon
library(DoubletDecon)
dir.create("./DoubletDecon/")

Idents(sub.object) <- sub.object$seurat_clusters
newFiles <- Improved_Seurat_Pre_Process(sub.object, num_genes=30, write_files=FALSE)
DoubletDecon_results <- Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, groupsFile=newFiles$newGroupsFile, filename=sub.name,
                                           location="./DoubletDecon/",fullDataFile=NULL,removeCC=FALSE,species="hsa", rhop=0.6, write=FALSE,PMF=TRUE,
                                           useFull=FALSE, heatmap=FALSE,centroids=TRUE,num_doubs=100,only50=FALSE,min_uniq=4,nCores=10)

doublets <- DoubletDecon_results[["Final_doublets_groups"]]; singlets <- DoubletDecon_results[["Final_nondoublets_groups"]]
sub_DoubleDeconEffect <- data.frame(groups = c(rep("doublets", dim(doublets)[1]),rep("singlets", dim(singlets)[1])))
rownames(sub_DoubleDeconEffect) <- c(rownames(doublets), rownames(singlets))

table(sub_DoubleDeconEffect)
rownames(sub_DoubleDeconEffect) <- gsub("[.]","-",rownames(sub_DoubleDeconEffect))
sub.object$sub_DoubleDeconEffect <- sub_DoubleDeconEffect
filepath = paste("./", sub.name, "_Umap_total_DoubleDecon_group.jpg", sep='')
ggsave(filepath,DimPlot(sub.object,reduction = "umap", group.by = 'sub_DoubleDeconEffect', cols = c("red","gray")), dpi = 300,width = 10,height = 8)

saveRDS(sub.object, file = paste0("./", sub.name, ".RDS"))

## combine all qc results from total all cells
# 1. doublets algorithms
p1 <- DimPlot(sub.object, reduction = "tsne", group.by = 'seurat_clusters', label = T)
p2 <- DimPlot(sub.object, reduction = "tsne", group.by = 'sample.patients')
p3 <- DimPlot(sub.object, reduction = "tsne", group.by = 'orig.ident')
p4 <- DimPlot(sub.object, reduction = "umap", group.by = 'seurat_clusters', label = T)
p5 <- DimPlot(sub.object, reduction = "umap", group.by = 'sample.patients')
p6 <- DimPlot(sub.object, reduction = "umap", group.by = 'orig.ident')
p7 <- FeaturePlot(sub.object, reduction = "umap", features = "total.srublets.scores")
p8 <- DimPlot(sub.object, reduction = "umap", group.by = 'Total_HighSrubeltsEffect', cols = c("red","gray"))
p9 <- DimPlot(sub.object,reduction = "umap", group.by = 'Total_DoubleDeconEffect', cols = c("red","gray"))
p10 <- FeaturePlot(sub.object, reduction = "umap", features = "sub.srublets.scores")
p11 <- DimPlot(sub.object, reduction = "umap", group.by = 'Sub_HighSrubeltsEffect', cols = c("red","gray"))
p12 <- DimPlot(sub.object,reduction = "umap", group.by = 'sub_DoubleDeconEffect', cols = c("red","gray"))
p <- (p1|p2|p3) / (p4|p5|p6) / (p7|p8|p9) / (p10|p11|p12)
ggsave(paste0("./", sub.name, "_Umap_doublets_QC_summary.jpg"), p, dpi = 300, width = 17, height = 18)

# 2. markers, diss and mito effect
markers <- c("EPCAM","KRT19","KRT18","CLDN5","CDH5","COL1A1","COL13A1","PTPRC","MS4A1","CD79A","CD3E","CD4","CD8A","LYZ","CD14","CST3","percent.diss","percent.mito")
filepath = paste("./", sub.name, "_Umap_Dotplot_QC_checking.jpg", sep='')
ggsave(filepath,DotPlot(sub.object, features = markers), dpi = 300,width = 20,height = 10)
filepath = paste("./", sub.name, "_Umap_featureplot_QC_checking1.jpg", sep='')
ggsave(filepath,FeaturePlot(sub.object, reduction = "umap", features = markers, ncol = 6), dpi = 300,width = 20,height = 10)
filepath = paste("./", sub.name, "_Umap_featureplot_QC_checking2.jpg", sep='')
ggsave(filepath,FeaturePlot(sub.object, reduction = "tsne", features = markers, ncol = 6), dpi = 300,width = 20,height = 10)



## comphrehensive chekcing
table(sub.object$seurat_clusters);table(sub.object$sample.patients); table(sub.object$orig.ident)
markers <- c("EPCAM","KRT19","KRT18","COL1A1","CLDN5","CDH5","PTPRC","MS4A1","CD79A","CD3D","CD8A","LYZ")
MAR_F <- c("COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","COL5A2","COL6A2","COL8A1","COL14A1","C1R","VIM","DCN", # common markers，1-12
           "ACTA2","MYL9","HOPX", # myofibroblast (myCAF), rCAF, 13-15
           "FAP","PDPN","PDGFRA","CXCL12","CFD","DPT","CD74","SLPI") # inflammatory fibroblast (iCAF), pCAF, 最后两个可能是apCAF, 16-21, 22-23
p1 <- FeaturePlot(sub.object, reduction = "umap", features = MAR_F[1:12], ncol = 6)
p2 <- FeaturePlot(sub.object, reduction = "umap", features = MAR_F[13:15], ncol = 3)
p3 <- FeaturePlot(sub.object, reduction = "umap", features = MAR_F[16:23], ncol = 4)
ggsave(paste0("./", sub.name, "_Umap_featureplot_QC_checking3.jpg"),p1/p2/p3, dpi = 300,width = 20,height = 15)

# cluster 16 checking, 180 cells
temp.object <- subset(sub.object, idents = 16)
dim(temp.object); table(temp.object$sample.patients); table(temp.object$orig.ident)
p1 <- FeaturePlot(temp.object, reduction = "umap", features = markers, ncol = 6)
p2 <- DotPlot(temp.object, features = markers) + coord_flip()
p3 <- DotPlot(temp.object, features = MAR_F[1:12]) + coord_flip()
p4 <- DotPlot(temp.object, features = MAR_F[13:15]) + coord_flip()
p5 <- DotPlot(temp.object, features = MAR_F[15:23]) + coord_flip()
p6 <- DimPlot(temp.object, reduction = "umap", group.by = 'sample.patients')
p7 <- DimPlot(temp.object, reduction = "umap", group.by = 'orig.ident')
p8 <- FeaturePlot(temp.object, reduction = "umap", features = "total.srublets.scores")
p9 <- FeaturePlot(temp.object, reduction = "umap", features = "sub.srublets.scores")
p10 <- DimPlot(temp.object, reduction = "umap", group.by = 'Total_DoubleDeconEffect')
p11 <- DimPlot(temp.object, reduction = "umap", group.by = 'sub_DoubleDeconEffect')
p12 <- ggplot(temp.object@meta.data, aes(total.srublets.scores)) + geom_histogram(bins = 50)
p13 <- ggplot(temp.object@meta.data, aes(sub.srublets.scores)) + geom_histogram(bins = 50)
p14 <- ggplot(temp.object@meta.data, aes(percent.diss)) + geom_histogram(bins = 50)
p15 <- ggplot(temp.object@meta.data, aes(percent.mito)) + geom_histogram(bins = 50)
p <- p1 / (p2|p3|p4|p5) / (p6|p7) / (p8|p9|p10|p11) / (p12|p13|p14|p15)
ggsave(paste0("./", sub.name, "_Umap_featureplot_QC_checkingC16.jpg"), p, width = 20, height = 33, dpi = 300, limitsize = F)
high_endo_cells <- colnames(subset(temp.object, subset = CDH5 > 0)); length(high_endo_cells)  # 130/180 cells
high_Epi_cells <- colnames(subset(temp.object, subset = KRT19 > 0)); length(high_Epi_cells)  # 18/180 cells

# cluster 17 checking, 140 cells
temp.object <- subset(sub.object, idents = 17)
dim(temp.object); table(temp.object$sample.patients); table(temp.object$orig.ident)
p1 <- FeaturePlot(temp.object, reduction = "umap", features = markers, ncol = 6)
p2 <- DotPlot(temp.object, features = markers) + coord_flip()
p3 <- DotPlot(temp.object, features = MAR_F[1:12]) + coord_flip()
p4 <- DotPlot(temp.object, features = MAR_F[13:15]) + coord_flip()
p5 <- DotPlot(temp.object, features = MAR_F[15:23]) + coord_flip()
p6 <- DimPlot(temp.object, reduction = "umap", group.by = 'sample.patients')
p7 <- DimPlot(temp.object, reduction = "umap", group.by = 'orig.ident')
p8 <- FeaturePlot(temp.object, reduction = "umap", features = "total.srublets.scores")
p9 <- FeaturePlot(temp.object, reduction = "umap", features = "sub.srublets.scores")
p10 <- DimPlot(temp.object, reduction = "umap", group.by = 'Total_DoubleDeconEffect')
p11 <- DimPlot(temp.object, reduction = "umap", group.by = 'sub_DoubleDeconEffect')
p12 <- ggplot(temp.object@meta.data, aes(total.srublets.scores)) + geom_histogram(bins = 50)
p13 <- ggplot(temp.object@meta.data, aes(sub.srublets.scores)) + geom_histogram(bins = 50)
p14 <- ggplot(temp.object@meta.data, aes(percent.diss)) + geom_histogram(bins = 50)
p15 <- ggplot(temp.object@meta.data, aes(percent.mito)) + geom_histogram(bins = 50)
p <- p1 / (p2|p3|p4|p5) / (p6|p7) / (p8|p9|p10|p11) / (p12|p13|p14|p15)
ggsave(paste0("./", sub.name, "_Umap_featureplot_QC_checkingC17.jpg"), p, width = 20, height = 33, dpi = 300, limitsize = F)
high_endo_cells <- colnames(subset(temp.object, subset = CDH5 > 0)); length(high_endo_cells)  # 4/140 cells
high_Epi_cells <- colnames(subset(temp.object, subset = KRT19 > 0)); length(high_Epi_cells)  # 77/140 cells

# cluster 11/12/14 checking, 902 cells from T5-2
temp.object <- subset(sub.object, idents = c(11,12,14))
dim(temp.object); table(temp.object$sample.patients); table(temp.object$orig.ident)
p1 <- FeaturePlot(temp.object, reduction = "umap", features = markers, ncol = 6)
p2 <- DotPlot(temp.object, features = markers) + coord_flip()
p3 <- DotPlot(temp.object, features = MAR_F[1:12]) + coord_flip()
p4 <- DotPlot(temp.object, features = MAR_F[13:15]) + coord_flip()
p5 <- DotPlot(temp.object, features = MAR_F[15:23]) + coord_flip()
p6 <- DimPlot(temp.object, reduction = "umap", group.by = 'sample.patients')
p7 <- DimPlot(temp.object, reduction = "umap", group.by = 'orig.ident')
p8 <- FeaturePlot(temp.object, reduction = "umap", features = "total.srublets.scores")
p9 <- FeaturePlot(temp.object, reduction = "umap", features = "sub.srublets.scores")
p10 <- DimPlot(temp.object, reduction = "umap", group.by = 'Total_DoubleDeconEffect')
p11 <- DimPlot(temp.object, reduction = "umap", group.by = 'sub_DoubleDeconEffect')
p12 <- ggplot(temp.object@meta.data, aes(total.srublets.scores)) + geom_histogram(bins = 50)
p13 <- ggplot(temp.object@meta.data, aes(sub.srublets.scores)) + geom_histogram(bins = 50)
p14 <- ggplot(temp.object@meta.data, aes(percent.diss)) + geom_histogram(bins = 50)
p15 <- ggplot(temp.object@meta.data, aes(percent.mito)) + geom_histogram(bins = 50)
p <- p1 / (p2|p3|p4|p5) / (p6|p7) / (p8|p9|p10|p11) / (p12|p13|p14|p15)
ggsave(paste0("./", sub.name, "_Umap_featureplot_QC_checkingC11C12C14.jpg"), p, width = 20, height = 33, dpi = 300, limitsize = F)
# ==》细胞基本都受到了LYZ基因的影响，删除细胞

saveRDS(sub.object, file = paste0("./", sub.name, ".RDS"))







#############################################
########## detailed QC Fib cells1 ###########
#############################################
cells_unqualified1 <- colnames(subset(sub.object, idents = c(10,16,17))) # doublets
cells_unqualified2 <- colnames(subset(sub.object, idents = c(18))) # small clusters far away
cells_unqualified4 <- colnames(subset(sub.object, idents = c(11,12,14))) # sample T5-2

# adding infomation 
KeepInfo <- rep("Qualified",dim(sub.object)[2])
KeepInfo[colnames(sub.object) %in% c(cells_unqualified1,cells_unqualified2,cells_unqualified3,cells_unqualified4)] <- "Unqualified"
sub.object$KeepInfo <- factor(KeepInfo); table(sub.object$KeepInfo)
p1 <-  DimPlot(sub.object, reduction = "umap", group.by = 'seurat_clusters', label = T)
p2 <-  DimPlot(sub.object, reduction = "umap", group.by = 'KeepInfo', cols = c("gray","red"))
ggsave(paste0("./", sub.name, "_Umap_QualityMarked.jpg"),p1/p2, dpi = 300,width = 6,height = 10)
saveRDS(sub.object, file = paste0("./", sub.name, ".RDS"))  

## subset Qualified cells 
subQC.object <- subset(sub.object, subset = KeepInfo == "Qualified"); dim(subQC.object)
subQC.name <- "AllSample_Fib_QC"
subQC.object$KeepInfo <- NULL
saveRDS(subQC.object, file = paste0("./", subQC.name, ".RDS"))  
rm(sub.object);rm(sub.name)

## basic steps
# exclusing some genes expressed less than min cells
dim(subQC.object)
subQC.object <- ObjectUpdate(subQC.object, min.cells = 10);dim(subQC.object)
# seurat pipeline
subQC.object <- NormalizeData(subQC.object, normalization.method = "LogNormalize", scale.factor = 10000)
subQC.object <- FindVariableFeatures(subQC.object, selection.method = "vst", nfeatures = max(round(dim(subQC.object)[1] * 0.1), 2000))
subQC.object <- ScaleData(subQC.object, features = rownames(subQC.object))
subQC.object <- RunPCA(subQC.object, features = VariableFeatures(object = subQC.object))
subQC.object <- RunUMAP(subQC.object, dims = 1:20)
subQC.object <- RunTSNE(subQC.object, dims = 1:20)
subQC.object <- FindNeighbors(subQC.object, dims = 1:20)
subQC.object <- FindClusters(subQC.object, resolution = 0.5)

p1 <- DimPlot(subQC.object, reduction = "umap", label = TRUE)
p2 <- DimPlot(subQC.object, reduction = "tsne", label = TRUE)
ggsave(paste0("./", subQC.name, "_Reduction_raw.jpg"), p1|p2, dpi = 300, width = 12,height = 5)

p1 <- DimPlot(subQC.object, reduction = "umap", group.by = 'sample.locations')
p2 <- DimPlot(subQC.object, reduction = "umap", group.by = 'sample.cancers')
p3 <- DimPlot(subQC.object, reduction = "umap", group.by = 'sample.results')
p4 <- DimPlot(subQC.object, reduction = "umap", group.by = 'sample.times')
p5 <- DimPlot(subQC.object, reduction = "umap", group.by = 'sample.patients')
p6 <- DimPlot(subQC.object, reduction = "umap", group.by = 'orig.ident')
p <- (p1|p2|p3) / (p4|p5|p6)
filepath = paste0("./", subQC.name, "_Umap_samples_summary.jpg")
ggsave(filepath, p, dpi = 300, width = 17, height = 10)

p1 <- DimPlot(subQC.object, reduction = "tsne", group.by = 'sample.locations')
p2 <- DimPlot(subQC.object, reduction = "tsne", group.by = 'sample.cancers')
p3 <- DimPlot(subQC.object, reduction = "tsne", group.by = 'sample.results')
p4 <- DimPlot(subQC.object, reduction = "tsne", group.by = 'sample.times')
p5 <- DimPlot(subQC.object, reduction = "tsne", group.by = 'sample.patients')
p6 <- DimPlot(subQC.object, reduction = "tsne", group.by = 'orig.ident')
p <- (p1|p2|p3) / (p4|p5|p6)
filepath = paste0("./", subQC.name, "_Tsne_samples_summary.jpg")
ggsave(filepath, p, dpi = 300, width = 17, height = 10)

saveRDS(subQC.object, file = paste0("./", subQC.name, ".RDS"))







