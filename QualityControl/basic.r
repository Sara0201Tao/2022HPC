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
annotations <- getAnnotations("/data/mztao/ref/scRNA/GRCh38/refdata-gex-GRCh38-2020-A/genes/genes.gtf", is.mouse =FALSE)

data_root_dir <- "/home/mztao/scRNA_analysis/data_analysis/HPC/results/raw_counts/"
sample.dir <- c("HPC_T0_CD45/", "HPC_T0_NonCD45/", "HPC_T1_1/", "HPC_T1_2/", "HPC_T2_1/", "HPC_T3_1/", "HPC_T3_2/", "HPC_T4_1/", "HPC_T4_2/", "HPC_NT4/", "HPC_T5_1/", "HPC_T5_2/", "HPC_NT5/", "HPC_T6_1/", "HPC_T6_2/", "HPC_T7/")
sample.object <- c("T0_CD45", "T0_NonCD45", "T1_1", "T1_2", "T2_1", "T3_1", "T3_2", "T4_1", "T4_2", "NT4", "T5_1", "T5_2", "NT5", "T6_1", "T6_2", "T7")
sample.name <- c("T0_CD45", "T0_NonCD45", "T1_1", "T1_2", "T2_1", "T3_1", "T3_2", "T4_1", "T4_2", "NT4", "T5_1", "T5_2", "NT5", "T6_1", "T6_2", "T7")



###########################################
# create object and do basic QC by seurat
###########################################
dir.create("./step1_MergeMainTypes/BasicQC", recursive = T)
setwd("./step1_MergeMainTypes/BasicQC")

RawObjects <- vector(mode = "list", length = length(sample.name));
names(RawObjects) <- sample.name
for (i in (1:length(sample.dir))){
  # load the dataset 
  i <- i + 1
  current.dir <- sample.dir[i]
  current.object <- sample.object[i]
  current.name <- sample.name[i]
  
  data.dir <- paste(data_root_dir,current.dir,"outs/filtered_feature_bc_matrix/",sep="")
  load.data <- Read10X(data.dir = data.dir)
  dim(load.data) 
  
  current.object <- scStat(load.data, current.name, annotations, mini_cells = 5);dim(current.object)
  current.object <- subset(current.object, subset = nFeature_RNA >= 200 & nFeature_RNA <= 7500 & percent.mito <= 25)
  dim(current.object)
  
  # basic seurat steps
  current.object <- NormalizeData(current.object, normalization.method = "LogNormalize", scale.factor = 10000)
  current.object <- FindVariableFeatures(current.object, selection.method = "vst", nfeatures = 2000)
  current.object <- ScaleData(current.object, features = rownames(current.object))
  current.object <- RunPCA(current.object, features = VariableFeatures(object = current.object))
  current.object <- FindNeighbors(current.object, dims = 1:15)
  current.object <- FindClusters(current.object, resolution = 0.4)
  current.object <- RunUMAP(current.object, dims = 1:15)
  
  # visualization
  markers <- c("EPCAM","COL1A1","CLDN5","PTPRC","CD19","CD79A","CD3E","CD4","CD8A","LYZ","CD14","CST3")
  p1 <- DimPlot(current.object, reduction = "umap", label = TRUE)
  p2 <- FeaturePlot(current.object, reduction = "umap", features = markers, ncol = 3)
  
  filepath = paste("./", current.name, "_Umap_PrimaryChecking.jpg", sep='')
  ggsave(filepath, p1|p2, dpi = 300, width = 18, height = 8)

  # save the object into list
  RawObjects[[current.name]] <- current.object
}
saveRDS(RawObjects, file = paste0("./RawObjects.RDS"))





###########################################
############# Merging data ################
###########################################
current.object <- merge(x = RawObjects[[1]], y = RawObjects[2:length(sample.name)], add.cell.ids = sample.name, project = "AllSamples")
dim(current.object) # 28038 104094
table(current.object$orig.ident)

## save the object
setwd(wkdir);setwd("./step1_MergeMainTypes/")
current.name <- "AllSamples"
saveRDS(current.object, file = paste0("./", current.name, ".RDS"))
rm(RawObjects)

## basic seurat steps
current.object <- NormalizeData(current.object, normalization.method = "LogNormalize", scale.factor = 10000)
current.object <- FindVariableFeatures(current.object, selection.method = "vst", nfeatures = max(round(dim(current.object)[1] * 0.1), 2000))
current.object <- ScaleData(current.object, features = rownames(current.object))
current.object <- RunPCA(current.object, features = VariableFeatures(object = current.object))
filepath = paste("./", current.name, "_ElbowPlot.jpg", sep='')
ggsave(filepath,ElbowPlot(current.object,ndims = 40), dpi = 300)
# cluster the cells and visualization
current.object <- FindNeighbors(current.object, dims = 1:20)
current.object <- FindClusters(current.object, resolution = 0.5)
current.object <- RunUMAP(current.object, dims = 1:20)
filepath = paste("./", current.name, "_Umap_raw.jpg", sep='')
ggsave(filepath,DimPlot(current.object, reduction = "umap", label = TRUE), dpi = 300)

saveRDS(current.object, file = paste0("./", current.name, ".RDS"))




###########################################
# advanced QC in all sample cells
###########################################
##### basic ones
# 1. plot
filepath = paste("./", current.name, "_BasicStats.jpg", sep='')
ggsave(filepath, FeatureScatter(current.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"), dpi = 300)

# 2. disssociation effect calculating for further use
DissGene.stats <- current.object$percent.diss
OverRatio_DissGene <- length(DissGene.stats[DissGene.stats > 5.75 ]) / length(DissGene.stats)
pdf(file=paste("./", current.name, "_DissliHist1.pdf", sep=''))
hist(DissGene.stats, breaks = 200, xlim = c(0,15))
abline(v = 5.75, col = "red", lwd =2)
text(10,2000,paste0(as.character(round(OverRatio_DissGene*100,2)),"%"))
dev.off()

filepath = paste0("./", current.name, "_Umap_DissEffectionScores.jpg")
ggsave(filepath, FeaturePlot(current.object,reduction = "umap", features = "percent.diss"), dpi = 300,width = 10,height = 8)

subset_diss <- rep("NoEffected",dim(current.object)[2])
subset_diss[current.object$percent.diss > 5.75] <- "Effected"
current.object$diss_info <- subset_diss
filepath = paste0("./", current.name, "_Umap_DissEffection.jpg")
ggsave(filepath, DimPlot(current.object,reduction = "umap", group.by = 'diss_info', cols = c("red","gray")), dpi = 300,width = 10,height = 8)

# 3. proliferation effect checking
ProliGene.stats <- current.object$percent.proli
pdf(file=paste("./", current.name, "_ProliHist.pdf", sep=''))
hist(ProliGene.stats, breaks = 200, xlim = c(0,0.5))
dev.off()

filepath = paste("./", current.name, "_ProliEffection_Umap1.jpg", sep='')
ggsave(filepath,FeaturePlot(current.object, reduction = "umap", features = c("MKI67","NUSAP1","PLK1","CDC20","CDK1")), dpi = 300)
filepath = paste("./", current.name, "_ProliEffection_Umap2.jpg", sep='')
ggsave(filepath,FeaturePlot(current.object, reduction = "umap", features = c("CDKN3","CENPA","BIRC5","PCNA","CCNA2")), dpi = 300)

# 5.1 Srublets
dir.create("./Srublets/")

# split into 5 parts
x  <- runif(dim(current.object)[2], min=0, max=5) #uniform distribution
sub1.object <- current.object[,colnames(current.object)[x<=1]] ; dim(sub1.object)
sub2.object <- current.object[,colnames(current.object)[x>1 & x <= 2]] ; dim(sub2.object)
sub3.object <- current.object[,colnames(current.object)[x>2 & x <= 3]] ; dim(sub3.object)
sub4.object <- current.object[,colnames(current.object)[x>3 & x <= 4]] ; dim(sub4.object)
sub5.object <- current.object[,colnames(current.object)[x>4]] ; dim(sub5.object)
Out10x(sub1.object, "./Srublets/", "AllSample_sub1")
Out10x(sub2.object, "./Srublets/", "AllSample_sub2")
Out10x(sub3.object, "./Srublets/", "AllSample_sub3")
Out10x(sub4.object, "./Srublets/", "AllSample_sub4")
Out10x(sub5.object, "./Srublets/", "AllSample_sub5")

sub1.scores <- as.data.frame(read.table(paste0("./Srublets/","AllSample_sub1","_doublet_scores.txt")))
rownames(sub1.scores) <- colnames(sub1.object) ; colnames(sub1.scores) <- "scores"
sub2.scores <- as.data.frame(read.table(paste0("./Srublets/","AllSample_sub2","_doublet_scores.txt")))
rownames(sub2.scores) <- colnames(sub2.object) ; colnames(sub2.scores) <- "scores"
sub3.scores <- as.data.frame(read.table(paste0("./Srublets/","AllSample_sub3","_doublet_scores.txt")))
rownames(sub3.scores) <- colnames(sub3.object) ; colnames(sub3.scores) <- "scores"
sub4.scores <- as.data.frame(read.table(paste0("./Srublets/","AllSample_sub4","_doublet_scores.txt")))
rownames(sub4.scores) <- colnames(sub4.object) ; colnames(sub4.scores) <- "scores"
sub5.scores <- as.data.frame(read.table(paste0("./Srublets/","AllSample_sub5","_doublet_scores.txt")))
rownames(sub5.scores) <- colnames(sub5.object) ; colnames(sub5.scores) <- "scores"

# add scores info into the object
total.srublets.scores <- rbind(sub1.scores,sub2.scores,sub3.scores,sub4.scores,sub5.scores)
current.object$total.srublets.scores <- total.srublets.scores

# plotting for visualization
pdf(file=paste("./", current.name, "_total_srublets_scores.pdf", sep=''))
hist(current.object$total.srublets.scores, breaks = 200, xlim = c(0,1))
dev.off()

filepath = paste("./", current.name, "_Umap_total_srublets_scores.jpg", sep='')
ggsave(filepath,FeaturePlot(current.object, reduction = "umap",features = "total.srublets.scores"), dpi = 300, width = 10, height = 8)

threshold <- 0.2
pdf(file=paste("./", current.name, "_total_srublets_scores_threshold.pdf", sep=''))
hist(current.object$total.srublets.scores, breaks = 100, xlim = c(0,1))
abline(v = threshold, col = "red", lwd =2)
text(0.6,5000,paste0(as.character(round(length(which(current.object$total.srublets.scores > threshold)) / dim(current.object)[2]*100, 2)),"%"))
dev.off()

subset_Sdoublets <- rep("NoEffected",dim(current.object)[2])
subset_Sdoublets[current.object$total.srublets.scores > threshold] <- "Effected"
current.object$Total_HighSrubeltsEffect <- as.factor(subset_Sdoublets)
filepath = paste0("./", current.name, "_Umap_total_srublets_scores_group.jpg")
ggsave(filepath,DimPlot(current.object,reduction = "umap", group.by = 'Total_HighSrubeltsEffect', cols = c("red","gray")), dpi = 300,width = 10,height = 8)

rm(sub1.object);rm(sub2.object);rm(sub3.object);rm(sub4.object);rm(sub5.object)
rm(sub1.scores);rm(sub2.scores);rm(sub3.scores);rm(sub4.scores);rm(sub5.scores)
saveRDS(current.object, file = paste0("./", current.name, ".RDS"))

# 5.2 DoubletDecon
library(DoubletDecon)
dir.create("./DoubletDecon/")

Idents(current.object) <- current.object$seurat_clusters
x  <- runif(dim(current.object)[2], min=0, max=5) #uniform distribution
sub1.object <- current.object[,colnames(current.object)[x<=1]] ; dim(sub1.object)
sub2.object <- current.object[,colnames(current.object)[x>1 & x <= 2]] ; dim(sub2.object)
sub3.object <- current.object[,colnames(current.object)[x>2 & x <= 3]] ; dim(sub3.object)
sub4.object <- current.object[,colnames(current.object)[x>3 & x <= 4]] ; dim(sub4.object)
sub5.object <- current.object[,colnames(current.object)[x>4]] ; dim(sub5.object)

sub.objects <- list(sub1.object,sub2.object,sub3.object,sub4.object,sub5.object)
names(sub.objects) <- c("sub1","sub2","sub3","sub4","sub5")
for (i in 1:5){
  temp.name <- names(sub.objects)[i];print(temp.name)
  newFiles <- Improved_Seurat_Pre_Process(sub.objects[[temp.name]], num_genes=50, write_files=FALSE)
  DoubletDecon_results <- Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, groupsFile=newFiles$newGroupsFile, filename=temp.name,
                                             location="./DoubletDecon/",fullDataFile=NULL,removeCC=FALSE,species="hsa", rhop=0.7, write=FALSE,PMF=TRUE,
                                             useFull=FALSE, heatmap=FALSE,centroids=TRUE,num_doubs=100,only50=FALSE,min_uniq=4,nCores=30)
  saveRDS(DoubletDecon_results, file = paste0("./DoubletDecon/DoubletDeconResults_",temp.name,".RDS"))
}

Total_DoubleDeconEffect <- data.frame()
for (i in 1:5){
  print(i)
  DoubletDecon_results <- readRDS(file = paste0("./DoubletDecon/DoubletDeconResults_sub",i,".RDS"))
  doublets <- DoubletDecon_results[["Final_doublets_groups"]]; singlets <- DoubletDecon_results[["Final_nondoublets_groups"]]
  temp_cell_info <- data.frame(groups = c(rep("doublets", dim(doublets)[1]),rep("singlets", dim(singlets)[1])))
  rownames(temp_cell_info) <- c(rownames(doublets), rownames(singlets))
  Total_DoubleDeconEffect <- rbind(Total_DoubleDeconEffect, temp_cell_info)
} 
table(Total_DoubleDeconEffect)
rownames(Total_DoubleDeconEffect) <- gsub("[.]","-",rownames(Total_DoubleDeconEffect))
current.object$Total_DoubleDeconEffect <- Total_DoubleDeconEffect
filepath = paste("./", current.name, "_Umap_total_DoubleDecon_group.jpg", sep='')
ggsave(filepath,DimPlot(current.object,reduction = "umap", group.by = 'Total_DoubleDeconEffect', cols = c("red","gray")), dpi = 300,width = 10,height = 8)





###########################################
# primary main types identification
###########################################
## labelling by markers
MAR <- c("PTPRC")
MAR_Epithelial <- c("EPCAM","KRT19","KRT18","KRT15","KRT5") 
MAR_Fibroblasts <- c("COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","COL5A2","COL6A2","COL8A1","COL14A1","C1R","VIM","DCN", # common markers，
           "ACTA2","MYL9","HOPX", # myofibroblast (myCAF), rCAF, 13-15
           "FAP","PDPN","PDGFRA","CXCL12","CFD","DPT","CD74","SLPI") # inflammatory fibroblast (iCAF), pCAF,
MAR_Endothelial <- c("CLDN5","FLT1","CDH5","RAMP2","PECAM1")
MAR_B <- c("MS4A1","CD19","CD22","CD79A","CD79B") 
MAR_T <- c("CD2","CD3D", "CD3E", "CD3G","CD4","CD8A","GNLY","NKG7")
MAR_myeloid <- c("LYZ","CD68","CD14","CST3","MS4A7","MS4A6A","FCGR3A")

filepath = paste("./", current.name, "_FeaturePlot1.jpg", sep='')
ggsave(filepath,FeaturePlot(current.object, reduction = "umap", features = MAR), dpi = 300,width = 10,height = 8)
filepath = paste("./", current.name, "_FeaturePlot2.jpg", sep='')
ggsave(filepath,FeaturePlot(current.object, reduction = "umap", features = MAR_epithelial, ncol = 3), dpi = 300,width = 10, height = 6)
filepath = paste("./", current.name, "_FeaturePlot3.jpg", sep='')
ggsave(filepath,FeaturePlot(current.object, reduction = "umap", features = MAR_endothelial, ncol = 3), dpi = 300,width = 10,height = 6)
filepath = paste("./", current.name, "_FeaturePlot4.jpg", sep='')
ggsave(filepath,FeaturePlot(current.object, reduction = "umap", features = MAR_Fibroblasts, ncol = 6), dpi = 300,width = 25, height = 16)
filepath = paste("./", current.name, "_FeaturePlot5.jpg", sep='')
ggsave(filepath,FeaturePlot(current.object, reduction = "umap", features = MAR_B, ncol = 3), dpi = 300,width = 10, height = 6)
filepath = paste("./", current.name, "_FeaturePlot6.jpg", sep='')
ggsave(filepath,FeaturePlot(current.object, reduction = "umap", features = MAR_T, ncol = 4), dpi = 300,width = 13, height = 6)
filepath = paste("./", current.name, "_FeaturePlot7.jpg", sep='')
ggsave(filepath,FeaturePlot(current.object, reduction = "umap", features = MAR_myeloid, ncol = 4), dpi = 300,width = 13, height = 6)


## singleR mapping
library(SingleR)
library(BiocParallel)

SourceFile.dir <- "/home/mztao/scRNA_analysis/sup_info/singleR/"
RefFiles <- list.files(SourceFile.dir, pattern = ".rds");RefFiles
ref_data <- readRDS(file = paste0(SourceFile.dir, RefFiles[3]))
table(ref_data$label.main)

test_data <- GetAssayData(current.object, slot="data")
# 1. main types with single mode
MainType.pred <- SingleR(test = test_data, ref = ref_data, labels = ref_data$label.main,
                         de.method = "classic", method = "single",
                         assay.type.test = "logcounts", assay.type.ref = "logcounts", BPPARAM=MulticoreParam(30))
current.object$MainType.pred_HPCA_single <- MainType.pred$labels
filepath = paste("./", current.name, "_Umap_MainType_HPCA_single.jpg", sep='')
ggsave(filepath, DimPlot(current.object, reduction = "umap", group.by = "MainType.pred_HPCA_single", label = T), dpi = 300,width = 10,height = 8)
# 2. main type with cluster mode
MainType.pred <- SingleR(test = test_data, ref = ref_data, labels = ref_data$label.main,
                         de.method = "classic", method = "cluster", clusters = current.object$seurat_clusters,
                         assay.type.test = "logcounts", assay.type.ref = "logcounts", BPPARAM=MulticoreParam(15))
clusters <- current.object$seurat_clusters
levels(clusters) <- MainType.pred$labels
current.object$MainType.pred_HPCA_cluster <- clusters
filepath = paste("./", current.name, "_Umap_MainType_HPCA_cluster.jpg", sep='')
ggsave(filepath, DimPlot(current.object, reduction = "umap", group.by = "MainType.pred_HPCA_cluster", label = T), dpi = 300,width = 10,height = 8)

saveRDS(current.object, file = paste0("./", current.name, ".RDS"))












###############################################
# detailed check and maintype assignment primary
###############################################
# epi
Epi_cells1 <- colnames(subset(current.object, idents = c(7,8,11,13,18,19,20))) 
Epi_cells2 <- colnames(subset(current.object, idents = c(16), subset = UMAP_2 > -8 & UMAP_1 < 8))
Epi_cells <- c(Epi_cells1, Epi_cells2)
# endo
Endo_cells1 <- colnames(subset(current.object, idents = c(3,23))) 
Endo_cells2 <- colnames(subset(current.object, idents = c(17), subset = UMAP_2 > 0))
Endo_cells3 <- colnames(subset(current.object, idents = c(4), subset = UMAP_2 > 0))
Endo_cells4 <- colnames(subset(current.object, idents = c(14), subset = UMAP_2 > 0))
Endo_cells <- c(Endo_cells1, Endo_cells2, Endo_cells3, Endo_cells4)
# fib
Fib_cells1 <- colnames(subset(current.object, idents = c(6,10))) 
Fib_cells2 <- setdiff(colnames(subset(current.object, idents = c(16))) ,Epi_cells2)
Fib_cells3 <- setdiff(colnames(subset(current.object, idents = c(14))) ,Endo_cells4)
Fib_cells <- c(Fib_cells1, Fib_cells2, Fib_cells3)
# Bcell
B_cells1 <- colnames(subset(current.object, idents = c(2,5,26,24))) 
B_cells2 <- colnames(subset(current.object, idents = c(21), subset = UMAP_1 < 0 & UMAP_2 > 0)) 
B_cells3 <- colnames(subset(temp.object, idents = c(1,3,6))) 
B_cells <- c(B_cells1, B_cells2, B_cells3)
# Tcell
T_cells1 <- colnames(subset(current.object, idents = c(0,1,12))) 
T_cells2 <- colnames(subset(current.object, idents = c(21), subset = UMAP_1 > 0)) 
T_cells3 <- colnames(subset(temp.object, idents = c(0,2,4,5))) 
T_cells <- c(T_cells1, T_cells2, T_cells3)
# Mcell
M_cells1 <- colnames(subset(current.object, idents = c(9,22,25)))
M_cells2 <- colnames(subset(current.object, idents = c(21), subset = UMAP_1 < 0 & UMAP_2 < 0))
M_cells3 <- setdiff(colnames(subset(current.object, idents = c(17))), Endo_cells2)
M_cells4 <- setdiff(colnames(subset(current.object, idents = c(4))), Endo_cells3)
M_cells <- c(M_cells1, M_cells2, M_cells3, M_cells4)
# checking
length(Epi_cells) + length(Endo_cells) + length(Fib_cells) + length(B_cells) + length(T_cells) + length(M_cells) 

# try to assign
maintypes <- rep("NoEffected",dim(current.object)[2])
maintypes[colnames(current.object) %in% Epi_cells] <- "Epithelial"
maintypes[colnames(current.object) %in% Endo_cells] <- "Endothelial"
maintypes[colnames(current.object) %in% Fib_cells] <- "Fibroblast"
maintypes[colnames(current.object) %in% B_cells] <- "B cell"
maintypes[colnames(current.object) %in% T_cells] <- "T cell"
maintypes[colnames(current.object) %in% M_cells] <- "Myeloid"

current.object$maintypes <- factor(maintypes)
filepath = paste("./", current.name, "_Umap_MainTypes_try1.jpg", sep='')
ggsave(filepath,DimPlot(current.object,reduction = "umap", group.by = 'maintypes'), dpi = 300,width = 10,height = 8)

saveRDS(current.object, file = paste0("./", current.name, ".RDS"))
