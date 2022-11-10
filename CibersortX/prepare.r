suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(scales)
  library(ggthemr)
  library(patchwork)
  library(Matrix)
  library(stats)
  library(stringr)
  library(reshape2)
  library(scales)
  
  library(Seurat)
  library(pheatmap)
  library(circlize)
  library(randomcoloR)
  library(RColorBrewer)
  library(UpSetR)
  library(survival)
  library(survminer)
  library(GSVA)
})


#############################################
################# preparation ###############
#############################################
set.seed(1)
wkdir = "/home/mztao/scRNA_analysis/data_analysis/HPC/results/"; setwd(wkdir)
source("/home/mztao/scRNA_analysis/data_analysis/HPC/scripts/common_functions.r")

setwd("./TME/")





#############################################
## prepare seurat object
#############################################
AllSample <- readRDS("../0_data/AllSample_final__Normout.RDS")
table(AllSample$detailed1); table(AllSample$detailed2)

Idents(AllSample) <- AllSample$detailed1
current.object <- subset(AllSample, idents = c("Epithelial","Bcell"), invert = T)
Idents(current.object) <- current.object$detailed2
current.object <- subset(current.object, idents = c("EndoLym","proFib","CAF3","CD8T_memory","CD8T_Ebo","CD4Treg_Ebo","neutrophil"), invert = T)

current.object$detailed1 <- factor(as.character(current.object$detailed1))
current.object$detailed2 <- factor(as.character(current.object$detailed2))
table(current.object$detailed1); table(current.object$detailed2)
saveRDS(current.object, "./temp_AllSample.RDS")



#############################################
# test for input cell subtypes
#############################################
## bulk info 
bulkdata <- read.table("../0_data/bulkRNA/bulkRNA_Survival44_exp.csv", header = T, sep = ",", stringsAsFactors = F)
bulkinfo <- read.table("../0_data/bulkRNA/bulkRNA_Survival44_SampleInfoUp.csv", header = T, sep = ",", stringsAsFactors = F)

ann <- bulkinfo; rownames(ann) <- bulkinfo$SampleID; ann <- ann[,c(2,7)]; table(ann$Info); table(ann$Status_up)
ann2 <- ann[rownames(ann) %in% (bulkinfo[bulkinfo$Info %in% c("NaiveGood","NaiveBad"),"SampleID"]),]
ann4 <- ann

## scRNA data load
AllSample <- readRDS("./temp_AllSample.RDS")
table(AllSample$detailed1); table(AllSample$detailed2)

## 1, all cellp3
Idents(AllSample) <- AllSample$detailed2
test_name <- "all19"
subs <-c("EndoBlood1","EndoBlood2","MyoFib","CAF1","CAF2",
         "CD8T_naive","CD8T_cytoxic","CD8T_exhaust","CD4Treg","CD4Th","NK",
         "mast","monocyte","pDC","moDC","cDC1","cDC2","macrophageM1","macrophageM2")

current.object <- subset(AllSample, idents = subs);dim(current.object)
current.object$detailed1 <- factor(as.character(current.object$detailed1)); print(table(current.object$detailed1))
current.object$detailed2 <- factor(as.character(current.object$detailed2)); print(table(current.object$detailed2))


temp.cells <- union(sample(colnames(current.object), 20000), colnames(subset(current.object, subset = detailed2 == "CD8T_naive")))
temp.object <- current.object[,temp.cells]; dim(temp.object); table(temp.object$detailed2)
count_raw <- Matrix(as.matrix(GetAssayData(temp.object, slot = "counts")), sparse = T); dim(count_raw)
count_raw <- count_raw[rownames(count_raw) %in% rownames(bulkdata),]; dim(count_raw)
count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*1e5)
colnames(count_norm) <- as.character(temp.object@meta.data[,"detailed2"])
singlRef <- cbind(Gene = rownames(count_norm), count_norm)
write.table(singlRef, file = paste0(test_name, "_singleRef.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

# Naive30 testing
temp.info1 <- ann2[ann2$Info %in% "NaiveGood",]; temp.info1 <- temp.info1[order(temp.info1$Status_up),]
temp.info2 <- ann2[ann2$Info %in% "NaiveBad",]; temp.info2 <- temp.info2[order(temp.info2$Status_up, decreasing = T),]
SampleOrder <- rownames(rbind(temp.info1, temp.info2))

TypeOrder <- c("EndoBlood1","CAF2","CD8T_exhaust","CD4Treg","cDC2","macrophageM2",
               "EndoBlood2","MyoFib","CAF1","CD8T_naive","CD8T_cytoxic","CD4Th","NK","mast","monocyte","pDC","moDC","cDC1","macrophageM1"); rn <- 6

FileType <- "Sbatch"

fraction.re <- read.table(paste0("Naive30_",test_name, "_", FileType,"_Results.csv"), header = T, stringsAsFactors = F, sep = ",")
rownames(fraction.re) <- fraction.re[,1]; fraction.re <- fraction.re[, c(-1,-((ncol(fraction.re)-2):ncol(fraction.re)))]
pheatmap(t(fraction.re[SampleOrder, TypeOrder]), color = colorRampPalette(c("blue", "white", "red"))(100), scale = "column",
         cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = T, annotation_col = ann2, gaps_col = c(15), gaps_row = c(rn), 
         fontsize_row = 15, fontsize_col = 6, filename = paste0("Naive30_", test_name, "_", FileType, "_Results.jpg"), width = 10, height = 7)

## 2
Idents(AllSample) <- AllSample$detailed2
test_name <- "sub15"
subs <-c("EndoBlood1","EndoBlood2","CAF1","CAF2",
         "CD8T_naive","CD8T_cytoxic","CD8T_exhaust","CD4Treg","CD4Th",
         "monocyte","pDC","moDC","cDC1","cDC2","macrophageM2")

current.object <- subset(AllSample, idents = subs);dim(current.object)
current.object$detailed1 <- factor(as.character(current.object$detailed1)); print(table(current.object$detailed1))
current.object$detailed2 <- factor(as.character(current.object$detailed2)); print(table(current.object$detailed2))

temp.cells <- union(sample(colnames(current.object), 20000), colnames(subset(current.object, subset = Cellp3 == "CD8T_naive")))
temp.object <- current.object[,temp.cells]; dim(temp.object); table(temp.object$detailed2)
count_raw <- Matrix(as.matrix(GetAssayData(temp.object, slot = "counts")), sparse = T); dim(count_raw)
count_raw <- count_raw[rownames(count_raw) %in% rownames(bulkdata),]; dim(count_raw)
count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*1e5)
colnames(count_norm) <- as.character(temp.object@meta.data[,"detailed2"])
singlRef <- cbind(Gene = rownames(count_norm), count_norm)
write.table(singlRef, file = paste0(test_name, "_singleRef.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

# Naive30 testing
temp.info1 <- ann2[ann2$Info %in% "NaiveGood",]; temp.info1 <- temp.info1[order(temp.info1$Status_up),]
temp.info2 <- ann2[ann2$Info %in% "NaiveBad",]; temp.info2 <- temp.info2[order(temp.info2$Status_up, decreasing = T),]
SampleOrder <- rownames(rbind(temp.info1, temp.info2))

TypeOrder <- c("EndoBlood1","CAF2","CD8T_exhaust","CD4Treg","cDC2","macrophageM2",
               "EndoBlood2","CAF1","CD8T_naive","CD8T_cytoxic","CD4Th","monocyte","pDC","moDC","cDC1"); rn <- 6

FileType <- "Sbatch"

fraction.re <- read.table(paste0("Naive30_",test_name, "_", FileType,"_Results.csv"), header = T, stringsAsFactors = F, sep = ",")
rownames(fraction.re) <- fraction.re[,1]; fraction.re <- fraction.re[, c(-1,-((ncol(fraction.re)-2):ncol(fraction.re)))]
pheatmap(t(fraction.re[SampleOrder, TypeOrder]), color = colorRampPalette(c("blue", "white", "red"))(100), scale = "column",
         cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = T, annotation_col = ann2, gaps_col = c(15), gaps_row = c(rn), 
         fontsize_row = 15, fontsize_col = 6, filename = paste0("Naive30_", test_name, "_", FileType, "_Results.jpg"), width = 10, height = 7)
