suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsci)
  library(scales)
  library(ggthemr)
  library(patchwork)
  library(Matrix)
  library(stats)
  
  library(Seurat)
  library(GSVA)
  library(pheatmap)
})



#############################################
################# preparation ###############
#############################################
set.seed(1)
wkdir = "/home/mztao/scRNA_analysis/data_analysis/HPC/results/"; setwd(wkdir)
source("/home/mztao/scRNA_analysis/data_analysis/HPC/scripts/common_functions.r")

current.dir <- paste0(getwd(), "/", "MergeTotal_whole/")
setwd(current.dir)


###########################################
# new object updating and subtypes assign
###########################################
Allsample <- readRDS("~/scRNA_analysis/data_analysis/HPC/results/0data/AllSample.RDS"); dim(Allsample)

AllSample_Epi <- readRDS("../0data/AllSample_Epi.RDS")
AllSample_Endo <- readRDS("../0_data/AllSample_Endo_final.RDS")
AllSample_Fib <- readRDS("../0_data/AllSample_Fib_final.RDS")
AllSample_Tcell <- readRDS("../0_data/AllSample_Tcell.RDS")
AllSample_Bcell <- readRDS("../0_data/AllSample_Bcell.RDS")
AllSample_Mcell <- readRDS("../0_data/AllSample_Mcell.RDS")

### assign the cell types, six main cell types
detailed1 <- rep("Unqualified",dim(Allsample)[2])
detailed1[colnames(Allsample) %in% colnames(AllSample_Epi)] <- "Epithelial"  
detailed1[colnames(Allsample) %in% colnames(AllSample_Endo)] <- "Endothelial"
detailed1[colnames(Allsample) %in% colnames(AllSample_Fib)] <- "Fibroblast"
detailed1[colnames(Allsample) %in% colnames(AllSample_Tcell)] <- "Tcell"  
detailed1[colnames(Allsample) %in% colnames(AllSample_Bcell)] <- "Bcell"
detailed1[colnames(Allsample) %in% colnames(AllSample_Mcell)] <- "Myeloid"
Allsample$detailed1 <- factor(detailed1); table(Allsample$detailed1)

Idents(Allsample) <- Allsample$detailed1
Allsample_update <- subset(Allsample, ident = c("Unqualified"), invert = T); table(Allsample_update$detailed1)
Allsample_update$detailed1 <- factor(as.character(Allsample_update$detailed1), levels = c("Epithelial","Endothelial","Fibroblast","Tcell","Bcell","Myeloid")); table(Allsample_update$detailed1)
DimPlot(Allsample_update, reduction = "tsne", group.by = "detailed1", label = F)


### 样本分组标记update
table(Allsample_update$orig.ident); length(table(Allsample_update$orig.ident))
sample.samples <- plyr::revalue(Allsample_update$orig.ident, c("T0_CD45" = "T0_1", "T0_NonCD45" = "T0_1", "NT4" = "NT4", "NT5" = "NT5", "T1_1" = "T1_1", "T1_2" = "T1_2", "T2_1" = "T2_1", "T3_1" = "T3_1", "T3_2" = "T3_2", "T4_1" = "T4_1", "T4_2" = "T4_2", "T5_1" = "T5_1", "T5_2" = "T5_2", "T6_1" = "T6_1", "T6_2" = "T6_2", "T7" = "T7_2")); table(sample.samples); length(sample.samples)
sample.patients <- plyr::revalue(Allsample_update$orig.ident, c("T0_CD45" = "P0", "T0_NonCD45" = "P0", "NT4" = "P4", "NT5" = "P5", "T1_1" = "P1", "T1_2" = "P1", "T2_1" = "P2", "T3_1" = "P3", "T3_2" = "P3", "T4_1" = "P4", "T4_2" = "P4", "T5_1" = "P5", "T5_2" = "P5", "T6_1" = "P6", "T6_2" = "P6", "T7" = "P7")); table(sample.patients); length(sample.patients)
sample.locations <- plyr::revalue(Allsample_update$orig.ident, c("T0_CD45" = "Tumor", "T0_NonCD45" = "Tumor", "NT4" = "NT", "NT5" = "NT", "T1_1" = "Tumor", "T1_2" = "Lym", "T2_1" = "Tumor", "T3_1" = "Tumor", "T3_2" = "Tumor", "T4_1" = "Tumor", "T4_2" = "Tumor", "T5_1" = "Tumor", "T5_2" = "Tumor", "T6_1" = "Tumor", "T6_2" = "Tumor", "T7" = "Tumor")); table(sample.locations); length(sample.locations)
sample.times <- plyr::revalue(Allsample_update$orig.ident, c("T0_CD45" = "naive", "T0_NonCD45" = "naive", "NT4" = "naive", "NT5" = "naive", "T1_1" = "naive", "T1_2" = "treat", "T2_1" = "naive", "T3_1" = "naive", "T3_2" = "treat", "T4_1" = "naive", "T4_2" = "treat", "T5_1" = "naive", "T5_2" = "treat", "T6_1" = "naive", "T6_2" = "treat", "T7" = "treat")); table(sample.times); length(sample.times)
sample.results <- plyr::revalue(Allsample_update$orig.ident, c("T0_CD45" = "Good", "T0_NonCD45" = "Good", "NT4" = "Bad", "NT5" = "Bad", "T1_1" = "Good", "T1_2" = "Good", "T2_1" = "Bad", "T3_1" = "Good", "T3_2" = "Good", "T4_1" = "Bad", "T4_2" = "Bad", "T5_1" = "Bad", "T5_2" = "Bad", "T6_1" = "Bad", "T6_2" = "Bad", "T7" = "Bad")); table(sample.results); length(sample.results)
sample.fourgroups <- plyr::revalue(Allsample_update$orig.ident, c("T0_CD45" = "NaiveGood", "T0_NonCD45" = "NaiveGood", "NT4" = "NaiveBad", "NT5" = "NaiveBad", "T1_1" = "NaiveGood", "T1_2" = "TreatGood", "T2_1" = "NaiveBad", "T3_1" = "NaiveGood", "T3_2" = "TreatGood", "T4_1" = "NaiveBad", "T4_2" = "TreatBad", "T5_1" = "NaiveBad", "T5_2" = "TreatBad", "T6_1" = "NaiveBad", "T6_2" = "TreatBad", "T7" = "TreatBad")); table(sample.fourgroups); length(sample.fourgroups) 
sample.sixgroups <- plyr::revalue(Allsample_update$orig.ident, c("T0_CD45" = "NaiveGood", "T0_NonCD45" = "NaiveGood", "NT4" = "NT", "NT5" = "NT", "T1_1" = "NaiveGood", "T1_2" = "Lym", "T2_1" = "NaiveBad", "T3_1" = "NaiveGood", "T3_2" = "TreatGood", "T4_1" = "NaiveBad", "T4_2" = "TreatBad", "T5_1" = "NaiveBad", "T5_2" = "TreatBad", "T6_1" = "NaiveBad", "T6_2" = "TreatBad", "T7" = "TreatBad")); table(sample.sixgroups); length(sample.sixgroups)

Allsample_update$sample.samples <- factor(sample.samples, levels = unique(sample.samples)); table(Allsample_update$sample.samples)
Allsample_update$sample.patients <- factor(sample.patients, levels = unique(sample.patients)); table(Allsample_update$sample.patients)
Allsample_update$sample.locations <- factor(sample.locations, levels = c("Lym","NT","Tumor")); table(Allsample_update$sample.locations)
Allsample_update$sample.times <- factor(sample.times, levels = c("naive","treat")); table(Allsample_update$sample.times)
Allsample_update$sample.results <- factor(sample.results, levels = c("Good","Bad")); table(Allsample_update$sample.results)
Allsample_update$sample.fourgroups <- factor(sample.fourgroups, levels = c("NaiveGood","NaiveBad","TreatGood","TreatBad")); table(Allsample_update$sample.fourgroups)
Allsample_update$sample.sixgroups <- factor(sample.sixgroups, levels = c("Lym","NT","NaiveGood","NaiveBad","TreatGood","TreatBad")); table(Allsample_update$sample.sixgroups)

### assign the detailed subtypes
table(AllSample_Epi$CancerInfo); DimPlot(AllSample_Epi, reduction = "tsne", group.by = "CancerInfo")
table(AllSample_Endo$detailed3); DimPlot(AllSample_Endo, reduction = "tsne", group.by = "detailed3")
table(AllSample_Fib$detailed3); DimPlot(AllSample_Fib, reduction = "tsne", group.by = "detailed3")
table(AllSample_Tcell$detailed3); DimPlot(AllSample_Tcell, reduction = "tsne", group.by = "detailed3")
table(AllSample_Bcell$detailed3); DimPlot(AllSample_Bcell, reduction = "tsne", group.by = "detailed3")
table(AllSample_Mcell$detailed3); DimPlot(AllSample_Mcell, reduction = "tsne", group.by = "detailed3")

Epi_subs <- as.character(AllSample_Epi$CancerInfo); names(Epi_subs) <- colnames(AllSample_Epi)
Endo_subs <- as.character(AllSample_Endo$detailed3); names(Endo_subs) <- colnames(AllSample_Endo)
Fib_subs <- as.character(AllSample_Fib$detailed3); names(Fib_subs) <- colnames(AllSample_Fib)
T_subs <- as.character(AllSample_Tcell$detailed3); names(T_subs) <- colnames(AllSample_Tcell)
B_subs <- as.character(AllSample_Bcell$detailed3); names(B_subs) <- colnames(AllSample_Bcell)
M_subs <- as.character(AllSample_Mcell$detailed3); names(M_subs) <- colnames(AllSample_Mcell)
all <- c(Epi_subs, Endo_subs, Fib_subs, B_subs, M_subs, T_subs)
Allsample_update$detailed2 <- factor(all); table(Allsample_update$detailed2)

detailed2 <- as.character(Allsample_update$detailed2)
detailed2[colnames(Allsample_update) %in% colnames(subset(Allsample_update, subset = detailed2 == "Normal_cell"))] <- "NormalEpi"
detailed2[colnames(Allsample_update) %in% colnames(subset(Allsample_update, subset = detailed2 == "Maligant_Cell"))] <- "MaligantEpi"
Allsample_update$detailed2 <- factor(detailed2)
Allsample_update$detailed2 <- factor(Allsample_update$detailed2, levels = c("NormalEpi","MaligantEpi","EndoLym","EndoBlood1","EndoBlood2",
                                                                            "proFib","MyoFib","CAF1","CAF2","CAF3",
                                                                            "CD8T_naive","CD8T_memory","CD8T_cytoxic1","CD8T_cytoxic2","CD8T_cytoxic3","CD8T_exhaust","CD8T_Ebo","CD4Treg_naive","CD4Treg_act","CD4Treg_Ebo","CD4Th_naive","CD4Th_act","CD4Th_exhaust","NK",
                                                                            "B_GC","B_MemoInter","B_Memory","B_PlasamIgA","B_PlasamIgG",
                                                                            "neutrophil","mast","monocyte","pDC","moDC","cDC1","cDC2","macrophageM1","macrophageM2"))
table(Allsample_update$detailed2)
saveRDS(Allsample_update, file = "../0_data/Allsample_final.RDS")


