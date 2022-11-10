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

setwd("./4_TME/")





#############################################
# bulk44
#############################################
## bulk44
bulkinfo <- read.table("../0_data/bulkRNA/bulkRNA_Survival44_SampleInfo.csv", header = T, sep = ",", stringsAsFactors = F)

ann <- bulkinfo; rownames(ann) <- bulkinfo$SampleID; ann <- ann[,c(2,7)]; table(ann$Info); table(ann$Status_up)
ann4 <- ann; colnames(ann4) <- c("Group","Status"); head(ann4)

test_name <- "sub15"
TypeOrder <- c("EndoBlood2","CAF2","CD8T_exhaust","CD4Treg","cDC2","macrophageM2",
               "EndoBlood1","CAF1","CD8T_naive","CD8T_cytoxic","CD4Th","monocyte","pDC","moDC","cDC1"); rn <- 6
FileType <- "Sbatch"


fraction.re <- read.table(paste0("bulk44_",test_name, "_", FileType,"_Results.csv"), header = T, stringsAsFactors = F, sep = ",")
rownames(fraction.re) <- fraction.re[,1]; fraction.re <- fraction.re[, c(-1,-((ncol(fraction.re)-2):ncol(fraction.re)))]
temp.info1 <- ann4[ann4$Group %in% "NaiveGood",]; temp.info1 <- temp.info1[order(temp.info1$Status),]
temp.info2 <- ann4[ann4$Group %in% "NaiveBad",]; temp.info2 <- temp.info2[order(temp.info2$Status),]
temp.info3 <- ann4[ann4$Group %in% "TreatGood",]; temp.info3 <- temp.info3[order(temp.info3$Status),]
temp.info4 <- ann4[ann4$Group %in% "TreatBad",]; temp.info4 <- temp.info4[order(temp.info4$Status),]
SampleOrder <- rownames(rbind(temp.info1, temp.info2, temp.info3, temp.info4))

ann_colors = list(Status = c(Alive= "#8B008B", Dead = "#FFE4B5"),
                  Group = c(NaiveGood= "#00468B", NaiveBad = "#ED0000", TreatGood = "#0099B4", TreatBad = "#42B540"))
pheatmap(t(fraction.re[SampleOrder, TypeOrder]), color = colorRampPalette(c("blue", "white", "red"))(100), scale = "column",annotation_colors = ann_colors,
         cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, annotation_col = ann4, gaps_col = c(15,30,37), gaps_row = c(rn), 
         fontsize_row = 16, fontsize_col = 6, filename = paste0("bulk44_", test_name, "_", FileType, "_样本物理取样分组_update_plotting.jpg"), width = 12, height = 8)


### validation 12 samples
fraction.re <- read.table("./bulkRNA_FinalTest12_sub15_Sbatch.csv", header = T, sep = ",")
ann <- data.frame(PredRe = c("Correct", "Correct", "Correct", "Correct", "Correct", "Incorrect", "Correct", "Correct", "Correct", "Correct", "Incorrect", "Correct"),
                  Group = c(rep("NaiveGood", 7), rep("NaiveBad", 5)))
rownames(ann) <- fraction.re$Mixture
ann_colors = list(Group = c("NaiveGood"= "#00468B", "NaiveBad" = "#ED0000"),
                  PredRe = c("Correct" = "#ff0066", "Incorrect"="black"))
TypeOrder <- c("EndoBlood2","CAF2","CD8T_exhaust","CD4Treg","cDC2","macrophageM2",
               "EndoBlood1","CAF1","CD8T_naive","CD8T_cytoxic","CD4Th","monocyte","pDC","moDC","cDC1"); rn <- 6
fraction.re <- fraction.re[,-1]; rownames(fraction.re) <- rownames(ann)
pheatmap(t(fraction.re[, TypeOrder]), color = colorRampPalette(c("blue", "white", "red"))(100), scale = "column", annotation_colors = ann_colors,
         cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, annotation_col = ann, gaps_col = c(7), gaps_row = c(rn),
         fontsize_row = 16, fontsize_col = 6, filename = "bulkRNA_FinalTest12_sub15_batch_update.jpg", width = 14, height = 8)


## NPC88 + 15sub
bulkinfo <- read.table("../0_data/bulkRNA/NPC88_SampleInfo.csv", header = T, sep = ",", stringsAsFactors = F)
colnames(bulkinfo)[1] <- "SampleID"; bulkinfo$Status_up <- bulkinfo$Status; bulkinfo$Status_up[bulkinfo$Status_up == 0] <- "DiseaseFree"; bulkinfo$Status_up[bulkinfo$Status_up == 1] <- "DiseaseProgress"
ann <- bulkinfo; rownames(ann) <- bulkinfo$SampleID; ann <- ann[,c(4,6)]; table(ann$ClinicalStage); table(ann$Status_up)
ann$Status_up <- factor(ann$Status_up, levels = c("DiseaseFree","DiseaseProgress")); table(ann$Status_up)

test_name <- "sub15"
TypeOrder <- c("EndoBlood2","CAF2","CD8T_exhaust","CD4Treg","cDC2","macrophageM2",
               "EndoBlood1","CAF1","CD8T_naive","CD8T_cytoxic","CD4Th","monocyte","pDC","moDC","cDC1"); rn <- 6
FileType <- "Sbatch"


fraction.re <- read.table(paste0("NPC88_",test_name, "_", FileType,"_Results.csv"), header = T, stringsAsFactors = F, sep = ",")
rownames(fraction.re) <- fraction.re[,1]; fraction.re <- fraction.re[, c(-1,-((ncol(fraction.re)-2):ncol(fraction.re)))]

p <- pheatmap(t(fraction.re[,TypeOrder]), color = colorRampPalette(c("blue", "white", "red"))(100), scale = "column",
              cluster_rows = F, cluster_cols = T, show_rownames = T, show_colnames = T, annotation_col = ann, gaps_row = c(rn),
              fontsize_row = 15, fontsize_col = 6, cutree_cols = 3)
col_groups <- cutree(p$tree_col,k=3); table(col_groups)

col_groups <- data.frame(SampleID = names(col_groups), col_groups, row.names = NULL)
col_groups$group <- col_groups$col_groups; col_groups$group[col_groups$col_groups == 1] <- "Group3";col_groups$group[col_groups$col_groups == 2] <- "Group1";col_groups$group[col_groups$col_groups == 3] <- "Group2";
info3 <- merge(col_groups[,c(1,3)], bulkinfo[,c(1,6)], by = "SampleID"); colnames(info3) <- c("SampleID", "group", "Status"); table(info3$group)
temp.info1 <- info3[info3$group == "Group1",]; temp.info1 <- temp.info1[order(temp.info1$Status, decreasing = F),]
temp.info2 <- info3[info3$group == "Group2",]; temp.info2 <- temp.info2[order(temp.info2$Status, decreasing = F),]
temp.info3 <- info3[info3$group == "Group3",]; temp.info3 <- temp.info3[order(temp.info3$Status, decreasing = T),]
info3 <- rbind(temp.info1, temp.info2, temp.info3)
ann <- info3; rownames(ann) <- ann$SampleID; ann <- ann[,c(2,3)]; colnames(ann) <- c("Group", "Event")

ann_colors <- list(Event = c(DiseaseFree= "#8B008B", DiseaseProgress = "#FFE4B5"), Group = c(Group1= "#374E55", Group2 = "#DF8F44", Group3 = "#00A1D5"))
pheatmap(t(fraction.re[, TypeOrder]), color = colorRampPalette(c("blue", "white", "red"))(100), scale = "column",annotation_colors = ann_colors,
         cluster_rows = F, cluster_cols = T, show_rownames = T, show_colnames = F, annotation_col = ann, gaps_row = c(rn), cutree_cols = 3, 
         fontsize_row = 15, fontsize_col = 6, filename = paste0("NPC88_", test_name, "_", FileType, "_样本均数学聚类_plotting.jpg"), width = 12, height = 8)






