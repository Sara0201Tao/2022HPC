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
  library(RColorBrewer)
  library(randomcoloR)
  library(UpSetR)
  library(survival)
  library(survminer)
})


#############################################
################# preparation ###############
#############################################
set.seed(1)
source("/home/mztao/scRNA_analysis/data_analysis/HPC/scripts/common_functions.r")
wkdir = "/home/mztao/scRNA_analysis/data_analysis/HPC/results/"; setwd(wkdir)


group.names <- c("NaiveGood", "NaiveBad", "TreatBad")
alltypes <- c("MaligantEpi","EndoLym","EndoBlood","proFib","MyoFib","CAF","CD8T","CD4Treg","CD4Th","NK","Bcell","neutrophil","mast","monocyte","pDC","DC","macrophage")
alltypes_focused <- alltypes[c(1,3,6,7,8,9,10,14:17)]; alltypes_focused


####################  
# 造血相关
df <- read.table(file = paste0("./cellp2_015", "_EC造血相关",".csv"), header = T, sep = ",")

g1 <- c("EndoBlood"); g2 <- setdiff(alltypes_focused, g1)
CTpairs1 <- c(); CTpairs2 <- c()
for (i in g1){for (j in g2){CTpairs1 <- c(CTpairs1, paste0(i,".",j))}}; print(CTpairs1)
for (i in g1){for (j in g2){CTpairs2 <- c(CTpairs2, paste0(j,".",i))}}; print(CTpairs2)
final_cell <- c(CTpairs1)

final_LR <- c("FLT1 complex_VEGFA","FLT1 complex_VEGFB","KDR_VEGFC","ADRB2_VEGFB", "FN1_a3b1 complex","LAMC1_a6b1 complex","LAMP1_FAM3C")

df$LR <- factor(df$LR, levels = final_LR); df$CellType_pairs <- factor(df$CellType_pairs, levels = final_cell)
df$Group <- factor(df$Group, levels = group.names); df$CellGroup <- factor(df$CellGroup, levels = lapply(final_cell, function(x) paste(x,combn(group.names,1),sep = "_")) %>% unlist)
colorMin <- min(df$logExpMeans); colorMax <- max(df$logExpMeans); sizeMin <- min(df$logP); sizeMax <- max(df$logP)
Nz <- length(final_cell); Ng <- length(group.names)

p <- ggplot(df, aes(x = CellGroup, y = LR, color = logExpMeans, size = logP)) + geom_point() + labs(x="",y="", title= "EC centered_Angiogenesis") + 
  scale_color_gradientn(colours = colorRampPalette(c("MidnightBlue", "Wheat1", "Firebrick4"))(200), limits = c(floor(colorMin),ceiling(colorMax)), name = "log2(MeanExp)") +  # viridis::magma(20)
  scale_size(name = "-log10(p)", limits = c(0, ceiling(sizeMax))) + 
  geom_vline(xintercept = (1:(Nz-1)) * Ng + 0.5, linetype="dashed", col="black", size=0.5) + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6), axis.text.y = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, size = 20), legend.position = "right")
ggsave(p, filename = "Final_EC.jpg", dpi = 500, width = 16.5, height = 7)


####################  
# matrix modeling
df <- read.table(file = paste0("./cellp2_015", "_Fib_ECM",".csv"), header = T, sep = ",")

g1 <- c("CAF"); g2 <- c("MaligantEpi","EndoBlood")
CTpairs1 <- c(); CTpairs2 <- c()
for (i in g1){for (j in g2){CTpairs1 <- c(CTpairs1, paste0(i,".",j))}}; print(CTpairs1)
for (i in g1){for (j in g2){CTpairs2 <- c(CTpairs2, paste0(j,".",i))}}; print(CTpairs2)
final_cell <- c(CTpairs1)

final_LR <- c("COL1A2_a1b1 complex","COL4A1_a2b1 complex","COL4A1_a10b1 complex","COL4A2_a2b1 complex","COL4A2_a10b1 complex")

df$LR <- factor(df$LR, levels = final_LR); df$CellType_pairs <- factor(df$CellType_pairs, levels = final_cell)
df$Group <- factor(df$Group, levels = group.names); df$CellGroup <- factor(df$CellGroup, levels = lapply(final_cell, function(x) paste(x,combn(group.names,1),sep = "_")) %>% unlist)
colorMin <- min(df$logExpMeans); colorMax <- max(df$logExpMeans); sizeMin <- min(df$logP); sizeMax <- max(df$logP)
Nz <- length(final_cell); Ng <- length(group.names)
p <- ggplot(df, aes(x = CellGroup, y = LR, color = logExpMeans, size = logP)) + geom_point() + labs(x="",y="", title= "Fib centered_ECM") + 
  scale_color_gradientn(colours = colorRampPalette(c("MidnightBlue", "Wheat1", "Firebrick4"))(200), limits = c(floor(colorMin),ceiling(colorMax)), name = "log2(MeanExp)") +  # viridis::magma(20)
  scale_size(name = "-log10(p)", limits = c(0, ceiling(sizeMax))) + 
  geom_vline(xintercept = (1:(Nz-1)) * Ng + 0.5, linetype="dashed", col="black", size=0.5) + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6), axis.text.y = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, size = 20), legend.position = "right")
ggsave(p, filename = "Final_Fib1.jpg", dpi = 500, width = 7.5, height = 5.5)


# immune cell homing
df <- read.table(file = paste0("./cellp2_015","_Fib_ImuHoming",".csv"), header = T, sep = ",")

g1 <- c("CAF"); g2 <- c("CD8T","CD4Treg","CD4Th","NK","DC")
CTpairs1 <- c(); CTpairs2 <- c()
for (i in g1){for (j in g2){CTpairs1 <- c(CTpairs1, paste0(i,".",j))}}; print(CTpairs1)
for (i in g1){for (j in g2){CTpairs2 <- c(CTpairs2, paste0(j,".",i))}}; print(CTpairs2)
final_cell <- c(CTpairs2)

final_LR <- c("CXCR3_CCL19","CXCR6_CXCL16","CCR7_CCL19","CXCR3_CXCL9")

df$LR <- factor(df$LR, levels = final_LR); df$CellType_pairs <- factor(df$CellType_pairs, levels = final_cell)
df$Group <- factor(df$Group, levels = group.names); df$CellGroup <- factor(df$CellGroup, levels = lapply(final_cell, function(x) paste(x,combn(group.names,1),sep = "_")) %>% unlist)
colorMin <- min(df$logExpMeans); colorMax <- max(df$logExpMeans); sizeMin <- min(df$logP); sizeMax <- max(df$logP)
Nz <- length(final_cell); Ng <- length(group.names)

p <- ggplot(df, aes(x = CellGroup, y = LR, color = logExpMeans, size = logP)) + geom_point() + labs(x="",y="", title= "Fib centered_ImmuHoming") + 
  scale_color_gradientn(colours = colorRampPalette(c("MidnightBlue", "Wheat1", "Firebrick4"))(200), limits = c(floor(colorMin),ceiling(colorMax)), name = "log2(MeanExp)") +  # viridis::magma(20)
  scale_size(name = "-log10(p)", limits = c(0, ceiling(sizeMax))) + 
  geom_vline(xintercept = (1:(Nz-1)) * Ng + 0.5, linetype="dashed", col="black", size=0.5) + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6), axis.text.y = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, size = 20), legend.position = "right")
ggsave(p, filename = "Final_Fib2.jpg",dpi = 500,  width = 9, height = 4.4)


####################  
# ImmuneUP, stimulatory
final_cell1 <- c("MaligantEpi.CD8T","MaligantEpi.CD4Treg","MaligantEpi.NK","MaligantEpi.monocyte","MaligantEpi.pDC","MaligantEpi.DC","MaligantEpi.macrophage","CD8T.MaligantEpi")
final_cell2 <- c("CD8T.CD4Treg","CD8T.CD4Th","CD8T.NK","CD8T.monocyte","CD8T.pDC","CD8T.DC","CD8T.macrophage")
final_cell <- c(final_cell1,final_cell2)
final_LR <- c("CD55_ADGRE5","IL6 receptor_IL6","CD48_CD244","CD28_CD80","CD28_CD86","IFNG_Type II IFNR")

df <- read.table(file = paste0("./cellp2_015", "_ImmuneUp",".csv"), header = T, sep = ",")
df$LR <- factor(df$LR, levels = final_LR); df$CellType_pairs <- factor(df$CellType_pairs, levels = final_cell)
df$Group <- factor(df$Group, levels = group.names); df$CellGroup <- factor(df$CellGroup, levels = lapply(final_cell, function(x) paste(x,combn(group.names,1),sep = "_")) %>% unlist)
colorMin <- min(df$logExpMeans); colorMax <- max(df$logExpMeans); sizeMin <- min(df$logP); sizeMax <- max(df$logP)
Nz <- length(final_cell); Ng <- length(group.names)
p <-  ggplot(df, aes(x = CellGroup, y = LR, color = logExpMeans, size = logP)) + geom_point() + labs(x="",y="", title= "ImmuneStimulatory") + 
  scale_color_gradientn(colours = colorRampPalette(c("MidnightBlue", "Wheat1", "Firebrick4"))(200), limits = c(floor(colorMin),ceiling(colorMax)), name = "log2(MeanExp)") +  # viridis::magma(20)
  scale_size(name = "-log10(p)", limits = c(0, ceiling(sizeMax))) + 
  geom_vline(xintercept = (1:(Nz-1)) * Ng + 0.5, linetype="dashed", col="black", size=0.5) + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6), axis.text.y = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, size = 20), legend.position = "right")
ggsave(p, filename = "Final_ImmuneUp.jpg", dpi = 500, width = 23, height = 6.5)


####################  
# ImmuneDown, ImmInhibitory
df <- read.table(file = paste0("./cellp2_015","_ImmuneDown1",".csv"), header = T, sep = ",")
final_cell <- c("MaligantEpi.CD8T","MaligantEpi.CD4Treg","MaligantEpi.CD4Th","MaligantEpi.NK","MaligantEpi.monocyte","MaligantEpi.pDC","MaligantEpi.DC","MaligantEpi.macrophage")
final_LR <- c("NECTIN3_TIGIT","SIRPA_CD47","PVR_CD96","PVR_TIGIT","FAM3C_PDCD1","CD99_PILRA")

df$LR <- factor(df$LR, levels = final_LR); df$CellType_pairs <- factor(df$CellType_pairs, levels = final_cell)
df$Group <- factor(df$Group, levels = group.names); df$CellGroup <- factor(df$CellGroup, levels = lapply(final_cell, function(x) paste(x,combn(group.names,1),sep = "_")) %>% unlist)
colorMin <- min(df$logExpMeans); colorMax <- max(df$logExpMeans); sizeMin <- min(df$logP); sizeMax <- max(df$logP)
Nz <- length(final_cell); Ng <- length(group.names)
p <- ggplot(df, aes(x = CellGroup, y = LR, color = logExpMeans, size = logP)) + geom_point() + labs(x="",y="", title= "ImmuneInhibition") + 
  scale_color_gradientn(colours = colorRampPalette(c("MidnightBlue", "Wheat1", "Firebrick4"))(200), limits = c(floor(colorMin),ceiling(colorMax)), name = "log2(MeanExp)") +  # viridis::magma(20)
  scale_size(name = "-log10(p)", limits = c(0, ceiling(sizeMax))) + 
  geom_vline(xintercept = (1:(Nz-1)) * Ng + 0.5, linetype="dashed", col="black", size=0.5) + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6), axis.text.y = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, size = 20), legend.position = "right")
ggsave(p, filename = "Final_ImmuneDown.jpg", dpi = 500, width = 14, height = 6)
