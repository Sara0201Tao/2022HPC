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
  
  library(Seurat)
})


#############################################
################# preparation ###############
#############################################
set.seed(1)
wkdir = "/home/mztao/scRNA_analysis/data_analysis/HPC/results/"; setwd(wkdir)
source("/home/mztao/scRNA_analysis/data_analysis/HPC/scripts/common_functions.r")

dir.create("./Cellp/")
setwd("./Cellp/")



###########################################
# load and update object
###########################################
current.object <- readRDS("../0_data/AllSample_final_ThreeGroups_Normout.RDS")
current.name <- "AllSample_final_cellp"

Cellp1 <- current.object$detailed1

Cellp2 <- plyr::revalue(current.object$detailed2, c("MaligantEpi" = "MaligantEpi","EndoLym" = "EndoLym","EndoBlood1" = "EndoBlood","EndoBlood2"="EndoBlood",
                                                    "proFib" = "proFib", "MyoFib" = "MyoFib", "CAF1" = "CAF", "CAF2" = "CAF", "CAF3" = "CAF",
                                                    "CD8T_naive" = "CD8T","CD8T_memory" = "CD8T", "CD8T_cytoxic1" = "CD8T","CD8T_cytoxic2" = "CD8T","CD8T_cytoxic3"="CD8T","CD8T_exhaust" = "CD8T","CD8T_Ebo"="CD8T",
                                                    "CD4Th_naive"="CD4Th","CD4Th_act"="CD4Th","CD4Th_exhaust"="CD4Th","CD4Treg_naive"="CD4Treg","CD4Treg_act"="CD4Treg","CD4Treg_Ebo"="CD4Treg","NK"="NK",
                                                    "B_GC"="Bcell","B_MemoInter"="Bcell","B_Memory"="Bcell","B_PlasamIgA"="Bcell","B_PlasamIgG"="Bcell",
                                                    "neutrophil"="neutrophil","mast"="mast","monocyte"="monocyte","pDC"="pDC","moDC"="DC","cDC1"="DC","cDC2"="DC","macrophageM1"="macrophage","macrophageM2"="macrophage"))
Cellp2 <- factor(Cellp2, levels = c("MaligantEpi","EndoLym","EndoBlood","proFib","MyoFib","CAF","CD8T","CD4Treg","CD4Th","NK","Bcell",
                                    "neutrophil","mast","monocyte","pDC","DC","macrophage"))


current.object$Cellp1 <- factor(Cellp1); table(current.object$Cellp1)
current.object$Cellp2 <- factor(Cellp2); table(current.object$Cellp2)
saveRDS(current.object, file = "../0_data/AllSample_final_ThreeGroups_Normout.RDS")




###########################################
# prepare for cellphoneDB running
###########################################
for (sub in as.character(unique(current.object$sample.sixgroups))){ 
  sub.object <- subset(current.object, subset = sample.sixgroups == sub)
  sub.name <- as.character(sub)
  print(sub.name); print(dim(sub.object)); print(table(sub.object$sample.samples))
  
  count_raw <- Matrix(as.matrix(GetAssayData(sub.object, slot = "counts")), sparse = TRUE); print(dim(count_raw))
  count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000);dim(count_norm)
  write.table(count_norm, paste0("./", sub.name, "_Ncount.csv"), sep=",", quote=F) # vim写入“Gene”

  meta_data2 <- data.frame(sub.object$Cellp2); colnames(meta_data2) <- c("cell_type")
  write.table(meta_data2, paste0("./", sub.name, "_meta_cellp2.csv"), sep=",", quote=F, row.names=T) # vim写入“Cell”
}

