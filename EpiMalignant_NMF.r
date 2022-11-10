library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggthemr)
library(patchwork)
library(Matrix)
library(stats)
library(stringr)

library(Seurat)
library(NMF)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db) 

set.seed(1)
source("./common_functions.r")
dir.create("./step3_Epi_NMF/")
setwd("./step3_Epi_NMF/")

# load data
Cancer.object <- readRDS("../AllSample_Mal.RDS")
dim(Cancer.object);table(Cancer.object$orig.ident)

################################################
##### getting All_MetaGenes for pre-test #######
################################################
dir.create("./rank_relative_All_MetaGenes/")

for (specified_rank in 4:7){
  dir.create(paste0("./rank_relative_All_MetaGenes/rank",specified_rank,"/"))
  for (sample in unique(Cancer.object$orig.ident)){
    SUBtumors.object <- subset(Cancer.object, subset = orig.ident == sample)
    SUBtumors.NormedData <- as.matrix(GetAssayData(SUBtumors.object, slot = "data"))
    SUBtumors.ReExpr <- ToReExpr(SUBtumors.NormedData, NoNeg = TRUE, clock_inter = 5000)
    
    # 1. fitering for genes (SD < 0.5)
    gene_SDs <- Row_SDs(SUBtumors.NormedData, clock_inter = 5000)
    
    SD_Threshold <- 0.5
    selected_genes <- rownames(SUBtumors.NormedData)[which(gene_SDs >= SD_Threshold)]
    SUBtumors.ReExpr.NMF <- SUBtumors.ReExpr[selected_genes,]
    dim(SUBtumors.ReExpr.NMF)
    
    # calaculating NMF
    res <- nmf(SUBtumors.ReExpr.NMF, rank = specified_rank, nrun = 100, .options='vtp5')
    saveRDS(res, file = paste0("./rank_relative_All_MetaGenes/rank",specified_rank,"/", sample, "_NMFres_rank",specified_rank,".RDS"))
  }
}



##### rank selection
samples <- c("T0-NonCD45","T1_1","T2","T3_1","T4_1","T4_2","T5_1","T5_2","T6_1","T6_2","T7-2")

file.dir <- "./Epi_NMF_RawResults/rank_relative_All_MetaGenes/"
for (sample in samples){# sample <- samples[1]
  print(sample)
  estim.c <- c()
  for (r in 4:7){ 
    temp <- readRDS(paste0(file.dir, "rank", r, "/",sample, "_NMFres_rank",r,".RDS"))
    estim.c <- c(estim.c,summary(temp)[["cophenetic"]])
  }
  estim.info <- data.frame(cophenetic = estim.c, rank = 4:7)
  write.table(estim.info, file = paste0("./", sample, "_NormedData_estimR.csv"),sep = ",", row.names = F, col.names = T, quote = F)
}

p <- list()
for (sample in samples){# sample <- samples[1]
  print(sample)
  estim.info <- read.table(file = paste0("../", sample, "_NormedData_estimR.csv"), sep = ",", header = T)
  temp_p <- ggplot(estim.info, aes(x = rank, y= cophenetic)) + geom_line() + labs(x = "rank",y = "cophenetic", title = sample) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  p <- c(p, list(temp_p))
}

final_p <- (p[[1]]|p[[2]]|p[[3]]|p[[4]]|p[[5]]) / (p[[6]]|p[[7]]|p[[8]]|p[[9]]|p[[10]]|p[[11]])
ggsave(filename = "./rank_summary.jpg", final_p , dpi = 300, width = 15, height = 5)








####################################
######### NMF modules defintion
####################################
##### preparing
dir.create("./Epi_NMF_PCCfunction/"); current.dir <- paste0(getwd(),"/")
NMF_types <- "relative"


NMF_types_current <- "relative"
Samples <- unique(CancerThreeGroups.object$orig.ident); print(Samples)
AllSampleRanks <- c(6,6,5,5,6,6,5,6,5,6,5) # ranking number
identical(length(Samples), length(AllSampleRanks))

for (i in 1:length(Samples)){
  ## load data
  print(Samples[i])
  current.res <- readRDS(paste0("./Epi_NMF_RawResults/rank_relative_All_MetaGenes/rank", AllSampleRanks[i], "/", Samples[i], "_NMFres_rank",AllSampleRanks[i], ".RDS"))
  raw_w <- basis(current.res); scaledR_w <- NMF:::scale_mat(raw_w, "r1")

  meta_genes <- Convert2Genes(raw_w, extractFeatures(current.res))
  genes <- c(); metainfo <- c(); Wscores_raw <- c(); Wscores_scaledR <- c()
  for (m in 1:length(meta_genes)){
    for (item in meta_genes[[m]]){
      genes <- c(genes, item)
      metainfo <- c(metainfo, paste0("MetaGene",m))
      Wscores_raw <- c(Wscores_raw, as.numeric(raw_w[item,m]))
      Wscores_scaledR <- c(Wscores_scaledR, as.numeric(scaledR_w[item,m]))}
  }
  MetaGenesInfo <- data.frame(genes = genes, metainfo = metainfo, Wscores_raw=Wscores_raw, Wscores_scaledR=Wscores_scaledR)
  write.table(MetaGenesInfo, file = paste0(current.dir, Samples[i], "_MetaGenes.csv"), sep = ",", row.names = F, col.names = T, quote = F)
}

AllSampleRanks_new <- c()
for (i in 1:length(Samples)){
  temp_metagenes <- read.table(file = paste0(current.dir, Samples[i], "_MetaGenes.csv"), sep = ",", header = T, stringsAsFactors = F)
  temp_rank <- length(unique(temp_metagenes$metainfo))
  AllSampleRanks_new <- c(AllSampleRanks_new, temp_rank)
}


AllMeta_tops_genes <- list()
for (i in 1:length(Samples)){
  temp_metagenes <- read.table(file = paste0(current.dir, Samples[i], "_MetaGenes.csv"), sep = ",", header = T, stringsAsFactors = F)
  sample_meta_tops_genes <- vector(mode = "list", length = AllSampleRanks_new[i]) 
  for (m in 1:nrow(temp_metagenes)){
    meta_group <- as.integer(substr(temp_metagenes$metainfo[m],9,9))
    sample_meta_tops_genes[[meta_group]] <- c(sample_meta_tops_genes[[meta_group]], temp_metagenes$genes[m])
  }
  names(sample_meta_tops_genes) <- paste0("Target_", Samples[i], "_MetaGene", 1:AllSampleRanks_new[i])
  AllMeta_tops_genes <- c(AllMeta_tops_genes, sample_meta_tops_genes)
}

### control gene sets
bins = 25; sizeN = 100
Idents(CancerThreeGroups.object) <- CancerThreeGroups.object$orig.ident
AllMeta_control_genes <- list()
for (i in 1:length(Samples)){
  temp_sample <- Samples[i]; print(temp_sample)
  sub.object <- subset(CancerThreeGroups.object, idents = temp_sample) 
  AllTumors.NormedData <- as.matrix(GetAssayData(sub.object, slot = "data"))

  data.avg <- Matrix::rowMeans(AllTumors.NormedData)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, n = bins, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)

  temp_metagenes <- read.table(file = paste0(current.dir, Samples[i], "_MetaGenes.csv"), sep = ",", header = T, stringsAsFactors = F)
  sample_meta_tops_genes <- vector(mode = "list", length = AllSampleRanks_new[i])
  for (m in 1:nrow(temp_metagenes)){
    meta_group <- as.integer(substr(temp_metagenes$metainfo[m],9,9))
    sample_meta_tops_genes[[meta_group]] <- c(sample_meta_tops_genes[[meta_group]], temp_metagenes$genes[m])}

  sample_control_genes <- vector(mode = "list", length = AllSampleRanks_new[i])
  for (m in 1:AllSampleRanks_new[i]) {
    features.use <- sample_meta_tops_genes[[m]]
    for (n in 1:length(features.use)) {
      random_genes <- sample(data.cut[which(data.cut == data.cut[features.use[n]])], size = sizeN, replace = FALSE)
      sample_control_genes[[m]] <- c(sample_control_genes[[m]], names(random_genes))}}
  sample_control_genes <- lapply(X = sample_control_genes, FUN = unique)
  names(sample_control_genes) <- paste0("Control_", temp_sample, "_MetaGene", 1:AllSampleRanks_new[i])
  
  AllMeta_control_genes <- c(AllMeta_control_genes, sample_control_genes)
}
saveRDS(AllMeta_tops_genes, file = paste0(current.dir, "AllMetaTopsGenes.RDS"))
saveRDS(AllMeta_control_genes, file = paste0(current.dir, "AllMetaControlGenes.RDS"))

### prepare data for PCCs
Idents(CancerThreeGroups.object) <- CancerThreeGroups.object$orig.ident
sub.object <- subset(CancerThreeGroups.object, idents = Samples) 
dim(sub.object);table(sub.object$orig.ident)
Sub_NormedData <- as.matrix(GetAssayData(sub.object, slot = "data"))
Sub_ReExpr <- ToReExpr(Sub_NormedData, NoNeg = TRUE, clock_inter = 4000)

# calculating cellscores from relative data
AllMetaCellScores <- data.frame()
all_metas <- length(AllMeta_tops_genes)
for (group in 1:all_metas){print(group);
  temp <- CellScores(Sub_ReExpr, AllMeta_tops_genes[[group]], AllMeta_control_genes[[group]])
  AllMetaCellScores <- rbind(AllMetaCellScores, temp)}
colnames(AllMetaCellScores) <- colnames(sub.object); rownames(AllMetaCellScores) <- names(AllMeta_tops_genes); print(dim(t(AllMetaCellScores)))
write.table(t(AllMetaCellScores), file = paste0(current.dir,"PCC_data_relative.csv"),sep = ",", row.names = T, col.names = T, quote = F)

### PCCs calculating for clustering
PCC_data <- as.matrix(read.table(file = paste0(current.dir,"PCC_data_relative.csv"), sep = ",", row.names = 1, header = T, stringsAsFactors = F))
cor_mat <- cor(PCC_data, method = "pearson")
annotations <- data.frame(Groups = c(rep(sub.names[1], sum(AllSampleRanks_new[1:3])), rep(sub.names[2], sum(AllSampleRanks_new[4:7])), rep(sub.names[3], sum(AllSampleRanks_new[8:11]))))
rownames(annotations) <- rownames(cor_mat)
pheatmap(cor_mat, show_rownames=1, show_colnames=F, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), annotation_row = annotations,
         scale = "none", fontsize_row=6, cluster_rows = T, cluster_cols = T, # display_numbers = T, fontsize_number = 7,
         filename=paste0(current.dir,"PCC_heatmap_relative.jpg"), width=12, height=8)
