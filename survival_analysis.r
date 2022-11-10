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
  
  library(Seurat)
  library(pheatmap)
})


#############################################
################# preparation ###############
#############################################
set.seed(1)
wkdir = "/home/mztao/scRNA_analysis/data_analysis/HPC/results/"; setwd(wkdir)
source("/home/mztao/scRNA_analysis/data_analysis/HPC/scripts/common_functions.r")


#############################################
### survival analysis
#############################################
current.dir <- paste0(getwd(), "/")

library(GSVA)
library(survival)
library(survminer)
library(patchwork)

## data prepare
# HPC
SurvDataDir <- "/home/mztao/scRNA_analysis/data_analysis/HPC/results/0_data/bulkRNA/"
SurvData_bulk44_gene <- read.table(paste0(SurvDataDir, "bulkRNA_Survival44_exp.csv"), header = T, sep = ",", stringsAsFactors = F)
SurvData_bulk44_gene <- as.data.frame(t(SurvData_bulk44_gene))
SurvData_bulk44_info <- read.table(paste0(SurvDataDir, "bulkRNA_Survival44_SampleInfo.csv"), header = T, sep = ",", stringsAsFactors = F)

Naive30 <- SurvData_bulk44_info[SurvData_bulk44_info$Info %in% c("NaiveGood","NaiveBad"), "SampleID"]
SurvData_bulk30_gene <- SurvData_bulk44_gene[Naive30,]
SurvData_bulk30_info <- SurvData_bulk44_info[SurvData_bulk44_info$SampleID %in% Naive30, ]

# NPC88. gene names in colnames with "-" is turn to be "."
SurvData_NPC <- read.table(file = "/home/mztao/scRNA_analysis/sup_info/PublicData/cohort113/NPC88_Survival.csv", sep = ",", header = T, stringsAsFactors = F)
info_col <- c(1,((ncol(SurvData_NPC))-3):(ncol(SurvData_NPC))); SurvData_NPC_info <- SurvData_NPC[,info_col]
gene_col <- setdiff(1:ncol(SurvData_NPC), info_col); SurvData_NPC_gene <- SurvData_NPC[,gene_col]
delete_col <- c(); for (i in 1:ncol(SurvData_NPC_gene)){temp <- SurvData_NPC_gene[,i];if (length(which(temp == 0)) >= 88-50 ){delete_col <- c(delete_col, i)}}
SurvData_NPC_gene <- SurvData_NPC_gene[,-delete_col]; rm(SurvData_NPC)


# gene CDKN2A
test_gene <- "CDKN2A"
suvData_bulk <- SurvivalData(datasetInfo = SurvData_bulk44_info, datasetExp = SurvData_bulk44_gene, test_gene = test_gene)
temp <- survdiff(Surv(Time, Status) ~ MedianType, data = suvData_bulk, rho = 0); km_p_median <- 1-pchisq(temp$chisq, length(temp$n) -1)
p <- ggsurvplot(survfit(Surv(Time, Status) ~ MedianType, data = suvData_bulk), data = suvData_bulk, pval = T, conf.int = F, xlab = 'Time(Days)', ylab = "Survival", title = paste(test_gene, "Bulk44", "median", sep = "_"),
                break.time.by = 400, legend.labs = c('low','High'), palette = c("blue","red"), 
                ggtheme = theme(panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, face = "bold")),
                risk.table = T, ncensor.plot = F)
ggsave(p$plot / plot_spacer()/ p$table + plot_layout(heights = c(4, 0.5, 2)), filename = paste0(current.dir, "CDKN2A_bulk44_median.jpg"), width = 6, height = 8)

suvData_bulk <- SurvivalData(datasetInfo = SurvData_NPC_info, datasetExp = SurvData_NPC_gene, test_gene = gsub("-",".",test_gene))
temp <- survdiff(Surv(Time, Status) ~ MedianType, data = suvData_bulk, rho = 0); km_p_median <- 1-pchisq(temp$chisq, length(temp$n) -1)
p <- ggsurvplot(survfit(Surv(Time, Status) ~ MedianType, data = suvData_bulk), data = suvData_bulk, pval = T, conf.int = F, xlab = 'Time(Months)', ylab = "Progress free survival", title = paste(test_gene, "NPC88", "median", sep = "_"),
                break.time.by = 10, legend.labs = c('low','High'), palette = c("blue","red"), 
                ggtheme = theme(panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, face = "bold")),
                risk.table = T, ncensor.plot = F)
ggsave(p$plot / plot_spacer()/ p$table + plot_layout(heights = c(4, 0.5, 2)), filename = paste0(current.dir, "CDKN2A_NPC88_median.jpg"), width = 6, height = 8)




# NMF module genes
# M3, EpiDev_modules
temp_genes <- c("GRHL3","OVOL1","DSC2","A2ML1","ZNF750","SERPINB2","RHCG","HOPX","PPL","SCEL","SPRR3","KRT19","KRT15", "CSTB", "SRD5A1", "CSTA", "IVL", "EMP1", "TRIM29")
genelist <- list(temp_genes); names(genelist) <- "EpiDev"

temp <- as.matrix(t(SurvData_bulk30_gene)); colnames(temp) <- SurvData_bulk30_info$SampleID
re_gsva <- gsva(temp, genelist, method = "gsva", kcdf="Gaussian", min.sz=2, max.sz = 500, parallel.sz=15, mx.diff=TRUE, verbose=TRUE) # gsva算法，表达值为data时，Gaussian分布

suvData2 <- SurvivalData(datasetInfo = SurvData_bulk30_info, datasetExp = t(re_gsva), test_gene = names(genelist))
p <- ggsurvplot(survfit(Surv(Time, Status) ~ MedianType, data = suvData2), data = suvData2, pval = T, conf.int = F, risk.table = T, ncensor.plot = F, censor = T, 
                  xlab = 'Time(Days)', ylab = "survival", title = paste(names(genelist), "Scores", "bulk30", "median", "gsva", sep = "_"),
                  break.time.by = 400, legend.labs = c('low','High'), palette = c("blue","red"),
                  ggtheme = theme(panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, size = 12, face = "bold")))
ggsave(p$plot / plot_spacer()/ p$table + plot_layout(heights = c(4, 0.5, 2)), filename = paste0(current.dir, "EpiDev_bulk30_median.jpg"), width = 6, height = 8)

p <- ggsurvplot(survfit(Surv(Time, Status) ~ MeanType, data = suvData2), data = suvData2, pval = T, conf.int = F, risk.table = T, ncensor.plot = F, censor = T, 
                  xlab = 'Time(Days)', ylab = "survival", title = paste(names(genelist), "Scores", "bulk30", "mean", "gsva", sep = "_"),
                  break.time.by = 400, legend.labs = c('low','High'), palette = c("blue","red"),
                  ggtheme = theme(panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, size = 12, face = "bold")))
ggsave(p$plot / plot_spacer()/ p$table + plot_layout(heights = c(4, 0.5, 2)), filename = paste0(current.dir, "EpiDev_bulk30_mean.jpg"), width = 6, height = 8)



# M2, CellCycle_module
temp_genes <- c("GMNN","SKA2","CKS1B","RRM1","CCNB1","ZWINT","PBK","G2E3","DEK","MLLT3","SMC2","NASP","USP1","SMC4","PSMA6","PSMD13","PSME2","CD44")
genelist <- list(temp_genes); names(genelist) <- "CellCycle"

temp <- as.matrix(t(SurvData_bulk30_gene)); colnames(temp) <- SurvData_bulk30_info$SampleID
re_gsva <- gsva(temp, genelist, method = "gsva", kcdf="Gaussian", min.sz=2, max.sz = 500, parallel.sz=15, mx.diff=TRUE, verbose=TRUE) # gsva算法，表达值为data时，Gaussian分布

suvData2 <- SurvivalData(datasetInfo = SurvData_bulk30_info, datasetExp = t(re_gsva), test_gene = names(genelist))
p <- ggsurvplot(survfit(Surv(Time, Status) ~ MedianType, data = suvData2), data = suvData2, pval = T, conf.int = F, risk.table = T, ncensor.plot = F, censor = T, 
                  xlab = 'Time(Days)', ylab = "survival", title = paste(names(genelist), "Scores", "bulk30", "median", "gsva", sep = "_"),
                  break.time.by = 400, legend.labs = c('low','High'), palette = c("blue","red"),
                  ggtheme = theme(panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, size = 12, face = "bold")))
ggsave(p$plot / plot_spacer()/ p$table + plot_layout(heights = c(4, 0.5, 2)), filename = paste0(current.dir, "Cellcycle_bulk30_median.jpg"), width = 6, height = 8)

p <- ggsurvplot(survfit(Surv(Time, Status) ~ MeanType, data = suvData2), data = suvData2, pval = T, conf.int = F, risk.table = T, ncensor.plot = F, censor = T, 
                  xlab = 'Time(Days)', ylab = "survival", title = paste(names(genelist), "Scores", "bulk30", "mean", "gsva", sep = "_"),
                  break.time.by = 400, legend.labs = c('low','High'), palette = c("blue","red"),
                  ggtheme = theme(panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, size = 12, face = "bold")))
ggsave(p$plot / plot_spacer()/ p$table + plot_layout(heights = c(4, 0.5, 2)), filename = paste0(current.dir, "Cellcycle_bulk30_mean.jpg"), width = 6, height = 8)


# M4, EMT_module
temp_genes <- c("AKR1B10","GSTO1","SEC61G","C12orf75","ATR","LAMC1","IGF1R","ITGA6","MET","PVR","PTEN","PAWR","HSPA4","SPP1","ITGA2","INHBB","MAP1B","LIFR","WNT5A","PGAM1")
genelist <- list(temp_genes); names(genelist) <- "EMT"

temp <- as.matrix(t(SurvData_bulk30_gene)); colnames(temp) <- SurvData_bulk30_info$SampleID
re_gsva <- gsva(temp, genelist, method = "gsva", kcdf="Gaussian", min.sz=2, max.sz = 500, parallel.sz=15, mx.diff=TRUE, verbose=TRUE) # gsva算法，表达值为data时，Gaussian分布

suvData2 <- SurvivalData(datasetInfo = SurvData_bulk30_info, datasetExp = t(re_gsva), test_gene = names(genelist))
p <- ggsurvplot(survfit(Surv(Time, Status) ~ MedianType, data = suvData2), data = suvData2, pval = T, conf.int = F, risk.table = T, ncensor.plot = F, censor = T, 
                  xlab = 'Time(Days)', ylab = "survival", title = paste(names(genelist), "Scores", "bulk30", "median", "gsva", sep = "_"),
                  break.time.by = 400, legend.labs = c('low','High'), palette = c("blue","red"),
                  ggtheme = theme(panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, size = 12, face = "bold")))
ggsave(p$plot / plot_spacer()/ p$table + plot_layout(heights = c(4, 0.5, 2)), filename = paste0(current.dir, "EMT_bulk30_median.jpg"), width = 6, height = 8)

p <- ggsurvplot(survfit(Surv(Time, Status) ~ MeanType, data = suvData2), data = suvData2, pval = T, conf.int = F, risk.table = T, ncensor.plot = F, censor = T, 
                  xlab = 'Time(Days)', ylab = "survival", title = paste(names(genelist), "Scores", "bulk30", "mean", "gsva", sep = "_"),
                  break.time.by = 400, legend.labs = c('low','High'), palette = c("blue","red"),
                  ggtheme = theme(panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, size = 12, face = "bold")))
ggsave(p$plot / plot_spacer()/ p$table + plot_layout(heights = c(4, 0.5, 2)), filename = paste0(current.dir, "EMT_bulk30_mean.jpg"), width = 6, height = 8)









