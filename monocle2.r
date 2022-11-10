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



current.object <- readRDS("./0_data/AllSample_CD8T.RDS")
current.name <- "AllSample_CD8T"
sample_info <- data.frame(id = c("T1_2",  "NT4","NT5",  "T0_1","T1_1","T3_1",  "T2_1","T4_1","T5_1","T6_1",  "T3_2", "T4_2","T5_2","T6_2","T7_2"),
                          group = c("Lym","NT","NT","NaiveGood","NaiveGood","NaiveGood","NaiveBad","NaiveBad","NaiveBad","NaiveBad","TreatGood","TreatBad","TreatBad","TreatBad","TreatBad"))

## monocle2
library(monocle)

data <- as(as.matrix(current.object@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = current.object@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = rownames(data), row.names = rownames(data)))
cds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size(), lowerDetectionLimit = 1)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, cores=20)

disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.05 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
plot1 <- plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State, nrow = 1)
plot2 <- plot_cell_trajectory(cds, color_by = "detailed3") + facet_wrap(~detailed3, nrow = 1)
p <- plot1/plot2; p

cds <- orderCells(cds, root_state = c(3)) # 确定起始点
plot1 <- plot_cell_trajectory(cds, color_by = "State")
plot2 <- plot_cell_trajectory(cds, color_by = "Pseudotime")
plot3 <- plot_cell_trajectory(cds, color_by = "detailed3")
plot4 <- plot_cell_trajectory(cds, color_by = "detailed3") + facet_wrap(~detailed3, nrow = 1)
p <- (plot1|plot2|plot3) / plot4; p

current.object$monocle2_state <- cds@phenoData@data[["State"]]
current.object$monocle2_time <- cds@phenoData@data[["Pseudotime"]]
current.object$monocle2_comp1 <- plot_cell_trajectory(cds)[["data"]][["data_dim_1"]]
current.object$monocle2_comp2 <- plot_cell_trajectory(cds)[["data"]][["data_dim_2"]]

saveRDS(current.object, file = paste0(current.dir, current.name, ".RDS"))
saveRDS(cds, file = paste0(current.dir, current.name, "_monocle.RDS"))


df1 <- data.frame(current.object@meta.data); df1$ID <- rownames(df1)
df2 <- as.data.frame(t(GetAssayData(current.object, slot = "data")[c("TCF7","GZMK","CTLA4","MKI67"),])); df2$ID <- rownames(df2)
df <- merge(df1, df2, by = "ID")
df$detailed3 <- factor(df$detailed3, levels = levels(current.object$detailed3))
df$sample.locations <- factor(df$sample.locations, levels = levels(current.object$sample.locations))
df$sample.sixgroups<- factor(df$sample.sixgroups, levels = levels(current.object$sample.sixgroups))


p1 <- ggplot(df, aes(x = monocle2_comp1, y = monocle2_comp2, color = detailed3)) + geom_point(size = 1.5) + scale_color_igv() + labs(x="Monocle2_Comp1", y="Monocle2_Comp2", title= "") +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5), axis.text = element_text(size = 15), legend.position = "top", plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) 

p2 <- ggplot(df, aes(x = monocle2_comp1, y = monocle2_comp2, color = TCF7)) + geom_point(size = 1) + labs(x="", y="", title= "TCF7") +
  scale_color_gradientn(colours = colorRampPalette(c("grey", "Firebrick3"))(50), name = "Exp") + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) 
p3 <- ggplot(df, aes(x = monocle2_comp1, y = monocle2_comp2, color = GZMK)) + geom_point(size = 1) + labs(x="", y="", title= "GZMK") +
  scale_color_gradientn(colours = colorRampPalette(c("grey", "Firebrick3"))(50), name = "Exp") + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) 
p4 <- ggplot(df, aes(x = monocle2_comp1, y = monocle2_comp2, color = CTLA4)) + geom_point(size = 1) + labs(x="", y="", title= "CTLA4") +
  scale_color_gradientn(colours = colorRampPalette(c("grey", "Firebrick3"))(50), name = "Exp") + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) 
p5 <- ggplot(df, aes(x = monocle2_comp1, y = monocle2_comp2, color = MKI67)) + geom_point(size = 1) + labs(x="", y="", title= "MKI67") +
  scale_color_gradientn(colours = colorRampPalette(c("grey", "Firebrick3"))(50), name = "Exp") + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) 
p <- p1 | ((p2|p3) / (p4|p5))
ggsave(p, filename = paste0(current.dir, current.name, "_monocle.jpg"), width = 12, height = 8)


