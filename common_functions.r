# Annotate genes and other RNA classes such as lncRNA, etc.
getAnnotations <- function(gtf.dir = NULL, is.mouse = FALSE) {
  library(rtracklayer)
  gtf <- import(gtf.dir)
  gtf <- gtf[gtf$gene_name!=""]
  gtf <- gtf[!is.na(gtf$gene_name)]
  biotypes <- unique(gtf$transcript_type)
  protein_coding <- gtf$gene_name[gtf$transcript_type %in% 
                                    c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", 
                                      "IG_M_gene", "IG_V_gene", "IG_Z_gene", 
                                      "nonsense_mediated_decay", "nontranslating_CDS", 
                                      "non_stop_decay", "polymorphic_pseudogene", 
                                      "protein_coding", "TR_C_gene", "TR_D_gene", "TR_gene", 
                                      "TR_J_gene", "TR_V_gene")]
  lincRNA <- gtf$gene_name[gtf$transcript_type %in% 
                             c("3prime_overlapping_ncrna", "ambiguous_orf", 
                               "antisense_RNA", "antisense", "lincRNA", "ncrna_host", "non_coding", 
                               "processed_transcript", "retained_intron", 
                               "sense_intronic", "sense_overlapping")]
  sncRNA <- gtf$gene_name[gtf$transcript_type %in% 
                            c("miRNA", "miRNA_pseudogene", "misc_RNA", 
                              "misc_RNA_pseudogene", "Mt_rRNA", "Mt_tRNA", 
                              "Mt_tRNA_pseudogene", "ncRNA", "pre_miRNA", 
                              "RNase_MRP_RNA", "RNase_P_RNA", "rRNA", "rRNA_pseudogene", 
                              "scRNA_pseudogene", "snlRNA", "snoRNA", "snRNA",
                              "snRNA_pseudogene", "SRP_RNA", "tmRNA", "tRNA",
                              "tRNA_pseudogene", "ribozyme", "scaRNA", "sRNA")]
  pseudogene <- gtf$gene_name[gtf$transcript_type %in% 
                                c("disrupted_domain", "IG_C_pseudogene", "IG_J_pseudogene", 
                                  "IG_pseudogene", "IG_V_pseudogene", "processed_pseudogene", 
                                  "pseudogene", "transcribed_processed_pseudogene",
                                  "transcribed_unprocessed_pseudogene", 
                                  "translated_processed_pseudogene", 
                                  "translated_unprocessed_pseudogene", "TR_J_pseudogene", 
                                  "TR_V_pseudogene", "unitary_pseudogene", 
                                  "unprocessed_pseudogene")]
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  diss.genes <- read.table("/home/mztao/scRNA_analysis/sup_info/diss-genes.txt", header = F, stringsAsFactors = F)$V1
  proliferation.genes <- c("MKI67","NUSAP1","PLK1","CDC20","CDK1","CDKN3","CENPA","BIRC5","PCNA","CCNA2")
  if (!is.mouse){
    mito.genes <- grep('^MT-', gtf$gene_name, value = TRUE)
    ribo.genes <- grep('^RPL|^RPS|^MRPL|^MRPS', gtf$gene_name, value = TRUE)
  } else {
    mito.genes <- grep('^mt-', gtf$gene_name, value = TRUE)
    ribo.genes <- grep('^Rpl|^Rps|^Mrpl|^Mrps', gtf$gene_name, value = TRUE)
    HB.genes <- getMouseGene(HB.genes)
    diss.genes <- getMouseGene(diss.genes)
    diss.genes.core <- getMouseGene(diss.genes.core)
    proliferation.genes <- getMouseGene(proliferation.genes)
  }
  annotations <- list(protein_coding=unique(protein_coding), 
                      pseudogene=unique(pseudogene),
                      lincRNA=unique(lincRNA),
                      sncRNA=unique(sncRNA),
                      hb=unique(HB.genes),
                      mito=unique(mito.genes),
                      ribo=unique(ribo.genes),
                      diss=unique(diss.genes),
                      diss.core=diss.genes.core,
                      proli=proliferation.genes)
  annotations <- annotations[lengths(annotations)>0]
}

# create an seurat object and calculating statiscis for plotting
scStat <- function(gex.data, name, annotations, mini_cells = 3) {
  gex <- CreateSeuratObject(counts = gex.data, project = name, min.cells = mini_cells)
  all_genes <- rownames(gex)
  gex[["percent.protein_coding"]] <- PercentageFeatureSet(gex, features = intersect(all_genes, annotations$protein_coding))
  gex[["percent.pseudogene"]] <- PercentageFeatureSet(gex, features = intersect(all_genes, annotations$pseudogene))
  gex[["percent.lincRNA"]] <- PercentageFeatureSet(gex, features = intersect(all_genes, annotations$lincRNA))
  gex[["percent.sncRNA"]] <- PercentageFeatureSet(gex, features = intersect(all_genes, annotations$sncRNA))
  gex[["percent.hb"]] <- PercentageFeatureSet(gex, features = intersect(all_genes, annotations$hb))
  gex[["percent.mito"]] <- PercentageFeatureSet(gex, features = intersect(all_genes, annotations$mito))
  gex[["percent.ribo"]] <- PercentageFeatureSet(gex, features = intersect(all_genes, annotations$ribo))
  gex[["percent.diss"]] <- PercentageFeatureSet(gex, features = intersect(all_genes, annotations$diss))
  gex[["percent.proli"]] <- PercentageFeatureSet(gex, features = intersect(all_genes, annotations$proli))
  
  filepath = paste("./", name, "_BasicStats1.jpg", sep='')
  ggsave(filepath, VlnPlot(gex, features = c("nFeature_RNA", "nCount_RNA","percent.protein_coding","percent.mito","percent.ribo", "percent.hb", "percent.diss","percent.proli"), ncol = 9), dpi = 300, width = 45)
  filepath = paste("./", name, "_BasicStats2.jpg", sep='')
  ggsave(filepath, FeatureScatter(gex, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"), dpi = 300)
  filepath = paste("./", name, "_BasicStats3.jpg", sep='')
  ggsave(filepath, FeatureScatter(gex, feature1 = "nCount_RNA", feature2 = "percent.mito"), dpi = 300)
  
  gex
}


# output for 10x stardard files from an seurat obecjt
Out10x <- function(seurat.object, wkdir = NULL, sample_name = NULL){
  raw_counts <- as.data.frame(GetAssayData(seurat.object, slot = "counts"))

  write.table(data.frame(rownames(raw_counts),rownames(raw_counts)),file = paste0(wkdir,sample_name,'_genes.tsv'), quote = F,sep = '\t',col.names = F,row.names = F)
  write.table(colnames(raw_counts), file = paste0(wkdir,sample_name,'_barcodes.tsv'), quote = F, col.names = F,row.names = F)
  
  file=paste0(wkdir,sample_name,'_matrix.mtx')
  sink(file)
  cat("%%MatrixMarket matrix coordinate integer general\n")
  cat("%\n")
  cat(paste(nrow(raw_counts),ncol(raw_counts),sum(raw_counts>0),"\n")) 
  sink()
  tmp=do.call(rbind,lapply(1:ncol(raw_counts),function(i){
    return(data.frame(row=1:nrow(raw_counts),col=i,exp=raw_counts[,i]))
  }))
  tmp=tmp[tmp$exp>0,]
  # head(tmp)
  write.table(tmp, file = paste0(wkdir,sample_name,'_matrix.mtx'), quote = F,col.names = F,row.names = F,append = T )
}


# Pie plot
PiePlotting <- function(data,SaveDir,name,LegendAdd){
  library(RColorBrewer)
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  jpeg(file=paste0(SaveDir, name, "_PiePlot.jpeg"))
  
  per <- paste(round(100 * data / sum(data),2),"%")
  ele.names <- names(data)
  ele.cols <- getPalette(length(data))
  
  pie(data, labels=per, col=ele.cols)
  if (LegendAdd == T) {legend("bottomleft", ele.names, fill = ele.cols)}
  
  dev.off()
}


# Veen plotting
VennPlot <- function(InList, InNames, mycolor, mydir, myfile){
  library(VennDiagram)
  library(grid)
  library(futile.logger)
  
  venn.diagram(x = InList, category = InNames, fill = mycolor, filename = paste0(mydir, myfile, "_venn.tiff"))
}


# convert a data.frame to a unique vector
GetAllVariables <- function(data){
  rownums <- nrow(data)
  colnums <- ncol(data)
  all_variables <- c()
  for (m in 1:rownums){
    for (n in 1:colnums){
      all_variables <- c(all_variables, as.character(data[m,n]))
    }
  }
  all_variables
}

# list convert to vector
List2Vector <- function(mylist){
  groups <- length(mylist)
  AllItems <- c()
  for (i in 1:groups){
    AllItems <- c(AllItems, mylist[[i]])
  }
  return(AllItems)
}

# calculating SD for each row of input
Row_SDs <- function(mymatrix, clock_inter = 3000){
  SDs <- c()
  genesN <- nrow(mymatrix)
  clock <- 0
  for (i in 1:genesN){
    clock <- clock + 1
    if (clock %% clock_inter == 0){print(clock)}
    
    SDs <- c(SDs, sd(as.numeric(mymatrix[i,])))
  }
  return(SDs)
}


# remove NAs in a list
ListUpdate <- function(mylist){
  NA_candidata <- which(lengths(mylist) == 1)
  NA_true <- c()
  
  for (m in NA_candidata){
    if (is.na(mylist[[m]]) == TRUE){
      NA_true <- c(NA_true, m)}
  }
  list_left <- setdiff((1:length(mylist)), NA_true)
  
  mylist_final <- list()
  for (m in list_left){
    mylist_final <- c(mylist_final, list(mylist[[m]]))}
  
  return(mylist_final)
}


# calculating cell scores with target and control gene groups per group one time
CellScores <- function(ExprMatrix, TargetGenes, ControlGenes){
  Gt_scores <- Matrix::colMeans(ExprMatrix[TargetGenes,,drop = FALSE])
  Gc_scores <- Matrix::colMeans(ExprMatrix[ControlGenes,,drop = FALSE])
  
  cell_scores <- Gt_scores - Gc_scores
  return(cell_scores)
}


# input a list of genes to do KEGG&GO in humans, gene symbols as input 
KeggGO_ORA <- function(mygenes){
  library(org.Hs.eg.db) 
  library(clusterProfiler)

  genes_df <- bitr(mygenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = TRUE) 
  colnames(genes_df) <- c("gene","EntrzID")
  
  # KEGG
  kegg.re <- enrichKEGG(gene = genes_df$EntrzID, organism  = 'hsa', keyType = "kegg",
                        pAdjustMethod = "fdr",pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500)
  if (is.null(kegg.re)) {} else {kegg.re <- setReadable(kegg.re, OrgDb = org.Hs.eg.db, keyType="ENTREZID")}
  print("kegg Done")
  
  # GO
  go.re1 <- enrichGO(gene = genes_df$EntrzID, keyType = "ENTREZID", OrgDb= org.Hs.eg.db, ont="BP", 
                     pAdjustMethod = "fdr", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, 
                     minGSSize = 10, maxGSSize = 500, readable = TRUE); print("GO1 Done")
  go.re2 <- enrichGO(gene = genes_df$EntrzID, keyType = "ENTREZID", OrgDb= org.Hs.eg.db, ont="CC", 
                     pAdjustMethod = "fdr", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, 
                     minGSSize = 10, maxGSSize = 500, readable = TRUE); print("GO2 Done")
  go.re3 <- enrichGO(gene = genes_df$EntrzID, keyType = "ENTREZID", OrgDb= org.Hs.eg.db, ont="MF", 
                     pAdjustMethod = "fdr", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, 
                     minGSSize = 10, maxGSSize = 500, readable = TRUE); print("GO3 Done")
  
  # summarize and paprare for output 
  mylist <- list(kegg.re, go.re1, go.re2, go.re3)
  names(mylist) <- c("KEGG","GO_BP","GO_CC","GO_MF")
  return(mylist)
}


## convert a gmt to list
gmt2list <- function(gmtfile){
  sets <- as.list(read_lines(gmtfile))
  for(i in 1:length(sets)){
    tmp = str_split(sets[[i]], '\t')
    n = length(tmp[[1]])
    names(sets)[i] = tmp[[1]][1]
    sets[[i]] = tmp[[1]][3:n]
    rm(tmp, n)
  }
  return(sets)
}



# top regulons sellection
GetTopRegulons <- function(data, topNum = 5){
  rownums <- nrow(data)
  colnums <- ncol(data)
  regulons <- colnames(data)
  all_tops <- c()
  for (i in 1:rownums){
    temp_data <- data[i,]
    temp_rank <- base::rank(temp_data, na.last = F)
    temp_top <- c()
    for (j in colnums:(colnums - topNum + 1)){
      temp_top <- c(temp_top, regulons[which(temp_rank == j)])}
    all_tops <- rbind(all_tops,temp_top)
  }
  rownames(all_tops) <- rownames(data)
  all_tops <- data.frame(all_tops)
  all_tops
}



# prepare geneset list for calculating 
PrepareGeneSets <- function(){
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(magrittr)
  
  hall_geneset <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(c("gs_name","gs_url","gene_symbol")) %>% as.data.frame()
  hall_geneset <- split(hall_geneset$gene_symbol, hall_geneset$gs_name)
  KEGG_geneset <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% dplyr::select(c("gs_name","gs_url","gene_symbol")) %>% as.data.frame()
  KEGG_geneset <- split(KEGG_geneset$gene_symbol, KEGG_geneset$gs_name)
  GO_geneset_BP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% dplyr::select(c("gs_name","gs_url","gene_symbol")) %>% as.data.frame()
  GO_geneset_BP <- split(GO_geneset_BP$gene_symbol, GO_geneset_BP$gs_name)
  GO_geneset_MF <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>% dplyr::select(c("gs_name","gs_url","gene_symbol")) %>% as.data.frame()
  GO_geneset_MF <- split(GO_geneset_MF$gene_symbol, GO_geneset_MF$gs_name)
  GO_geneset_CC <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>% dplyr::select(c("gs_name","gs_url","gene_symbol")) %>% as.data.frame()
  GO_geneset_CC <- split(GO_geneset_CC$gene_symbol, GO_geneset_CC$gs_name)
  meta_geneset <- read.gmt("/home/mztao/scRNA_analysis/sup_info/GeneSet/metastatic_pathways_update.gmt")
  colnames(meta_geneset) <- c("gs_name","gene_symbol")
  meta_geneset <- split(meta_geneset$gene_symbol, meta_geneset$gs_name)
  geneset_list <- list(hall_geneset,KEGG_geneset,GO_geneset_BP,GO_geneset_MF,GO_geneset_CC,meta_geneset)
  names(geneset_list) <- c("hallmark","KEGG","GO_BP","GO_MF","GO_CC","metabolic")
  
  geneset_list
}


## generate a data.frame to plot figs for cellphoneDB (dot plot)
CellPhoneDF <- function(ExpData, Pvalues, LR_pairs, CellType_pairs){
  m <- c()
  for (item in CellType_pairs) {m <- c(m, ExpData[LR_pairs, item])}
  df_exp <- data.frame(ExpMeans = m, LR = rep(LR_pairs, each = 1, time = length(CellType_pairs)), CellTypes = rep(CellType_pairs, each = length(LR_pairs), time = 1))
  
  m <- c()
  for (item in CellType_pairs) {m <- c(m, Pvalues[LR_pairs, item])}
  df_pvalue <- data.frame(pvalue = m, LR = rep(LR_pairs, each = 1, time = length(CellType_pairs)), CellTypes = rep(CellType_pairs, each = length(LR_pairs), time = 1))
  
  df <- cbind(df_exp, df_pvalue); df <- df[,-c(5,6)]
  colnames(df)[3] <- "CellType_pairs"
  
  df
}


BarPlotting_GroupType <- function(mydata, myGroup, myType, GroupLevel, TypeLevel, TypeCol, MainTitle = ""){
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  
  mydata[,myGroup] <- factor(mydata[,myGroup], levels = GroupLevel)
  mydata[,myType] <- factor(mydata[,myType], levels = TypeLevel)
  df <- as.data.frame(table(mydata[,myGroup], mydata[,myType])); colnames(df) <- c("Group","Type","Freq")
  if (is.numeric(TypeCol)) { TypeCol <- ggsci::pal_igv("default")(51)[1:TypeCol]}
  
  p1 <- ggplot(df, aes(x=Group, y=Freq, fill=Type)) + scale_fill_manual(values = TypeCol) + geom_bar(stat="identity", position = 'stack', width = 0.8) + labs(x = "", y = "Number", title = paste0(MainTitle, "1")) +  
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5), axis.text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), legend.position = "top")
  p2 <- ggplot(df, aes(x=Group, y=Freq, fill=Type)) + scale_fill_manual(values = TypeCol) + geom_bar(stat="identity", position = 'fill', width = 0.8) + labs(x = "", y = "Proportion", title = paste0(MainTitle, "2")) +  
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5), axis.text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), legend.position = "top")
  
  # return(p1|p2)
  return(list(p1, p2))
}


DensityPlot2V <- function(mydata, scoreX, scoreY, Xline, Yline, xMinMax, yMinMax, MainTitle = ""){
  library(ggplot2)
  library(tidyverse)
  
  Q1 <- sum((mydata[,scoreX] >= Xline) & (mydata[,scoreY] >= Yline))
  Q2 <- sum((mydata[,scoreX] < Xline) & (mydata[,scoreY] >= Yline))
  Q3 <- sum((mydata[,scoreX] < Xline) & (mydata[,scoreY] < Yline))
  Q4 <- sum((mydata[,scoreX] >= Xline) & (mydata[,scoreY] < Yline))
  Qs <- c(Q1, Q2, Q3, Q4); Qs <- round(Qs/sum(Qs), 4) * 100
  
  p <- ggplot(mydata, aes(mydata[,scoreX], mydata[,scoreY]) ) + geom_bin2d(bins = 150) + scale_fill_continuous(type = "viridis") + labs(x = scoreX, y = scoreY, title = MainTitle) + 
    geom_vline(xintercept = Xline, linetype="dashed", col="black", size=0.8) + geom_hline(yintercept = Yline, linetype="dashed", col="black", size=0.8) + 
    scale_x_continuous(limits = xMinMax) + scale_y_continuous(limits = yMinMax) + 
    annotate("text", x = quantile(xMinMax, probs = seq(0, 1, 0.01))[90] , y = quantile(yMinMax, probs = seq(0, 1, 0.01))[95], label = paste0("Q1: ", Qs[1], "%") , size= 5) + 
    annotate("text", x = quantile(xMinMax, probs = seq(0, 1, 0.01))[10], y = quantile(yMinMax, probs = seq(0, 1, 0.01))[95], label = paste0("Q2: ", Qs[2], "%") , size= 5) + 
    annotate("text", x = quantile(xMinMax, probs = seq(0, 1, 0.01))[10] , y = quantile(yMinMax, probs = seq(0, 1, 0.01))[5], label = paste0("Q3: ", Qs[3], "%") , size= 5) + 
    annotate("text", x = quantile(xMinMax, probs = seq(0, 1, 0.01))[90] , y = quantile(yMinMax, probs = seq(0, 1, 0.01))[5], label = paste0("Q4: ", Qs[4], "%") , size= 5) + 
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5), axis.text = element_text(size = 10), legend.position = "none")
  
  return(p)
}


