## gene expression with average value and its porportion
CellGeneConditionPert <- function(exp_data, mydata, myGroup, myType, GroupLevel, TypeLevel, MainTitle = "", style){
  library(ggplot2)
  library(tidyverse)
  
  if (style == "multiFactor"){ # there are two factors in x axis
    PertRe <- c(); PertMean <- c(); PertG <- c(); PertT <- c(); Gene <- c()
    table(mydata[,myGroup], mydata[,myType]); table(mydata[,myGroup]); table(mydata[,myType])
    for (i in GroupLevel){ # i <- GroupLevel[2]
      for (j in TypeLevel){ # j <- TypeLevel[1]
        temp_data <- as.data.frame(exp_data[,(mydata[,myGroup] %in% i) & (mydata[,myType] %in% j)])
        rownames(temp_data) <- rownames(exp_data); colnames(temp_data) <- colnames(exp_data)[(mydata[,myGroup] %in% i) & (mydata[,myType] %in% j)]
        if (dim(temp_data)[2] != 0){
          temp_Pert <- round(rowSums(temp_data > 0) / ncol(temp_data), 4) * 100
          temp_mean <- Matrix::rowMeans(temp_data) 
          Gene <- c(Gene, rownames(temp_data)); PertRe <- c(PertRe, temp_Pert); PertMean <- c(PertMean, temp_mean)
          PertG <- c(PertG, rep(i, nrow(temp_data))); PertT <- c(PertT, rep(j, nrow(temp_data)))
        } else {
          Gene <- c(Gene, rownames(temp_data)); PertRe <- c(PertRe, rep(0, nrow(temp_data))); PertMean <- c(PertMean, rep(0, nrow(temp_data)))
          PertG <- c(PertG, rep(i, nrow(temp_data))); PertT <- c(PertT, rep(j, nrow(temp_data)))
        }
      }
    }
    
    df <- data.frame(Gene, PertRe, PertMean, PertG, PertT)
    df$Gene <- factor(df$Gene, levels = rownames(exp_data)); df$PertT <- factor(df$PertT, levels = TypeLevel)
    df$TypeGroup <- paste0(df$PertT, "_", df$PertG); df$TypeGroup <- factor(df$TypeGroup, levels = lapply(TypeLevel, function(x) paste(x,combn(GroupLevel,1),sep = "_")) %>% unlist)
    Nz <- length(TypeLevel); Ng <- length(GroupLevel)
    p <- ggplot(df, aes(x=TypeGroup, y = Gene, color = PertMean, size = PertRe)) + geom_point() + labs(x="",y="", title= MainTitle) +
      scale_color_gradientn(colours = viridis::cividis(20), limits = c(0, ceiling(max(df$PertMean))), name = "ExpLevel") + 
      scale_size(name = "ExpPortion(%)", limits = c(0, 100)) + 
      geom_vline(xintercept = (1:(Nz-1)) * Ng + 0.5, linetype="dashed", col="black", size=0.5) + 
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 25), legend.position = "right")
    
  }
  
  if (style == "sigFactor"){ # single facotr on x axis，based on myType
    PertRe <- c(); PertMean <- c(); PertT <- c(); Gene <- c()
    for (j in TypeLevel){ # j <- TypeLevel[1]
      temp_data <- as.data.frame(exp_data[,(mydata[,myType] %in% j)])
      rownames(temp_data) <- rownames(exp_data); colnames(temp_data) <- colnames(exp_data)[(mydata[,myType] %in% j)]
      if (dim(temp_data)[2] != 0){
        temp_Pert <- round(rowSums(temp_data > 0) / ncol(temp_data), 4) * 100
        temp_mean <- Matrix::rowMeans(temp_data) 
        Gene <- c(Gene, rownames(temp_data)); PertRe <- c(PertRe, temp_Pert); PertMean <- c(PertMean, temp_mean)
        PertT <- c(PertT, rep(j, nrow(temp_data)))
      } else {
        Gene <- c(Gene, rownames(temp_data)); PertRe <- c(PertRe, rep(0, nrow(temp_data))); PertMean <- c(PertMean, rep(0, nrow(temp_data)))
        PertT <- c(PertT, rep(j, nrow(temp_data)))
      }
    }
    
    df <- data.frame(Gene, PertRe, PertMean, PertT); df$PertT <- factor(df$PertT, levels = TypeLevel); df$Gene <- factor(df$Gene, levels = rownames(exp_data))
    p <- ggplot(df, aes(x=PertT, y = Gene, color = PertMean, size = PertRe)) + geom_point() + labs(x="",y="", title= MainTitle) +
      scale_color_gradientn(colours = viridis::cividis(20), limits = c(0, ceiling(max(df$PertMean))), name = "ExpLevel") + 
      scale_size(name = "ExpPortion(%)", limits = c(0, 100)) + 
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 25), legend.position = "right")
  }
  
  return(p)
}
                                                                                                          

                                                                                                          
## stacked violin plots
# exp_data：Gene-Cell expression matrix(行列有标注)
# mydata是对细胞列的标注，是数据框，包括CellId(与exp_data列顺序一致且相同)
CellGeneStackViolin <- function(exp_data, mydata, myType, TypeLevel, TypeColor, MainTitle = "", style = "1"){
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(tidyverse)
  library(reshape2)
  
  exp_data$gene <- rownames(exp_data)
  df <- melt(exp_data, id="gene")
  colnames(df)[c(2,3)] <- c("ID","exp"); colnames(mydata)[1] <- "ID"
  df <- inner_join(df, mydata , by = "ID")
  df$gene <- factor(df$gene, levels =  rownames(exp_data)); df[,myType] <- factor(df[,myType] , levels =  TypeLevel)
  if (is.numeric(TypeColor)) { TypeColor <- ggsci::pal_igv("default")(51)[1:TypeColor]}
  
  if (style == 1){
    p <- ggplot(df, aes(x = df[,myType], y = exp)) + geom_violin(aes(fill=df[,myType])) + scale_fill_manual(values = TypeColor) + facet_grid(df$gene~.,) + labs(x="" ,y="", title= MainTitle) +
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),  
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y =  element_blank(), plot.title = element_text(hjust = 0.5, size = 25), legend.position = "none")
  } else{
    p <- ggplot(df, aes(x = exp, y = gene)) + geom_violin(aes(fill=df[,myType])) + scale_fill_manual(values = TypeColor) + facet_grid(.~df[,myType]) + labs(x="" ,y="", title= MainTitle) +
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), 
            axis.text.y = element_text(size = 15), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 25), legend.position = "none")
  }
  
  return(p)
}

                                                                                                          
## Block with differential color to show average gene expression
# exp_data：Gene-Cell expression matrix(行列有标注)
# mydata是对细胞列的标注，是数据框，包括CellId(与exp_data列顺序一致且相同)
CellGeneMeanBlock <- function(exp_data, mydata, myType, TypeLevel, scaleOpt = "no",  MainTitle = ""){
  library(ggplot2)
  library(tidyverse)
  library(reshape2)
  
  MeanRe <- data.frame()
  for (item in TypeLevel){ # item <- TypeLevel[1]
    temp_data <- exp_data[,which(mydata[,myType] %in% item)]
    temp_means <- Matrix::rowMeans(temp_data) 
    MeanRe <- rbind(MeanRe, temp_means)
  }; rownames(MeanRe) <- TypeLevel; colnames(MeanRe) <- rownames(exp_data)
  
  ifelse((scaleOpt == "yes"), df <- data.frame(scale(MeanRe)) , df <- data.frame(MeanRe) )
  df$Type <- rownames(df)
  df <- melt(df, id = "Type"); colnames(df)[c(2,3)] <- c("Gene","exp")
  df$Type <- factor(df$Type, levels = TypeLevel); df$Gene <- factor(df$Gene, levels = rownames(exp_data))
  
  p <- ggplot(df, aes(x = Type, y = Gene, fill = exp)) + geom_tile(colour = "white", size = 0.25) + 
    scale_fill_gradient2(low = "navy", mid = "white", high = "Firebrick3", name = "GeneExp") +
    scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) + labs(x="", y="", title= MainTitle) + # 移除多餘空白
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(), 
          axis.text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
          plot.title = element_text(hjust = 0.5, size = 25), legend.position = "right")
  
  return(p)
}                                                                                                        
