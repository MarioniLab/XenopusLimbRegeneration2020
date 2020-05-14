#Load reticulate and specify python directory
library("reticulate")
use_python("/Users/hiscoc01/miniconda3/bin/python3", required = T)
py_config()

#Load other packages
library(gridExtra)
library(gatepoints)
library(ggplot2)
library(Matrix)
library(igraph)
library(AUCell)
library(edgeR)
library(ggrepel)
library(reshape2)
library(Seurat)
library(scran)
library(DropletUtils)
library(gtools)
library(NNLM)
library(hexbin)
library(infotheo)
library(ggpubr)

#Define python variables and functions
umap <- import("umap")
np <- import("numpy", convert = FALSE)
scp <- import("scipy")
py_run_string("
import umap
import scipy as scp
import sklearn as sk
from sklearn.utils import check_random_state
from sklearn.utils import check_array
from scipy.sparse import csr_matrix, find
import numpy as np
              ")

py_run_string("
import umap
def get_umap(y, n, minDist,seed):
    return umap.UMAP(n_neighbors=n,min_dist=minDist,metric='cosine', n_components = 2, random_state = seed).fit_transform(y)
              ")

py_run_string("
import umap
def get_umap_euclidean(y, n, minDist,seed):
  return umap.UMAP(n_neighbors=n,min_dist=minDist,metric='euclidean', n_components = 2, random_state = seed).fit_transform(y)
              ")

py_run_string("
def get_umap_G(y,n, seed):
  random_state = check_random_state(seed)
  Y = check_array(y).astype(np.float64)
  G = umap.umap_.fuzzy_simplicial_set(X=Y,n_neighbors=n,random_state=random_state,metric='cosine')
  return find(G)
              ")
assembleCounts <- function(sce){
  counts <- counts(sce)
  rownames(counts) <- rowData(sce)$Symbol
  colnames(counts) <- sce$Barcode
  rownames(counts) <- make.names(rownames(counts), unique = TRUE)
  colnames(counts) <- make.names(colnames(counts), unique = TRUE)
  return(counts)
}

assembleMeta <- function(sce, labels){
  meta <- data.frame(cell = sce$Barcode)
  meta$cell <-  make.names(meta$cell, unique = TRUE)
  meta$condition <- labels$Condition[match(sce$Sample, labels$Sample)]
  meta$stage <- labels$Stage[match(sce$Sample, labels$Sample)]
  meta$lane <- labels$Lane[match(sce$Sample, labels$Sample)]
  meta$position <- labels$Position[match(sce$Sample, labels$Sample)]
  meta$batch <- labels$Batch..lib.prep.[match(sce$Sample, labels$Sample)]
  meta$sample <- sce$Sample
  meta$val <- paste0(meta$stage,meta$condition,"_", meta$position)
  return(meta)
}

#Normalize by total counts
normalize <- function(counts, pseudocount = 1e4){
  librarySize <- Matrix::colSums(counts)
  return(pseudocount * t(t(counts) / librarySize ))
}

#Select HVGs based on mean expression and fano factor
compute_hvg <- function(counts.normalized, lowmean = 0.05, highmean = 0.8, fano = 0.65){
  mean.expression <- Matrix::rowMeans(counts.normalized)
  meansq.expression <- Matrix::rowMeans(counts.normalized^2)
  fano.factor <- (meansq.expression - (mean.expression^2))/mean.expression
  hvg <- counts.normalized[which((mean.expression >= quantile(mean.expression, lowmean)) & (fano.factor >= quantile(na.omit(fano.factor),fano)) & (mean.expression <= quantile(mean.expression,highmean))),]
  return(hvg)
}

#Project data using UMAP
umap_project <- function(meta, hvg, neighbours = 10, dist = 0.3, seed = 42){
  y <- np$log2(as.matrix(hvg) + 1)$T
  out <- py$get_umap(y,as.integer(neighbours),dist, as.integer(seed))
  meta$x <- out[,1]
  meta$y <- out[,2]
  return(meta)
}

umap_project_custom <- function(meta, hvg, neighbours = 10, dist = 0.3, seed = 42){
  y <- np$abs(as.matrix(hvg))$T
  out <- py$get_umap(y,as.integer(neighbours),dist, as.integer(seed))
  meta$x <- out[,1]
  meta$y <- out[,2]
  return(meta)
}

#Cluster data using UMAP
umap_cluster <- function(meta, hvg, neighbours = 10, steps = 4, seed1 = 42, seed2 = 42){
  y <- np$log2(as.matrix(hvg) + 1)$T
  x <- py$get_umap_G(y,as.integer(neighbours),as.integer(seed1))
  sp <- sparseMatrix(i = (x[[1]]+1), j = (x[[2]]+1), x = as.vector(x[[3]]), dims = c(dim(hvg)[2],dim(hvg)[2]))
  g <- graph_from_adjacency_matrix(sp,weighted = TRUE,mode = "undirected")
  set.seed(seed2)
  meta$cluster <- cluster_walktrap(g,steps = steps)$membership
  return(meta)
}

#Cluster data using UMAP
umap_g <- function(meta, hvg, neighbours = 10, steps = 4, seed1 = 42){
  y <- np$log2(as.matrix(hvg) + 1)$T
  x <- py$get_umap_G(y,as.integer(neighbours),as.integer(seed1))
  sp <- sparseMatrix(i = (x[[1]]+1), j = (x[[2]]+1), x = as.vector(x[[3]]), dims = c(dim(hvg)[2],dim(hvg)[2]))
  g <- graph_from_adjacency_matrix(sp,weighted = TRUE,mode = "undirected")
  # set.seed(seed2)
  # meta$cluster <- cluster_walktrap(g,steps = steps)$membership
  return(g)
}

#Cluster data using UMAP, louvain clustering
umap_cluster_louvain <- function(meta, hvg, neighbours = 10, seed1 = 42, seed2 = 42){
  y <- np$log2(as.matrix(hvg) + 1)$T
  x <- py$get_umap_G(y,as.integer(neighbours),as.integer(seed1))
  sp <- sparseMatrix(i = (x[[1]]+1), j = (x[[2]]+1), x = as.vector(x[[3]]), dims = c(dim(hvg)[2],dim(hvg)[2]))
  g <- graph_from_adjacency_matrix(sp,weighted = TRUE,mode = "undirected")
  set.seed(seed2)
  meta$cluster <- cluster_louvain(g)$membership
  return(meta)
}

#Plot UMAP projection in 2D - with or without clusters colored differently
plotMeta <- function(meta, mode = "none"){
  if(mode == "none"){
    p <- ggplot( data = meta, mapping= aes(x = x, y = y)) +
      theme_bw() +
      geom_point(size=0.1, col = "black") +
      coord_equal() +
      theme(aspect.ratio = 1)+
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1) 
  }
  
  else if(mode == "cluster"){
    centroid <- aggregate(cbind(x,y) ~ cluster, data=meta, FUN=median)
    p <- ggplot( data = meta, mapping= aes(x = x, y = y, col = factor(meta$cluster))) +
      theme_bw() +
      geom_point(size=0.1) +
      coord_equal() +
      theme(aspect.ratio = 1)+
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1)  +
      theme(legend.position = "none") +
        geom_text(data = centroid, mapping=aes(x = x, y = y, label = cluster), col = "black")
      
  }
  else if(mode == "cycling"){
    p <- ggplot( data = meta, mapping= aes(x = x, y = y, col = factor(meta$CellCyclePhase))) +
      theme_bw() +
      geom_point(size=0.1) +
      coord_equal() +
      theme(aspect.ratio = 1)+
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1)  


  }
  
  return(p)
}

#Plot UMAP projection in 2D, with cluster names and colours specified by "clusters"
plotAnnotation <- function(meta, clusters){
   
  p <-  ggplot( data = meta, mapping= aes(x = x, y = y, col = factor(meta$cluster, levels = clusters$names))) +
      theme_bw() +
      geom_point(size=0.5) +
      scale_colour_manual(values = sapply(clusters$cols, as.character), labels= sapply(clusters$names, as.character), name = "") +
      coord_equal() +
      theme(aspect.ratio = 1)+
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1)  +
      guides(colour = guide_legend(override.aes = list(size=3,
                                                       alpha = 1))) + guides(fill=guide_legend(nrow=20, byrow = TRUE), label.theme = element_text(size = 1))
  
  print(p)
  return(p)
  
}


#Plot a boxplot of values (e.g. gene expression levels) separated by groups (e.g. cell types)
scBox <- function(values, cells = 1:length(values),groups, colours, title = "none", no_xaxis = FALSE, width = 1, xlabel = "none", ylabel = "none", return = FALSE, aspect.ratio = 1){
  m <- intersect(cells,which(groups %in% colours$names))
  cluster <- factor(groups[m], levels = sapply(colours$names, as.character), ordered = TRUE)
  df <- data.frame(score = as.vector(values)[m],cluster = cluster)
  p <- ggplot(df,aes(x = cluster, y=score, fill = cluster))+
    geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=0.4, notch=FALSE, width = width)+
    theme_classic() + theme(aspect.ratio = aspect.ratio)+
    scale_fill_manual(values = sapply(colours$cols, as.character), labels= sapply(colours$names, as.character), name = "") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), axis.line.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    if(no_xaxis) {p <- p + theme(axis.text.x = element_blank())}
  if(title != "none") { p <- p + ggtitle(title)}  
  if(xlabel != "none"){ p <- p + xlab(xlabel)}
  if(ylabel != "none"){ p <- p + ylab(ylabel)}
  p <- p + scale_y_continuous(expand = expand_scale(mult = c(0, 0)))
  
  if(return) {return(p)}
  else {print(p)}
}

#Plot a boxplot of values (e.g. gene expression levels) separated by groups (e.g. cell types)
scViolin <- function(values, cells = 1:length(values),groups, colours, title = "title", width = 1, aspect.ratio = 1, xlabel = "none", ylabel = "none"){
  m <- intersect(cells,which(groups %in% colours$names))
  cluster <- factor(groups[m], levels = sapply(colours$names, as.character), ordered = TRUE)
  df <- data.frame(score = as.vector(values)[m],cluster = cluster)
  p <- ggplot(df,aes(x = cluster, y=score, fill = cluster))+
    geom_violin(width = width)+
    theme_classic() +
    scale_fill_manual(values = sapply(colours$cols, as.character), labels= sapply(colours$names, as.character), name = "") +
    theme(aspect.ratio = aspect.ratio) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), axis.line.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(title)
  if(xlabel != "none"){ p <- p + xlab(xlabel)}
  if(ylabel != "none"){ p <- p + ylab(ylabel)}
  p <- p + scale_y_continuous(expand = c(0, 0))
  print(p)
}

#Plot cell cycle fractions for different groups (e.g. clusters)
plotCellCycle <- function(meta, clusters, title = "title"){
  output <- c()
  for (category in sort(unique(meta$CellCyclePhase))){
    output <- rbind(output,aggregate(meta$CellCyclePhase == category,list(meta$cluster),sum)$x)
  }
  rownames(output) <- sort(unique(meta$CellCyclePhase))
  colnames(output) <- aggregate((meta$CellCyclePhase == category),list(meta$cluster),sum)$Group.1
  output <- output[,match(clusters$names,colnames(output))]
  
  output <- output[c(1,3,2),]
  
  par(mfrow=c(1, 1), mar=c(10, 5, 4, 8)); barplot(prop.table(output,2),las=2,legend.text = TRUE,args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)), cex.names = 0.65, col = c("grey","#377eb8","#e41a1c")); title(title)
}

#Calculate the fraction of cycling cells for different conditions, for each cluster. We ignore clusters with low abundances for which this estimate will be highly noisy.
getCycling <- function(meta, clusters){
  output <- c()
  for (category in sort(unique(meta$CellCyclePhase))){
    output <- rbind(output,aggregate(meta$CellCyclePhase == category,list(meta$cluster),sum)$x)
  }
  rownames(output) <- sort(unique(meta$CellCyclePhase))
  colnames(output) <- aggregate((meta$CellCyclePhase == category),list(meta$cluster),sum)$Group.1
  output <- output[,match(clusters$names,colnames(output))]
  output <- output[c(1,3,2),]
  numbers <- apply(output,2,sum)
  proliferating <- 1-prop.table(output,2)[1,]
  proliferating[numbers < 10] <- NA
  return(proliferating)
}

#Overlay different samples on the UMAP projection
plotSample <- function(meta,m1, title = "none", mode = "none", size = 0.3, lims = "none", return = FALSE){
  umap <- as.data.frame(cbind(meta$x,meta$y))
  umap1 <- as.data.frame(cbind(meta$x[m1],meta$y[m1]))
  umap1$V3 <- meta$sample[m1]
  if(mode == "none"){ cols <- "red"}
  else if(mode == "batch"){ cols <- factor(umap1$V3)}
  p <- ggplot( data = umap, mapping= aes(x = V1, y = V2)) +
    theme_bw() +
    geom_point(size=size, col = "grey", stroke = 0) +
    scale_shape(solid = F) +
    coord_equal() +
    geom_point(data = umap1, mapping = aes(x = as.numeric(V1), y = as.numeric(V2), col = cols),size=size, fill = NA, stroke = 0) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    guides(colour = guide_legend(override.aes = list(size=3,
                                                     alpha = 1))) + theme(legend.position = "none")
  if(title != "none"){
    p <- p + ggtitle(title)
  }
  if(lims != "none"){
    p <- p + xlim(lims[1:2]) + ylim(lims[3:4])
  }
  if(return == TRUE){
    return(p)
  }
  else{
    print(p)
  }
}

plotSampleGrid <- function(meta,value,samples,size = 0.1,lims = "none", mode="none", name = "name"){
  samplePlots <- list()
  for(val in amputation){
    p <- plotSample(meta, which(meta$val == val), size = size, lims = lims, return = TRUE, mode = mode)
    samplePlots[[val]] <- p + theme(plot.margin=unit(rep(0.2,4), "cm")) + theme(panel.border = element_blank())
  }
  grid.arrange(grobs = samplePlots,widths = c(1,1,1), layout_matrix = rbind(c(1,3,NA),c(2,4,6)))
  g <- arrangeGrob(grobs = samplePlots,widths = c(1,1,1), layout_matrix = rbind(c(1,3,NA),c(2,4,6)))

}
    



#Plot gene expression level for a subset of cells, overlaid on the entire map
plotGeneSubset <- function(meta, countn, geneName,condition, title){
  XA <- meta$x
  YA <- meta$y
  X <- meta$x[condition]
  Y <- meta$y[condition]
  atlas <- as.data.frame(cbind(XA, YA))
  all <- as.data.frame(cbind(X, Y))
  gene <- countn[which(rownames(countn) == geneName), condition]
  if (length(gene) > 0) {
    log_counts = log10(gene + 1)
    Xsub <- X[log_counts > 0]
    Ysub <- Y[log_counts > 0]
    log_counts = log_counts[log_counts > 0]
    m <- order(log_counts)
    Xsub <- Xsub[m]
    Ysub <- Ysub[m]
    log_counts <- log_counts[m]
    if (length(Xsub) > 0) {
      subset <- as.data.frame(cbind(Xsub, Ysub))
      p <- ggplot(data = atlas, mapping = aes(x = XA, y = YA)) +
        geom_point(size = 0.1, col = alpha("grey", 1.0)) +
        geom_point(
          data = all,
          mapping = aes(x = X, y = Y),
          size = 1,
          shape = 21,
          colour = "black",
          stroke = 0.1,
          fill = "grey"
        ) +
        geom_point(
          data = subset,
          mapping = aes(x = Xsub, y = Ysub, fill = log_counts),
          size = 1,
          shape = 21,
          colour = "black",
          stroke = 0.1
        ) +
        scale_fill_gradient(low = "yellow",
                            high = "red",
                            limits = c(0, 1)) +
        theme_bw() +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1
        ) +
        ggtitle(paste0(geneName,":   ",title))
      print(p)
      
    }
  }
}




#Compute differential abundance statistics between conditions
differentialAbundance <- function(meta,labels,sample1,sample2,clusters, plot = T,minN = 5, pthresh = 0.05){
  clusterColours <- sapply(clusters$cols,as.character)
  clusterNames <- sapply(clusters$names,as.character)
  counts <- table(meta$cluster, meta$sample)
  counts <- counts[which(rowMeans(counts) > minN),]
  group <- factor(labels$val[match(colnames(counts),labels$Sample)])
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)

  y <- DGEList(counts=counts,group=group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)

  m1 <- match(sample1, colnames(design))
  m2 <- match(sample2, colnames(design))
  
  contrast <- rep(0,dim(design)[2])
  contrast[m2] <- 1
  contrast[m1] <- -1
  lrt <- glmQLFTest(fit, contrast = contrast)$table
  
  subset <- lrt[which(lrt$PValue < pthresh  ),]
  if (plot){
    p <- ggplot(lrt, aes(x=logFC, y=-log10(PValue), col = factor(rownames(lrt), levels = rownames(lrt)))) + 
      theme_bw() +
      geom_point(size=2) +
      geom_text_repel(data = subset, mapping= aes(x = logFC, y = -log10(PValue), label=rownames(subset),col = factor(rownames(subset), levels = rownames(subset))), size=3) +
      
      scale_colour_manual(values = clusterColours[match(rownames(lrt),(clusterNames))], labels= rownames(lrt), name = "") +
      theme(  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1, legend.position = "none")  +
      theme(aspect.ratio = 1) +
      guides(colour = guide_legend(override.aes = list(size=3,
                                                       alpha = 1))) +
      labs(x = "log fold change in abundance", y = "-log P-value") +
      ggtitle(paste0(sample2,"  vs.  ", sample1))
    print(p)
  }
  return(lrt)
}



#Plot gene expression values on the UMAP projection
plotGene <- function(meta, countn, geneName, lims = "none", size = 0.2){
  
  X <- meta$x
  Y <- meta$y
  atlas <- as.data.frame(cbind(X,Y))
  if(!(geneName %in% rownames(countn))){
    return()
  }
  gene <- countn[which(rownames(countn) == geneName),]
  if(max(gene) == 0){ p <- ggplot( data = atlas, mapping=aes(x=X,y=Y)) +
    geom_point(size=size,col=alpha("grey",1)) +
    theme_bw() +
    # scale_colour_gradient2(low="#ff6666", mid = "#e60000", high = "#990000") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", aspect.ratio = 1) +
    ggtitle(geneName)
  return(p)}
  log_counts = log10(gene+1)
  Xsub <- X[log_counts>0]
  Ysub <- Y[log_counts>0]
  log_counts = log_counts[log_counts>0]
  m <- order(log_counts)
  Xsub <- Xsub[m]
  Ysub <- Ysub[m]
  log_counts <- log_counts[m]
  subset <- as.data.frame(cbind(Xsub, Ysub))
  p <- ggplot( data = atlas, mapping=aes(x=X,y=Y)) +
    geom_point(size=size,col=alpha("grey",1)) +
    geom_point(data = subset,mapping = aes(x=Xsub, y=Ysub, col = log_counts), size = size) +
    theme_bw() +
    # scale_colour_gradient2(low="#ff6666", mid = "#e60000", high = "#990000") + 
    scale_colour_gradient(low="yellow", high = "red") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    coord_equal() + 
    
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  ggtitle(geneName)
  
  if(lims != "none"){
    p <- p + xlim(lims[1:2]) + ylim(lims[3:4])
  }

  return(p)
  
}


differentialAbundanceHex <- function(meta,labels,sample1,sample2,clusters, nbins, mode = "pvalue", clims){
  drhex <- hexbin(meta$x,meta$y, xbins = nbins, shape = 1, IDs = TRUE)
  cID <- drhex@cID
  drhex <- data.frame(x = as.numeric(hcell2xy(drhex)$x), y= as.numeric(hcell2xy(drhex)$y))
  tmp <- meta
  tmp$cluster <- cID
  da <- differentialAbundance(tmp,labels,sample1,sample2, minN = 0, pthresh = 1, clusters = clusters, plot = F)
  
  if(mode == "pvalue"){
  drhex$plot <- -log10(da$PValue)*sign(da$logFC)
  }
  else{
    drhex$plot <- da$logFC

  }
    
  
  
  
 p <-  ggplot(drhex, mapping = aes(x, y, fill = plot, colour = plot)) +
    geom_hex(stat = "identity") +
    scale_fill_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0,limits = clims, oob = scales::squish) +
    scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0,limits = clims, oob = scales::squish) +
  
   # scale_fill_gradient2(low = "blue", mid = "gray", high = "red") +
   # scale_color_gradient2(low = "blue", mid = "gray", high = "red") +

    
    theme_classic() + theme(legend.position = "none") +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1, axis.line = element_blank())  +
   theme(plot.margin=unit(rep(-0.2,4), "cm"))
 return(p)
}


#Plot arbitrary vector on the UMAP projection
plotArbitrary <- function(meta, input, title = "title", lims = "none", size = 0.2, return = FALSE){
  
  X <- meta$x
  Y <- meta$y
  atlas <- as.data.frame(cbind(X,Y))
  

  log_counts = input
  Xsub <- X[log_counts>0]
  Ysub <- Y[log_counts>0]
  log_counts = log_counts[log_counts>0]
  m <- order(log_counts)
  Xsub <- Xsub[m]
  Ysub <- Ysub[m]
  log_counts <- log_counts[m]
  subset <- as.data.frame(cbind(Xsub, Ysub))
  p <- ggplot( data = atlas, mapping=aes(x=X,y=Y)) +
    geom_point(size=size,col=alpha("grey",1)) +
    geom_point(data = subset,mapping = aes(x=Xsub, y=Ysub, col = log_counts), size = size) +
    theme_bw() + coord_equal() +
    # scale_colour_gradient2(low="#ff6666", mid = "#e60000", high = "#990000") + 
    scale_colour_gradient(low="yellow", high = "red") + 
    
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
    ggtitle(title)
  if(lims != "none"){
    p <- p + xlim(lims[1:2]) + ylim(lims[3:4])
  }
  if(return) {return(p)}
  else{print(p)}
  
  
}


#Plot UMI on the UMAP projection
plotUMI <- function(meta){
  
  X <- meta$x
  Y <- meta$y
  atlas <- as.data.frame(cbind(X,Y))
  log_counts = log10(meta$umi+1)
  Xsub <- X[log_counts>0]
  Ysub <- Y[log_counts>0]
  log_counts = log_counts[log_counts>0]
  m <- order(log_counts)
  Xsub <- Xsub[m]
  Ysub <- Ysub[m]
  log_counts <- log_counts[m]
  subset <- as.data.frame(cbind(Xsub, Ysub))
  p <- ggplot( data = atlas, mapping=aes(x=X,y=Y)) +
    geom_point(size=0.2,col=alpha("grey",1)) +
    geom_point(data = subset,mapping = aes(x=Xsub, y=Ysub, col = log_counts), size = 0.2) +
    theme_bw() +
    scale_colour_gradient(low="yellow", high = "red") + 
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", aspect.ratio = 1) +
    ggtitle("UMI")
  print(p)
  return()
  
}

#Given a set of gene names, create .L and .S versions (long and short chromosomes, since X. Laevis is pseudotetraploid)
xenifyGenes <- function(geneList){
  
  ligands.L<-sapply(tolower((geneList)),function(x) paste0(x,".L"))
  ligands.S<-sapply(tolower(geneList),function(x) paste0(x,".S"))
  return((c(ligands.L,ligands.S)))

}

#Construct a heatmap for the average expression of different genes across different conditions/clusters
heatmap <- function(countn, genes = rownames(countn), cells = 1:length(colnames(countn)), condition, norm = "max", conditionList = "none", aspect.ratio = 1, xAngle = 45, yAngle = 0, xSize = 10, ySize = 10, xAdj = 1, plot = TRUE, collow = "white", colhigh = "steelblue", title = "title", return = FALSE, normcells = 1:length(colnames(countn))){
  
  genes <- rev(genes)
  
  
  subset <- log10(countn[na.omit(match(genes,rownames(countn))),cells]+1)
  avgd <- (t(aggregate(t(as.matrix(subset)),list(condition[cells]),mean)))
  clusters <- avgd[1,]
  avgd <- avgd[-1,]
  avgd <- matrix(as.numeric(unlist(avgd)),nrow=nrow(avgd))
  rownames(avgd) <- rownames(subset)
  colnames(avgd) <- clusters
  
  
  
  if (norm == "max"){
    avgd_n <- (avgd / apply(avgd,1,function(x) max(x)))
  } else if (norm == "sum"){
    avgd_n <- avgd / rowSums(avgd)
  } 
  else if (norm == "totalMax"){
    subsetMax <- log10(countn[na.omit(match(genes,rownames(countn))),normcells]+1)
    avgdMax <- (t(aggregate(t(as.matrix(subsetMax)),list(condition[normcells]),mean)))
    clustersMax <- avgdMax[1,]
    avgdMax <- avgdMax[-1,]
    avgdMax <- matrix(as.numeric(unlist(avgdMax)),nrow=nrow(avgdMax))
    rownames(avgdMax) <- rownames(subsetMax)
    colnames(avgdMax) <- clustersMax
    avgd_n <- (avgd / apply(avgdMax,1,function(x) max(x)))
    
  } else {    avgd_n <- avgd
  }
  
  if(conditionList != "none"){
    avgd_n <- avgd_n[,match(conditionList,colnames(avgd_n))]
  }

  if(length(which(is.na(rowSums(avgd_n)))) > 0){
    avgd_n <- avgd_n[-which(is.na(rowSums(avgd_n))),]
  }
  
  melted_ <- melt(t(avgd_n))
  
  p <- ggplot(data = melted_, aes(x=Var1, y=Var2, fill=value)) + 
    theme_bw() +
    # geom_tile() + 
    geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = collow,high = colhigh)+
    theme(axis.text.x = element_text(angle = xAngle, hjust = xAdj, size = xSize), axis.text.y = element_text(angle = yAngle, size = ySize), aspect.ratio = aspect.ratio, axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(panel.background = element_blank()) + theme(panel.grid = element_blank())
  p <- p + ggtitle(title)
  if(plot == TRUE){
    print(p)
  }
  if (return) {return(avgd_n)}
}


#Construct a heatmap for the average expression of different genes across different conditions/clusters
dotplot <- function(countn, genes = rownames(countn), cells = 1:length(colnames(countn)), condition, norm = "max", conditionList = "none", aspect.ratio = 1, xAngle = 45, yAngle = 0, xSize = 10, ySize = 10, xAdj = 1, plot = TRUE, 
                    title = "title", return = FALSE, normcells = 1:length(colnames(countn)), dot.scale = 2, collow = "lightgrey", colhigh = "blue"){
  

  genes <- rev(genes)
  
  
  subset <- log10(countn[na.omit(match(genes,rownames(countn))),cells]+1)
  nonzero <- countn[na.omit(match(genes,rownames(countn))),cells] > 0
  
  avgd <- (t(aggregate(t(as.matrix(subset)),list(condition[cells]),mean)))
  clusters <- avgd[1,]
  avgd <- avgd[-1,]
  avgd <- matrix(as.numeric(unlist(avgd)),nrow=nrow(avgd))
  rownames(avgd) <- rownames(subset)
  colnames(avgd) <- clusters
  
  pct <- (t(aggregate(t(as.matrix(nonzero)),list(condition[cells]),mean)))
  clusters <- pct[1,]
  pct <- pct[-1,]
  pct <- matrix(as.numeric(unlist(pct)),nrow=nrow(pct))
  rownames(pct) <- rownames(nonzero)
  colnames(pct) <- clusters
  
  
  
  if (norm == "max"){
    avgd_n <- (avgd / apply(avgd,1,function(x) max(x)))
  } else if (norm == "sum"){
    avgd_n <- avgd / rowSums(avgd)
  } 
  else if (norm == "totalMax"){
    subsetMax <- log10(countn[na.omit(match(genes,rownames(countn))),normcells]+1)
    avgdMax <- (t(aggregate(t(as.matrix(subsetMax)),list(condition[normcells]),mean)))
    clustersMax <- avgdMax[1,]
    avgdMax <- avgdMax[-1,]
    avgdMax <- matrix(as.numeric(unlist(avgdMax)),nrow=nrow(avgdMax))
    rownames(avgdMax) <- rownames(subsetMax)
    colnames(avgdMax) <- clustersMax
    avgd_n <- (avgd / apply(avgdMax,1,function(x) max(x)))
    
  } else {    avgd_n <- avgd
  }
  
  if(conditionList != "none"){
    avgd_n <- avgd_n[,match(conditionList,colnames(avgd_n))]
    pct <- pct[,match(conditionList,colnames(pct))]
    
  }
  
  if(length(which(is.na(rowSums(avgd_n)))) > 0){
    pct <- pct[-which(is.na(rowSums(avgd_n))),]
    avgd_n <- avgd_n[-which(is.na(rowSums(avgd_n))),]
  }
  
  melted_ <- melt(t(avgd_n))
  pct_ <- melt(t(pct))
  melted_$pct <- 100*pct_$value
  
  p <- ggplot(data = melted_) + 
    theme_bw() +
    geom_point(mapping = aes(x = Var1, y = Var2,size = pct,col = value))+
    scale_color_gradient(low = collow,high = colhigh)+
    scale_size(range = c(0,dot.scale), limits = c(0, 100)) +
    theme(axis.text.x = element_text(angle = xAngle, hjust = xAdj, size = xSize), axis.text.y = element_text(angle = yAngle, size = ySize), aspect.ratio = aspect.ratio, axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(panel.background = element_blank()) + theme(panel.grid = element_blank())
  
  p <- p + ggtitle(title)
  if(plot == TRUE){
    print(p)
  }
  if (return) {return(avgd_n)}
}






heatmapGeneric <- function(X, features = rownames(X), cells = 1:length(colnames(X)), condition, norm = "max", conditionList = "none", aspect.ratio = 1, xAngle = 45, yAngle = 0, xSize = 10, ySize = 10, xAdj = 1, plot = TRUE, title = "title"){
    
  
  
  subset <- (X[na.omit(match(features,rownames(X))),cells]+1)
  avgd <- (t(aggregate(t(as.matrix(subset)),list(condition[cells]),mean)))
  clusters <- avgd[1,]
  avgd <- avgd[-1,]
  avgd <- matrix(as.numeric(unlist(avgd)),nrow=nrow(avgd))
  rownames(avgd) <- rownames(subset)
  colnames(avgd) <- clusters
  
  
  
  if (norm == "max"){
    avgd_n <- (avgd / apply(avgd,1,function(x) max(x)))
  } else if (norm == "sum"){
    avgd_n <- avgd / rowSums(avgd)
  } else if (norm == "minmax"){
    avgd_n <- (avgd - apply(avgd,1,function(x) min(x)))
    avgd_n <- (avgd_n / apply(avgd_n,1,function(x) max(x)))
  } else if (norm == "total"){
  } else {    avgd_n <- avgd
  }
  
  if(conditionList != "none"){
    avgd_n <- avgd_n[,match(conditionList,colnames(avgd_n))]
  }
  
  if(length(which(is.na(rowSums(avgd_n)))) > 0){
    avgd_n <- avgd_n[-which(is.na(rowSums(avgd_n))),]
  }
  
  melted_ <- melt(t(avgd_n))
  
  p <- ggplot(data = melted_, aes(x=Var1, y=Var2, fill=value)) + 
    theme_bw() +
    # geom_tile() + 
    geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")+
    theme(axis.text.x = element_text(angle = xAngle, hjust = xAdj, size = xSize), axis.text.y = element_text(angle = yAngle, size = ySize), aspect.ratio = aspect.ratio, axis.title.x = element_blank(), axis.title.y = element_blank())+
    ggtitle(title)
  
  if(plot == TRUE){
    print(p)
  }
}


#Same as heatmap above, but now selecting only one of each .L and .S pair, whichever has the highest expression
heatmap_xenify <- function(countn, genes = rownames(countn), cells = 1:length(colnames(countn)), condition, conditionList = "none", norm = "max", aspect.ratio = 1, xAngle = 45, yAngle = 0, xSize = 10, ySize = 10, xAdj = 1, plot = TRUE){
  
  genes <- rev(genes)
  conditionList_internal <- unique(condition[cells])
  
  markers <- unique((genes))
  markers.L<-sapply(tolower((markers)),function(x) paste0(x,".L"))
  markers.S<-sapply(tolower(markers),function(x) paste0(x,".S"))
  
  
  M.L <- heatmap(countn,genes = markers.L,cells = cells,condition = condition, conditionList = conditionList_internal, norm = "none", plot = FALSE, return = T)
  M.S <- heatmap(countn,genes = markers.S,cells = cells,condition = condition, conditionList = conditionList_internal, norm = "none", plot = FALSE, return = T)
  
  
  max.L <- apply(M.L,1,max)
  max.S <- apply(M.S,1,max)
  
  newMarkers <- c()
  for (i in 1:length(markers)){
    
    mL <- match(markers.L[i],names(max.L))
    mS <- match(markers.S[i],names(max.S))
    
    if(!is.na(mL) & is.na(mS)){ 
      newMarkers <- c(newMarkers,markers.L[i])
    } else if(is.na(mL) & !is.na(mS)) {
      newMarkers <- c(newMarkers,markers.S[i])
    } else if(!is.na(mL) & !is.na(mS)){
      if(max.S[mS] > max.L[mL]){
        newMarkers <- c(newMarkers,markers.S[i])
      }
      else{
        newMarkers <- c(newMarkers,markers.L[i])
      }
    }
  }
  
  
  M <- heatmap(countn,genes = newMarkers,cells = cells,condition = condition, conditionList = conditionList, norm = "max", aspect.ratio = aspect.ratio, xSize = xSize, ySize = ySize, xAngle = xAngle, yAngle = yAngle)
}


#Same as heatmap above, but now selecting only one of each .L and .S pair, whichever has the highest expression
higherGenesXenify <- function(countn, genes = rownames(countn), cells = 1:length(colnames(countn)), condition, norm = "max", conditionList = "none", aspect.ratio = 1, xAngle = 45, yAngle = 0, xSize = 10, ySize = 10, xAdj = 1, plot = TRUE, 
                           title = "title", return = FALSE, normcells = 1:length(colnames(countn))){
  
  
  conditionList_internal <- unique(condition[cells])
  
  markers <- unique((genes))
  markers.L<-sapply(tolower((markers)),function(x) paste0(x,".L"))
  markers.S<-sapply(tolower(markers),function(x) paste0(x,".S"))
  
  
  M.L <- heatmap(countn,genes = markers.L,cells = cells,condition = condition, conditionList = conditionList_internal, norm = "none", plot = FALSE, return = T)
  M.S <- heatmap(countn,genes = markers.S,cells = cells,condition = condition, conditionList = conditionList_internal, norm = "none", plot = FALSE, return = T)
  
  
  max.L <- apply(M.L,1,max)
  max.S <- apply(M.S,1,max)
  
  newMarkers <- c()
  for (i in 1:length(markers)){
    
    mL <- match(markers.L[i],names(max.L))
    mS <- match(markers.S[i],names(max.S))
    
    if(!is.na(mL) & is.na(mS)){ 
      newMarkers <- c(newMarkers,markers.L[i])
    } else if(is.na(mL) & !is.na(mS)) {
      newMarkers <- c(newMarkers,markers.S[i])
    } else if(!is.na(mL) & !is.na(mS)){
      if(max.S[mS] > max.L[mL]){
        newMarkers <- c(newMarkers,markers.S[i])
      }
      else{
        newMarkers <- c(newMarkers,markers.L[i])
      }
    }
  }
  
  return(newMarkers)
}



plotGeneHex <- function(meta, countn, geneName, nbins, return = F){
  if(!(geneName %in% rownames(countn))){
    return()
  }
  drhex <- hexbin(meta$x,meta$y, xbins = nbins, shape = 1, IDs = TRUE)
  cID <- drhex@cID
  drhex <- data.frame(x = as.numeric(hcell2xy(drhex)$x), y= as.numeric(hcell2xy(drhex)$y))
  count_hex <- aggregate(log10(countn[geneName,]+1),by = list(cID),FUN = mean)
  drhex$plot <- count_hex$x
  p <- ggplot(drhex, mapping = aes(x, y, fill = plot, colour = plot)) +
    geom_hex(stat = "identity") +
    scale_fill_gradient2(low = "blue", mid = "gray", high = "red") +
    scale_color_gradient2(low = "blue", mid = "gray", high = "red") +
    
    theme_classic() + theme(legend.position = "bottom") +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", aspect.ratio = 1, axis.line = element_blank())  +
    ggtitle(geneName)
  print(p)

  if(return){return(drhex)}
}



findDEGeneSets <- function(meta,countn,subsetClusters,sample1,sample2,FDR, geneSets){
  subset <- intersect(which(meta$val %in% c(sample1,sample2)),which(meta$cluster %in% subsetClusters))
  out <- findMarkers(x = log2(countn[,subset]+1),clusters = meta$val[subset])
  
  dat <- (out[[1]])
  dat <- data.frame(name = rownames(dat),PValue =dat[,2],FDR =dat[,3],logFC =dat[,4])
  
   
  for(pathway in names(geneSets)) {
    
    dat <- dat[which(dat$name %in% unlist(geneSets)),]
    datS <- dat[which(dat$name %in% geneSets[[pathway]]),]
    datS2 <- datS[which(datS$FDR < FDR),]
    
    p <- ggplot(dat, aes(x=logFC, y=-log10(PValue))) + 
      theme_bw() +
      geom_point(size=1, col = "black") +
      geom_text_repel(data = datS2, mapping= aes(x = logFC, y = -log10(PValue), label=name), size=5, col = "black") +
      geom_point(data = datS, mapping= aes(x = logFC, y = -log10(PValue)), size=1, col = "red") +
      
      theme(  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1, legend.position = "none")  +
      theme(aspect.ratio = 1) +
      labs(x = "log fold change", y = "-log P-value") +
      ggtitle(paste0(pathway,":     ",sample2,"  vs.  ", sample1))
    print(p)
  }
}

violinPlot <- function(values, groups, title = "title", mode = "box"){

  df <- data.frame(score = as.vector(values),cluster = factor(groups))
  if(mode == "box"){
    p <- ggplot(df,aes(x = cluster, y=score)) +
      geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=0.4, notch=FALSE)+
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.title = element_blank()) + ggtitle(title)
    print(p)
  } else if(mode == "violin") {
    p <- ggplot(df,aes(x = cluster, y=score)) +
      geom_violin(outlier.colour="black")+
      theme_bw() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.title = element_blank()) + ggtitle(title) +
      # geom_jitter(shape=16, position=position_jitter(0.2)) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=.2, col = alpha("black",0.5))
    print(p)
    # p <- ggplot(df,aes(x = cluster, y=score)) +
    #   geom_violin(outlier.colour="black")+
    #   theme_bw() +
    #   theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.title = element_blank()) + ggtitle(title) +
    #   geom_jitter(shape=16, position=position_jitter(0.2), col = alpha("black",0.5)) 
    #   # geom_dotplot(binaxis='y', stackdir='center', dotsize=.1)
    # print(p)
  }
 

}

computeNNMF <- function(hvg,k,n.threads=0,max.iter=100, seed = 42){
  
  set.seed(seed)
  nsclc <- as.matrix(cosineNorm(log2(hvg+1)))
  init <- list(W = matrix(runif(nrow(nsclc)*k), ncol = k),
               H <-  matrix(runif(ncol(nsclc)*k), nrow = k));
  factors  <- nnmf(nsclc, k, init = init, max.iter = max.iter, n.threads = n.threads)
  rownames(factors$W) <- rownames(hvg)
  colnames(factors$H) <- colnames(hvg)
  return(factors)

}



#Plot points in 2D
plotPoints <- function(x, y, group, colours, size = 0.1, return = FALSE, title = "none", grid = FALSE ,lims = "none"){
  dat <- data.frame(x = x, y = y, group = group)
  p <- ggplot( data = dat, mapping= aes(x = x, y = y, col = group)) +
      geom_point(size=size) +
      theme_bw() +
      scale_colour_manual(values = sapply(colours$cols, as.character), labels= sapply(colours$names, as.character), name = "") +
    # scale_colour_manual(values=setNames(paste0(colours), group))+
    
      coord_equal() +
      theme(legend.position = "none")
  if(grid == FALSE)
  { p <- p + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())}
  if(title != "none") {p <- p + ggtitle(title)}
  if(lims != "none"){
    p <- p + xlim(lims[1:2]) + ylim(lims[3:4])
  }
  if(return){return(p)} else {print(p)}
}

#Plot points in 2D
plotCells <- function(x, y, group, colours, size = 0.1, lims = "none", grid = FALSE){
  
  dat <- data.frame(x = x, y = y, group = group)
  p <- ggplot( data = dat, mapping= aes(x = x, y = y, col = group)) +
    geom_point(size=size, stroke = 0) +
    theme_bw() +
    scale_colour_manual(values=setNames(paste0(colours), group)) +
    coord_equal() +  theme(legend.position = "none")
  if(grid == FALSE)
  { p <- p + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())}
  
  
  if(lims != "none"){
    p <- p + xlim(lims[1:2]) + ylim(lims[3:4])
  }
  
  print(p)
}

savePDF <- function(filename,scale = 1, width = 5, height = 5, dpi = 300){
  ggsave(paste0(filename,".pdf"), plot = last_plot(), device = "pdf", scale = scale, width = width, height = height, units = "in", dpi = dpi,  useDingbats = FALSE)
}

saveTIFF <- function(filename,scale = 1, width = 5, height = 5, dpi = 300){
  ggsave(paste0(filename,".tiff"), plot = last_plot(), device = "tiff", scale = scale, width = width, height = height, units = "in", dpi = dpi)
}

saveFig <- function(name = "name",mode = "pdf",scale = 1, width = 5, height = 5, dpi = 300){ 
  name = paste0("panels/", name)
  if(plot){
    if(mode == "pdf"){savePDF(name,scale = 1, width = 5, height = 5, dpi = 300)}
    else if (mode == "tiff"){saveTIFF(name,scale = 1, width = 5, height = 5, dpi = 300)}
  }
}
  
  
  increment <- function()
  {
    
  }
  
  




plotFactorGenes <- function(factors,fac,aspect.ratio = 1,nGenes = 20, ySize = 10, return = FALSE){
  
  
  factorGenes <- factors$W[,fac]
  
  factorGenes <- factorGenes[is.na(sapply(names(factorGenes),function(x) pmatch("Xelaev",x)))]
  factorGenes <- factorGenes[is.na(sapply(names(factorGenes),function(x) pmatch("unnamed",x)))]
  factorGenes <- factorGenes[is.na(sapply(names(factorGenes),function(x) pmatch("Xetrov",x)))]
  factorGenes <- factorGenes[is.na(sapply(names(factorGenes),function(x) pmatch("loc1",x)))]
  
  
  plotGenes <- sort(factorGenes, decreasing = T)[1:nGenes]
  
  dat <- data.frame(name = names(plotGenes), val = plotGenes)
  
  
  p <- ggplot(dat,aes(x = reorder(name, val),y = val)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_bw() +
    theme(aspect.ratio = aspect.ratio) +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.title = element_blank())  + theme(axis.text.y = element_text(size = ySize)) + ggtitle(fac)
  
  if(return){return(p)} else {print(p)}
  
}

plotStackedBar <- function(variable,groups,conditionList,colours, aspect.ratio = 1, title = "title", return = FALSE, width = 0.9){

  cellCycleProportions <- table(variable,groups)
  cellCycleProportions <- cellCycleProportions[,match(conditionList,colnames(cellCycleProportions))]
  
  dat <- melt(cellCycleProportions)
  dat$phase <- factor(dat$variable, sapply(colours$names, as.character))
  
  p <- ggplot(dat, aes(fill=phase, y=value, x=groups)) + 
    geom_bar(position="fill", stat="identity", width = width) +
    scale_fill_manual(values = sapply(colours$cols, as.character), labels= sapply(colours$names, as.character), name = "") +
    theme_bw() +
    theme(aspect.ratio = aspect.ratio) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.title = element_blank()) +
    ggtitle(title)
  if(return){return(p)} else{print(p)}
}

plotCellCycleHeatmap <- function(meta,conditionList, valueList, aspect.ratio = 1, xAngle = 45, yAngle = 0, xSize = 10, ySize = 10, xAdj = 1){
  g1 <- which(meta$CellCyclePhase == "G1")
  prolif <- which(meta$CellCyclePhase != "G1")
  
  G1 <- table(meta$cluster[g1],meta$val[g1])
  Prolif <- table(meta$cluster[prolif],meta$val[prolif])
  
  Frac <- Prolif/(G1 + Prolif)
  
  Frac <- Frac[match(conditionList,rownames(Frac)),match(valueList,colnames(Frac))]
  
  Frac[is.na(Frac)] <- 0
  
  
  
  dat <- melt(t(Frac))
  
  p <- ggplot(data = dat, aes(x=Var1, y=Var2, fill=value)) + 
    theme_bw() +
    geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")+
    theme(axis.text.x = element_text(angle = xAngle, hjust = xAdj, size = xSize), axis.text.y = element_text(angle = yAngle, size = ySize), aspect.ratio = aspect.ratio, axis.title.x = element_blank(), axis.title.y = element_blank()) +
    ggtitle("Fraction of proliferating cells")
  print(p)
}

heatmapTwoGroups <- function(values, var1, var2,conditionList, valueList,aspect.ratio = 1, xAngle = 45, yAngle = 0, xSize = 10, ySize = 10, xAdj = 1, title = "title"){
  Agg <- aggregate(values, by = list(var1,var2), FUN = mean)
  Mean <- dcast(Agg, Group.1 ~ Group.2)
  Mean <- Mean[match(conditionList,Mean$Group.1),match(valueList,colnames(Mean))]
  rownames(Mean) <- conditionList
  Mean[is.na(Mean)] <- 0
  dat <- melt(t(Mean))
  p <- ggplot(data = dat, aes(x=Var1, y=Var2, fill=value)) +
    theme_bw() +
    geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")+
    theme(axis.text.x = element_text(angle = xAngle, hjust = xAdj, size = xSize), axis.text.y = element_text(angle = yAngle, size = ySize), aspect.ratio = aspect.ratio, axis.title.x = element_blank(), axis.title.y = element_blank()) +
    ggtitle(title)
  print(p)
}

dotplotTwoGroups <- function(values, var1, var2,conditionList, valueList,aspect.ratio = 1, xAngle = 45, yAngle = 0, xSize = 10, ySize = 10, xAdj = 1, title = "title", plot = T, return = F, collow = "lightgrey", colhigh = "blue", dot.scale = 2){
  Agg <- aggregate(values, by = list(var1,var2), FUN = mean)
  Mean <- dcast(Agg, Group.1 ~ Group.2)
  Mean <- Mean[match(conditionList,Mean$Group.1),match(valueList,colnames(Mean))]
  rownames(Mean) <- conditionList
  Mean[is.na(Mean)] <- 0
  
  Agg2 <- aggregate(values, by = list(var1,var2), FUN = function(x) {return(mean(x>0))})
  Pct <- dcast(Agg2, Group.1 ~ Group.2)
  Pct <- Pct[match(conditionList,Pct$Group.1),match(valueList,colnames(Pct))]
  rownames(Pct) <- conditionList
  Pct[is.na(Mean)] <- 0
  
  
  dat <- melt(t(Mean))
  
  pct <- melt(t(Pct))
  dat$pct <- 100*pct$value
  
  
  p <- ggplot(data = dat) + 
    theme_bw() +
    geom_point(mapping = aes(x = Var1, y = Var2,size = pct,col = value))+
    scale_color_gradient(low = collow,high = colhigh)+
    scale_size(range = c(0,dot.scale), limits = c(0, 100)) +
    theme(axis.text.x = element_text(angle = xAngle, hjust = xAdj, size = xSize), axis.text.y = element_text(angle = yAngle, size = ySize), aspect.ratio = aspect.ratio, axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(panel.background = element_blank()) + theme(panel.grid = element_blank())
  
  p <- p + ggtitle(title)
  if(plot == TRUE){
    print(p)
  }
  if (return) {return(avgd_n)}
  

}

plotCoexpressedFactors <- function(factors,meta,cells,conditionList,fac1,fac2,plot)
  
  for(val in amputation){
    gene1 <- paste0(factorNames$gene[factorNames$name == fac1])
    gene2 <- paste0(factorNames$gene[factorNames$name == fac2])
    maxDat <- data.frame(fac1 = max(factors$H[fac1,cells]),fac2 = max(factors$H[fac2,cells]),gene1 = max(log10(countn[gene1,cells]+1)),gene2 = max(log10(countn[gene1,cells]+1)))
    
    cells_val <- intersect(cells, which(meta$val == val))
    
    
    
    dat <- data.frame(fac1 = factors$H[fac1,cells_val], fac2 = factors$H[fac2,cells_val], gene1 = log10(countn[gene1,cells_val]+1), gene2 = log10(countn[gene2,cells_val]+1))  
    
    
    
    p <- ggplot(dat,mapping = aes(x = fac1,y = fac2)) +
      geom_point(size = .2, col = alpha("black",0.5)) +
      theme_bw() +
      theme(  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1, legend.position = "none")  + 
      theme(aspect.ratio = 1) + ggtitle(val) + xlab(fac1) + ylab(fac2) + xlim(c(0,maxDat$fac1)) + ylim(c(0,maxDat$fac2))
    print(p)
    saveFig(mode = "tiff")
    
    p <- ggplot(dat,mapping = aes(x = gene1,y = gene2)) +
      geom_point(size = .2, col = alpha("black",0.5)) +
      theme_bw() +
      theme(  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1, legend.position = "none")  + 
      theme(aspect.ratio = 1) + ggtitle(val) + xlab(gene1) + ylab(gene2) + xlim(c(0,maxDat$gene1)) + ylim(c(0,maxDat$gene2))
    print(p)
    saveFig(mode = "tiff")
  }


plotCoexpressedFactorsGrid <- function(factors,meta,cells,conditionList,fac1,fac2,plot){
  genePlots <- list()
  factorPlots <- list()
  for(val in amputation){
    gene1 <- paste0(factorNames$gene[factorNames$name == fac1])
    gene2 <- paste0(factorNames$gene[factorNames$name == fac2])
    maxDat <- data.frame(fac1 = max(factors$H[fac1,cells]),fac2 = max(factors$H[fac2,cells]),gene1 = max(log10(countn[gene1,cells]+1)),gene2 = max(log10(countn[gene1,cells]+1)))
    cells_val <- intersect(cells, which(meta$val == val))
    dat <- data.frame(fac1 = factors$H[fac1,cells_val], fac2 = factors$H[fac2,cells_val], gene1 = log10(countn[gene1,cells_val]+1), gene2 = log10(countn[gene2,cells_val]+1))  
    factorPlots[[val]] <- ggplot(dat,mapping = aes(x = fac1,y = fac2)) +
      geom_point(size = .2, col = alpha("black",0.5)) +
      theme_bw() +
      theme(  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1, legend.position = "none")  +
      theme(axis.title = element_blank(), axis.text = element_blank()) +
      theme(aspect.ratio = 1) + xlim(c(0,maxDat$fac1)) + ylim(c(0,maxDat$fac2)) +
      theme(plot.margin=unit(rep(0.2,4), "cm"))
    genePlots[[val]] <- ggplot(dat,mapping = aes(x = gene1,y = gene2)) +
      geom_point(size = .2, col = alpha("black",0.5)) +
      theme_bw() +
      theme(  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio=1, legend.position = "none")  + 
      theme(aspect.ratio = 1) + xlim(c(0,maxDat$gene1)) + ylim(c(0,maxDat$gene2)) +
      theme(axis.title = element_blank(), axis.text = element_blank()) +
      theme(plot.margin=unit(rep(0.2,4), "cm"))

  }

grid.arrange(grobs = factorPlots,widths = c(1,1,1), layout_matrix = rbind(c(1,3,NA),c(2,4,6)), top = paste0(fac1, " (x-axis) and ", fac2, " (y-axis)"))
g <- arrangeGrob(grobs = factorPlots,widths = c(1,1,1), layout_matrix = rbind(c(1,3,NA),c(2,4,6)), top = paste0(fac1, " (x-axis) and ", fac2, " (y-axis) "))

ggsave(paste0("panels/S2C_C",fac2,".tiff"), plot = g, device = "tiff", width = 5, height = 5, dpi = 300)




grid.arrange(grobs = genePlots,widths = c(1,1,1), layout_matrix = rbind(c(1,3,NA),c(2,4,6)), top = paste0("Coexpression of ", gene1, " (x-axis) and ", gene2, " (y-axis) across: ", unique(meta$group[cells])))
g <- arrangeGrob(grobs = genePlots,widths = c(1,1,1), layout_matrix = rbind(c(1,3,NA),c(2,4,6)), top = paste0("Coexpression of ", gene1, " (x-axis) and ", gene2, " (y-axis) across: ", unique(meta$group[cells])))

ggsave(paste0("panels/S2C_C",fac2,"_genes.tiff"), plot = g, device = "tiff", width = 5, height = 5, dpi = 300)



}



differentialExpression <- function(countn,condition, sample1, sample2, subset, genes = rownames(countn)){
  
  
  cells <- intersect(c(sample1,sample2),subset)
  
  
  out <- findMarkers(x = log2(countn[,cells]+1),clusters = condition[cells])
  dat <- (out[[1]])
  dat <- data.frame(name = rownames(dat),PValue =dat[,2],FDR =dat[,3],logFC =dat[,4])
  
  return(dat[which(dat$name %in% genes),])
}