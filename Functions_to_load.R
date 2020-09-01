###Packages and Functions to laod for figures

###packages

require(RaceID)
require(FateID)
require(RColorBrewer)
require(stringr)
library(clusterProfiler)
library("org.Mm.eg.db")
library(ReactomePA)
library(scran)
library(scater)
library(Matrix)
library(Seurat)
library(cluster)

##library size normalization seurat count object
getRC_object <- function(object) {
  colsum <- Matrix::colSums(object@assays$RNA@data)
  object@assays[["RC"]] <- object@assays$RNA
  object@assays$RC@data <- t(t(as.matrix(object@assays$RNA@data))/colsum)
  return(object)
}  
getRC_object <- function(object) {
  colsum <- Matrix::colSums(object@assays$RNA@data)
  object@assays[["RC"]] <- object@assays$RNA
  object@assays$RC@data <- t(t(as.matrix(object@assays$RNA@data))/colsum)
  return(object)
}  
##fraction dotplot
fracdotplot <- function( object, genes, cluster=NULL, sampclus = NULL, samples=NULL, limup, limdown, zsc=T, map=T, leg=T, mintotal=NULL, seurat=F) {
  library(ggplot2)
  library(RColorBrewer)
  if (seurat) {
    if ( is.null(object@assays$RC@data)) { 
      stop("compute relative counts matrix for seurat object")}
    ndata <- as.matrix(object@assays$RC@data * mintotal) + 0.1
    ndata <- data.frame(ndata)
  }
  else {
    if ( length(object@cpart) == 0) {
      object@cpart <- object@cluster$kpart
    }
    if ( ! identical(names(object@cpart), colnames(object@ndata)) ) { names(object@cpart) <- colnames(object@ndata)}
    ndata <- as.matrix(object@ndata * min(object@counts)) + 0.1
    ndata <- data.frame(ndata)}
  
  if ( !is.null(cluster) & !is.null(samples)) {
    stop("define either clusters OR samples")
  }
  if (is.null(cluster) & is.null(samples)) {
    stop("define either clusters OR samples")
  }
  if (!is.null(cluster)) {
    genevec <- c()
    clustervec <- c()
    fraction <- c()
    scaled_mean <- c()
    log2mean <- c()
    if (seurat) {
      cpart <- object$seurat_clusters
    }
    else {
      cpart <- object@cpart
    }
    for ( i in 1:length(genes)) {
      repgene <- rep(genes[i], length(cluster))
      
      if (!is.null(sampclus)) {
        meang <- mean(as.numeric(ndata[genes[i], grep(sampclus, colnames(ndata))]))
        sdg <- sd(ndata[genes[i],grep(sampclus, colnames(ndata))])
      }
      else {
        meang <- mean(as.numeric(ndata[genes[i],]))
        sdg <- sd(ndata[genes[i],]) }
      repclus <- c()
      frac <- c()
      cent_mean <- c()
      log2_mean <- c()
      for ( n in 1:length(cluster)) {
        if (!is.null(sampclus)) {
          clus <- names(cpart[as.numeric(cpart) == cluster[n]])
          clus <- clus[grep(sampclus, clus)]
        }
        else {
          clus <- names(cpart[as.numeric(cpart) == cluster[n]])}
        leng_clus <- length(clus)
        leng_gene_in_clus <- length(which(ndata[genes[i], clus] > 0.1))
        frac <- c(frac, leng_gene_in_clus/leng_clus)
        #repclus <- c(repclus, paste("cl",cluster[n], sep="_"))
        repclus <- c(repclus, cluster[n])
        if (zsc) { cent_mean <- c(cent_mean, (mean(as.numeric(ndata[genes[i], clus])) - meang)/sdg)}
        else { log2_mean <- c(log2_mean, log2(mean(as.numeric(ndata[genes[i], clus]))))}
      }
      genevec <- c(genevec, repgene) 
      clustervec <- c(clustervec, repclus)
      fraction <- c(fraction, frac)
      if (zsc) { scaled_mean <- c(scaled_mean, cent_mean)}
      else { log2mean <- c(log2mean, log2_mean) }
    }
    if ( zsc==T ) {
      data <- data.frame(Gene = factor(genevec, levels = genes) , Cluster = factor(clustervec, levels = cluster), Fraction = fraction, Expression = scaled_mean )
    }
    else {
      data <- data.frame(Gene = factor(genevec, levels = genes) , Cluster = factor(clustervec, levels = cluster), Fraction = fraction, Expression = log2mean )  
    }
    data[which(data$Expression > limup), "Expression"] <- limup
    data[which(data$Expression < limdown), "Expression"] <- limdown
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
    
    frac <- ggplot(data, aes(x = Gene, y = Cluster))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank()) + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())                                                                    
    
    if ( map==T && leg==F) {
      print(frac + geom_point(aes(size = Fraction, color = Expression))  + scale_colour_gradientn(colours = ColorRamp) + theme(legend.position="none"))
    }
    if (leg==T && map==F){
      print(frac + theme(axis.line = element_line(colour = "black")) + theme(axis.title = element_text(color = "black"), axis.text = element_text(color = "black"),axis.ticks = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1)))
    }
    if (map==T && leg==T) {
      print(frac + geom_point(aes(size = Fraction, color = Expression))  + scale_colour_gradientn(colours = ColorRamp) + theme(axis.line = element_line(colour = "black")) + theme(axis.title = element_text(color = "black"), axis.text = element_text(color = "black"),axis.ticks = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1)) )
    }
  }
  
  if (!is.null(samples)) {
    genevec <- c()
    samplevec <- c()
    fraction <- c()
    fraction_log <- c()
    scaled_mean <- c()
    log2mean <- c()
    for ( i in 1:length(genes)) {
      repgene <- rep(genes[i], length(samples))
      meang <- mean(as.numeric(ndata[genes[i],]))
      sdg <- sd(ndata[genes[i],])
      repsamp <- c()
      frac <- c()
      cent_mean <- c()
      log2_mean <- c()
      for ( n in 1:length(samples)) {
        samp <- colnames(ndata)[grep(samples[n], colnames(ndata))]
        leng_samp <- length(samp)
        leng_gene_in_samp <- length(which(ndata[genes[i], samp ]> 0.1))
        frac <- c(frac, leng_gene_in_samp/leng_samp)
        repsamp <- c(repsamp, samples[n])
        if ( zsc ) { cent_mean <- c(cent_mean, (mean(as.numeric(ndata[genes[i], samp])) - meang)/sdg) } 
        else { log2_mean <- c(log2_mean, log2(mean(as.numeric(ndata[genes[i], samp])))) } 
      }
      genevec <- c(genevec, repgene) 
      samplevec <- c(samplevec, repsamp)
      fraction <- c(fraction, frac)
      if ( zsc ) { scaled_mean <- c(scaled_mean, cent_mean) }
      else { log2mean <- c(log2mean, log2_mean) } 
    }
    if (zsc==T) {
      data <- data.frame(Gene = factor(genevec, levels = genes) , Sample = factor(samplevec, levels = samples), Fraction = fraction, Expression = scaled_mean )
    }
    else {
      data <- data.frame(Gene = factor(genevec, levels = genes) , Sample = factor(samplevec, levels = samples), Fraction = fraction, Expression = log2mean )  
    }
    data[which(data$Expression > limup), "Expression"] <- limup
    data[which(data$Expression < limdown), "Expression"] <- limdown
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
    
    frac <- ggplot(data, aes(x = Gene, y = Sample))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank()) + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())                                                                    
    
    if ( map==T && leg==F) {
      print(frac + geom_point(aes(size = Fraction, color = Expression))  + scale_colour_gradientn(colours = ColorRamp) + theme(legend.position="none"))
    }
    if (leg==T && map==F){
      print(frac + theme(axis.line = element_line(colour = "black")) + theme(axis.title = element_text(color = "black"), axis.text = element_text(color = "black"),axis.ticks = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1)))
    }
    if (map==T && leg==T) {
      print(frac + geom_point(aes(size = Fraction, color = Expression))  + scale_colour_gradientn(colours = ColorRamp) + theme(axis.line = element_line(colour = "black")) + theme(axis.title = element_text(color = "black"), axis.text = element_text(color = "black"),axis.ticks = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1)) )
    }
  }
}


##enrichment and piecharts 
get_enrichment <- function(object, sample, clustsize = 20, alternat = "greater") { 
  sample <- sample[order(sample)]
  if (is.null(object@cpart)) {
    object@cpart <- object@cluster$kpart
  }
  if ( length(sample) > 9)
    stop( "too many samples, use x <= 9")
  types <- unique(sub("\\_.+", "", names(object@cpart)))
  if ( sum(sample %in% types) != length(sample))
    stop("samples must be element of colnames")
  pop <- length(object@cpart)
  cluster <- as.numeric(names(table(object@cpart)[table(object@cpart) > clustsize]))
  pvals <- list()
  for ( k in 1:length(sample)) {
    size_samp <- sum(grepl(sample[k], names(object@cpart)))
    cat(paste(size_samp, "\n",sep=""))
    pvalues <- c()
    for ( i in 1:length(cluster)) {
      clusteri <- object@cpart[object@cpart == cluster[i]]
      size_cluster <- length(clusteri)
      size_samp_clust <- sum(grepl(sample[k], names(clusteri)))
      fisher <- fisher.test(matrix(c(pop - size_samp, size_samp,size_cluster - size_samp_clust, size_samp_clust),ncol=2), alternative = alternat)
      pv <- fisher$p.value
      pvalues <- append(pvalues, pv)
      names(pvalues)[i] <- paste(k, "_cl", cluster[i], sep="")
      
    }
    pvals[[k]] <- pvalues
    
  }
  pvals2 <- pvals
  names(pvals2) <- sample
  
  if ( length(sample) >= 3 ){
    marker_col <- brewer.pal(length(sample), "Set1")}
  
  else {
    marker_col <- brewer.pal(3, "Set1")
    marker_col <- marker_col[1:length(sample)]
  }
  
  types <- sub("(\\_|\\.).+","", colnames(object@ndata))
  types_num <- as.numeric(table(types))
  types <- names(table(types))
  
  rel_size <- types_num/pop
  rel_size_abs <- rel_size *100
  
  values <- list()
  
  
  a <- cluster
  
  values <- list()
  for ( i in 1:length(a)){
    enrichment <- c()
    x <- table(as.character(sapply(names(object@cpart[object@cpart==a[i]]), function(x) str_split(x, pattern="[_]")[[1]][1])))
    
    lbls <- names(x)
    cat(paste(lbls, "\n", sep=""))
    y <- which(types %in% lbls)  ## index for color
    slices <- as.numeric(x)/types_num[y] ### normalize with number of total cells of sample
    pct <- round(slices/sum(slices)*100)
    lbls <- paste(pct, " % ", lbls,sep=" ")
    enrichment <- pct/rel_size_abs[y]
    names(enrichment) <- lbls
    values[[i]] <- enrichment
    names(values)[i] <- paste("cluster", a[i], sep="")
    #pdf(file.path(".", paste("cluster_norm_", alternat, "_", a[i],"distribution",".pdf", sep="")))
    print(pie(slices, labels=lbls,cex.main=1.5, cex.lab=0.75,col=marker_col[y], main=paste("Cluster ", a[i], " distribution ", " N=", length(object@cpart[object@cpart == a[i]]), sep=" ")))
    print(barplot(enrichment,  names.arg=c(round(-log10(pvals2[[1]][i]), 2), round(-log10(pvals2[[2]][i]),2)), main=paste("Cluster ", a[i], " enrichment ", sep=""), col = marker_col, ylim = range(0:2)))
    legend("topright", names(x), fill=marker_col[y], cex=1, bty="n")
    #dev.off()
  }
}

##tSNE map for weights
setGeneric("plotexptsne", function(object,g,n="",logsc=FALSE, map=T, leg=T, cex=1) standardGeneric("plotexptsne"))

setMethod("plotexptsne",
          signature = "SCseq",
          definition = function(object,g,n="",logsc=FALSE, map=T, leg=T, cex=1){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
            if ( length(intersect(g,rownames(object@ndata))) < length(unique(g)) ) stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
            if ( !is.numeric(logsc) & !is.logical(logsc) ) stop("argument logsc has to be logical (TRUE/FALSE)")
            if ( n == "" ) n <- g[1]
            l <- apply(object@ndata[g,] - .1,2,sum) + .1
            if (logsc) {
              f <- l == 0
              l <- log2(l)
              l[f] <- NA
            }
            cells <- colnames(object@ndata)
            h <- colnames(object@ndata) %in% cells
            mi <- min(l,na.rm=TRUE)
            ma <- max(l,na.rm=TRUE)
            ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
            ColorLevels <- seq(mi, ma, length=length(ColorRamp))
            v <- round((l - mi)/(ma - mi)*99 + 1,0)
            d <- object@tsne
            layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
            
            
            
            pardefault <- par()
            if (!leg) 
              n <- NA
            plot(c(min(d[, 1]), max(d[, 1])), c(min(d[, 2]), max(d[, 
                                                                   2])), xlab = NA, ylab = NA, main = n, pch = 20, cex = 0, 
                 col = "lightgrey", axes = FALSE)
            if (map) {
              v <- v[h]
              d <- d[h, ]
              kk <- order(v, decreasing = F)
              points(d[kk, 1], d[kk, 2], col = ColorRamp[v[kk]], pch = 20, 
                     cex = cex)
            }
            if (leg) {
              par(mar = c(10, 2.5, 2.5, 4))
              image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), 
                                           nrow = 1), col = ColorRamp, xlab = "", ylab = "", 
                    xaxt = "n")
              layout(1)
              par(mar = pardefault$mar)
            }
          }
)

##plotexpression mutiple genes for pseudo-temoral ordered cells
plotexpression <- function (x, y, g, n, col = NULL, name = NULL, cluster = FALSE, 
                            alpha = 0.5, types = NULL, cex = 3, ylim = NULL, leg = T,
                            lwd = 2.5, legendpos="topleft", samp=F, samppart=NULL, sampcol=NULL) {
  library(RColorBrewer)
  
  if (length(g) <= 9) {
    set1 <- brewer.pal(length(g), "Set1")
    if (length(g) >= 6) {
      set2 <- brewer.pal(length(g), "Set2") 
      set1[6] <- set2[6]} 
  }
  else { 
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    set1 <- sample(col_vector, length(g))
  }
  
  for ( i in 1:length(g)){
    cl <- unique(y[n])
    set.seed(111111)
    if (is.null(col)) 
      col <- sample(rainbow(max(y)))
    xlim <- c(1, length(n))
    
    z <-  x[g[i], n]
    if (leg) {
      ylab = "Expression"
    }
    else {
      ylab = NA
    }
    if (!is.null(name)) {
      #main = name
      main = name
    }
    else {
      main = NA
    }
    if (i == 1){
      if (is.null(ylim)) {
        plot(c(1, length(n)), c(min(z), max(z)), cex = 0, axes = FALSE, 
             xlab = "", ylab = ylab, main = main, xlim = xlim)
      }
      else {
        plot(c(1, length(n)), c(min(z), max(z)), cex = 0, axes = T, 
             xlab = "", ylab = ylab, main = main, xlim = xlim, 
             ylim = ylim, xaxt='n')
      }}
    
    
    u <- 1:length(n)
    v <- as.vector(t(z))
    zc <- predict(loess(v ~ u, span = alpha))
    zc[zc < 0] <- 0.1
    lines(u, zc, lwd = lwd, col=set1[i])
    
  }
  width <- .45
  
  k <- 1:length(y)
  
  rect(xleft = k-width, ybottom = rep(-0.2, length(y)), xright =  k + width, ytop = rep(0, length(y)), col = col[y], border = NA)
  
  if ( samp == T) {
    if ( is.null(samppart)){
      stop("set samp part ")
    }
    rect(xleft = k-width, ybottom = rep(-0.5, length(y)), xright =  k + width, ytop = rep(-0.3, length(y)), col = sampcol[samppart], border = NA)
  }
  
  if (!leg) 
    box(col = "white")
  else legend(legendpos, g, col=set1, pch=20)
}

## plot3marker 
plot3marker <- function(object, logsc=T, cluster=cluster, gene1, gene2, gene3, leg=T, map=T, xlim=NULL, ylim=NULL){
  if (!(length(cluster) == 1 || length(cluster) == 2)){
    stop("set either one or two clusters")
  }
  mat <- as.matrix(object@ndata* min(object@counts)) + 0.1
  if(length(cluster) == 1) {
    cl <- names(object@cpart)[object@cpart == cluster]
    l <- as.numeric(mat[gene3, cl])
    m1 <- data.frame(cbind(as.numeric(mat[gene1,cl]), as.numeric(mat[gene2,cl])))
    if (logsc==T) {
      f <- l == 0
      l <- log2(l)
      l[f] <- NA
    }
    
    mi <- min(l, na.rm = TRUE)
    ma <- max(l, na.rm = TRUE)
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    ColorLevels <- seq(mi, ma, length = length(ColorRamp))
    v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
    kk <- order(v, decreasing = F)
  }
  else {
    plot3_gene3 <- list()
    plot3_mat <- list()
    mi <- c()
    ma <- c()
    xlim <- c()
    ylim <- c()
    for (i in 1:length(cluster)) {
      cl <- names(object@cpart)[object@cpart == cluster[i]]
      xlim <- c(xlim, max(as.numeric(mat[gene1,cl])))
      ylim <- c(ylim, max(as.numeric(mat[gene2,cl])))
      plot3_gene3[[i]] <- as.numeric(mat[gene3,cl])
      if ( logsc==T) {
        f <- plot3_gene3[[i]] == 0
        plot3_gene3[[i]] <- log2(plot3_gene3[[i]])
        plot3_gene3[[i]][f] <- NA}
      plot3_mat[[i]] <- data.frame(cbind(as.numeric(mat[gene1,cl]), as.numeric(mat[gene2,cl])))
      mi <- c(mi, min(plot3_gene3[[i]],  na.rm=T))
      ma <- c(ma, max(plot3_gene3[[i]],  na.rm=T))
    }
    mi <- min(mi)
    ma <- max(ma)
    xlim <- max(xlim) + 1
    ylim <- max(ylim) + 1
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    ColorLevels <- seq(mi, ma, length = length(ColorRamp))
    list_v <- list()
    list_kk <- list()
    for ( i in 1:length(cluster)) {
      list_v[[i]] <- round((plot3_gene3[[i]] - mi)/(ma - mi) * 99 + 1, 0)
      list_kk[[i]] <- order(list_v[[i]], decreasing = F)
    }
  }
  if ( length(cluster) == 1) {
    pardefault <- par()
    layout(matrix(data = c(1, 3, 2, 4), nrow = 2, ncol = 2), 
           widths = c(5, 1, 5, 1), heights = c(5, 1, 1, 1))
    par(mar = c(4, 5, 2.5, 2))
    plot(m1[kk,1], m1[kk,2], xlab = gene1, ylab = gene2, main = paste("expression in cluster", cluster, sep="") ,col = ColorRamp[v[kk]],pch =20, cex=6) #}
  }
  else {
    pardefault <- par()
    layout(matrix(c(0, 0, 0,
                    0, 0, 0,
                    1, 2, 3), nrow=3, byrow=TRUE), heights = c(1, 1, 1.5), widths = c(3,3,1))
    par(mar = c(4, 5, 2.5, 2))
    for ( i in 1:length(cluster)) { 
      kk <- list_kk[[i]]
      v <- list_v[[i]]
      plot(plot3_mat[[i]][kk,1], plot3_mat[[i]][kk,2], xlab = gene1, ylab = gene2, ylim = range(0:ylim), xlim = range(0:xlim), main = paste("expression in cluster", cluster[i], sep="") ,col = ColorRamp[v[kk]],pch =20, cex=3)
      
    }
    
  }
  par(mar = c(10, 2.5, 2.5, 4))
  image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), 
                               nrow = 1), col = ColorRamp, xlab = gene3, ylab = "", 
        xaxt = "n")
  layout(1)
  par(mar = pardefault$mar)
}

## cluster and sample composition 

get_cluster_composition <- function(object, cluster_size, norm=T, color = NULL, map = T, leg=T, final=F, order=NULL, seurat=F, symbols=NULL, samples=NULL, fcol=NULL){
  if ( seurat ) {
    types <- unique(sub("\\_.+", "", names(object$seurat_clusters)))
    types <- types[order(types)]
    cluster <- as.numeric(names(table(as.numeric(object$seurat_clusters))[table(as.numeric(object$seurat_clusters)) >= cluster_size]))
    f <- as.numeric(object$seurat_clusters) %in% cluster
    cpart <- object$seurat_clusters[f]
  }
  else {
    if ( final==F) {
      object@cpart <- object@cluster$kpart
    }
    if (!is.null(symbols)) {
      colnames(object@ndata) <- symbols
    }
    if ( !identical(colnames(object@ndata), names(object@cpart))) {
      names(object@cpart) <- colnames(object@ndata)
    }  
    
    cluster <- as.numeric(names(table(object@cpart)[table(object@cpart) >= cluster_size]))
    types <- unique(sub("\\_.+", "", colnames(object@ndata)))
    types <- types[order(types)]
    
    cpart <- object@cpart[as.numeric(object@cpart) %in% cluster ]
  }
  if (is.null(color)) {
    color <- rainbow(length(types))
  }
  iels <- cbind(cluster=as.numeric(cpart),sample=sub("\\_.+", "", names(cpart)))
  rownames(iels) <- names(cpart)
  iels <- data.frame(iels)
  if (!is.null(samples)) {
    counts <- as.matrix(table( iels$cluster,iels$sample))
  }
  else { counts <- as.matrix(table( iels$sample,iels$cluster))}
  #
  if ( norm == T){
    if (seurat) { counts <- counts/as.numeric(table(sub("\\_.+", "", colnames(object@assays$RNA))))}
    else { counts <- counts/as.numeric(table(sub("\\_.+", "", colnames(object@ndata))))}  # normalize
  }
  
  rel_counts <- t(t(counts)/apply(counts,2,sum))
  if (is.null(samples)) {
    rel_counts <- rel_counts[,order(as.numeric(colnames(rel_counts)))] }
  if ( !is.null(order)) {
    rel_counts <- rel_counts[,as.character(order)]
  }
  #rel_counts <- rel_counts[,cluster]
  
  if (map == F && leg==T) {
    #pdf(file.path(".", paste(name, ".pdf")))
    barplot(rel_counts, col = NA, border = NA, axes = FALSE, axisnames = F)
    if (is.null(samples)) {
      legend( "topright", pch=20, bty="n", cex=1, legend=types, col=color) }
    else{legend( "topright", pch=20, bty="n", cex=1, legend=paste("cluster", rownames(rel_counts), sep="."), col=color) }
    #dev.off()
  }
  
  else {
    #pdf(file.path(".", paste(name, ".pdf")))
    barplot(rel_counts, main="Sample Contribution to RaceID3 cluster", ylab="% of cluster",
            xlab="", col=color, names.arg=as.character(colnames(rel_counts)), cex.names=1, las=2)
    if ( leg ==T) {
      
      legend( "topright", pch=20, bty="n", cex=1, legend=types, col=color) 
      
    }
    #dev.off()
  }
}

## GO term genes of trajectory modules
go_analysis <- function (nodes, uni, node=NULL) {
  y <- nodes
  if (!is.null(node)){
    if (!length(node) == 1 ) stop("please set 1 node")
    uni <- bitr(uni, fromType = "SYMBOL", toType = c("ENTREZID", "GO"),OrgDb = org.Mm.eg.db, drop=F)
    uni <- uni$ENTREZID
    uni <- uni[!is.na(uni)]
    uni <- unique(uni)
    y2 <- y[y==node]
    z <- names(y2)
    g <- z
    g <- bitr(g, fromType = "SYMBOL", toType = c("ENTREZID", "GO"),OrgDb = org.Mm.eg.db, drop=F)
    g <- g$ENTREZID
    g <- g[!is.na(g)]
    g <- unique(g)
    for (l in c("BP", "CC", "MF")) {
      ego <- enrichGO(gene          = g,
                      universe      = uni,
                      OrgDb         = org.Mm.eg.db,
                      ont           = l,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
      if ( ! nrow(ego) == 0) {
        print(dotplot(ego))
      }
      else{ cat(paste("no go enrichment", l, "\n"))}
    }
    ereactome <- enrichPathway(gene = g,
                               universe = uni,
                               organism = "mouse",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05,
                               pAdjustMethod = "BH")
    if (!nrow(ereactome) == 0) {
      print(dotplot(ereactome))
    }
    else{ cat(paste("no reactome enrichment", "\n"))}
    ekeg <- enrichKEGG(gene  = g, 
                       universe = uni,
                       organism = "mmu",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       pAdjustMethod = "BH"
    )
    if (! nrow(ekeg) == 0) {
      print(dotplot(ekeg))
    }
    else{cat(paste("no Kegg enrichment", "\n"))}
  }
  if ( is.null(node)) {
    for ( i in 1:length(table(nodes))) {  
      
      
      ##for GO apply p-value cutoff
      
      uni <- bitr(uni, fromType = "SYMBOL", toType = c("ENTREZID", "GO"),OrgDb = org.Mm.eg.db, drop=F)
      uni <- uni$ENTREZID
      uni <- uni[!is.na(uni)]
      uni <- unique(uni)
      y2 <- y[y==i]
      cat(paste("nodes", i, sep="_"))
      z <- names(y2)
      g <- z
      g <- bitr(g, fromType = "SYMBOL", toType = c("ENTREZID", "GO"),OrgDb = org.Mm.eg.db, drop=F)
      g <- g$ENTREZID
      g <- g[!is.na(g)]
      g <- unique(g)
      for (l in c("BP", "CC", "MF")) {
        ego <- enrichGO(gene          = g,
                        universe      = uni,
                        OrgDb         = org.Mm.eg.db,
                        ont           = l,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)
        if (! nrow(ego) == 0) {
          pdf(file.path(".", paste("GOenrich_nodes", i, "_", l, ".pdf")), 12, 12)
          print(dotplot(ego))
          dev.off() }
        else{ cat(paste("no go enrichment", l, "\n"))}
        
      }
      ereactome <- enrichPathway(gene = g,
                                 universe = uni,
                                 organism = "mouse",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05,
                                 pAdjustMethod = "BH")
      if (!nrow (ereactome) == 0) {
        pdf(file.path(".", paste("Reactomeenrich_", i, ".pdf")), 12, 12)
        print(dotplot(ereactome))
        dev.off()}
      else{ paste(cat("no reactome enrichment", "\n"))}
      
      ekeg <- enrichKEGG(gene  = g, 
                         universe = uni,
                         organism = "mmu",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         pAdjustMethod = "BH"
      )
      if ( ! nrow(ekeg) == 0) {
        pdf(file.path(".", paste("Keggenrich_", i, ".pdf")), 12, 12)
        print(dotplot(ekeg))
        dev.off()}
      else{ paste(cat("no Kegg enrichment", "\n"))}
      
      
    }
  }
}

## plotexpmap for Seurat

plotexpmap2 <- function (object, g, n = NULL, logsc = FALSE, imputed = FALSE, 
                         fr = FALSE, um = FALSE, cells = NULL, cex = 1, map = TRUE, 
                         leg = TRUE, noise = FALSE, seurat=F, mintotal=NULL, reduction="umap" ) 
{
  if (seurat) {
    if ( is.null(mintotal)){
      stop("if seurat = T, mintotal as scaling factor has to be set")
    }
    DefaultAssay(object) <- "RC"
    
    
    l <- as.vector(object@assays$RC@data[g, ] * mintotal + 
                     0.1)
    if (is.null(cells)) {
      cells <- colnames(object@assays$RC@data)}
    if (is.null(n)) 
      n <- g[1]
    if (logsc) {
      f <- l == 0
      l <- log2(l)
      l[f] <- NA
    }
    h <- colnames(object@assays$RC@data) %in% cells
  }
  else {
    if (length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 
        0) 
      stop("run comptsne/compfr/compumap before plotlabelsmap")
    if (!is.logical(fr)) 
      stop("fr has to be TRUE or FALSE")
    if (!is.logical(um)) 
      stop("um has to be TRUE or FALSE")
    if (length(intersect(g, rownames(object@ndata))) < length(unique(g))) 
      stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
    if (!is.numeric(logsc) & !is.logical(logsc)) 
      stop("argument logsc has to be logical (TRUE/FALSE)")
    if (!is.null(cells)) {
      if (sum(!cells %in% colnames(object@ndata)) > 0) 
        stop("cells has to be a subset of cell ids, i.e. column names of slot ndata")
    }
    if (fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0) {
      if (dim(object@fr)[1] != 0) {
        fr <- TRUE
      }
      else if (dim(object@umap)[1] != 0) {
        um <- TRUE
      }
    }
    if (imputed & length(object@imputed) == 0) 
      stop("imputing needs to be done by running compdist with knn > 0")
    if (is.null(n)) 
      n <- g[1]
    if (is.null(cells)) 
      cells <- colnames(object@ndata)
    knn <- object@imputed$knn
    if (!noise) {
      if (length(g) == 1) {
        l <- as.vector(object@ndata[g, ] * min(object@counts) + 
                         0.1)
      }
      else {
        l <- apply(as.data.frame(as.matrix(object@ndata)[g, 
                                                         ]) * min(object@counts), 2, sum) + 0.1
      }
      if (imputed) {
        l <- apply(rbind(object@imputed$nn, object@imputed$probs), 
                   2, function(y) {
                     ind <- y[1:(knn + 1)]
                     p <- y[(knn + 2):(2 * knn + 2)]
                     sum(l[ind] * p)
                   })
      }
    }
    else {
      if (is.null(object@noise)) 
        stop("run noise analysis first!")
      if (length(g) == 1) {
        l <- as.vector(object@noise[g, ] + 0.1)
      }
      else {
        l <- apply(as.data.frame(as.matrix(object@noise)[g, 
                                                         ]), 2, sum) + 0.1
      }
    }
    if (logsc) {
      f <- l == 0
      l <- log2(l)
      l[f] <- NA
    }
    h <- colnames(object@ndata) %in% cells
  }
  mi <- min(l, na.rm = TRUE)
  ma <- max(l, na.rm = TRUE)
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  ColorLevels <- seq(mi, ma, length = length(ColorRamp))
  v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
  if ( seurat) {
    if (reduction=="umap") {d <- object@reductions$umap@cell.embeddings}
    if (reduction=="pca") { d <- object@reductions$pca@cell.embeddings}
    if (reduction=="tsne") {d <- object@reductions$pca@cell.embeddings}
    
  }
  else {
    if (fr) {
      d <- object@fr
    }
    else if (um) {
      d <- object@umap
    }
    else {
      d <- object@tsne
    }
  }
  pardefault <- par()
  layout(matrix(data = c(1, 3, 2, 4), nrow = 2, ncol = 2), 
         widths = c(5, 1, 5, 1), heights = c(5, 1, 1, 1))
  par(mar = c(3, 5, 2.5, 2))
  if (!leg) 
    n <- NA
  plot(c(min(d[, 1]), max(d[, 1])), c(min(d[, 2]), max(d[, 
                                                         2])), xlab = NA, ylab = NA, main = n, pch = 20, cex = 0, 
       col = "lightgrey", axes = FALSE)
  if (map) {
    v <- v[h]
    d <- d[h, ]
    kk <- order(v, decreasing = F)
    points(d[kk, 1], d[kk, 2], col = ColorRamp[v[kk]], pch = 20, 
           cex = cex)
  }
  if (leg) {
    par(mar = c(10, 2.5, 2.5, 4))
    image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), 
                                 nrow = 1), col = ColorRamp, xlab = "", ylab = "", 
          xaxt = "n")
    layout(1)
    par(mar = pardefault$mar)
  }
}

## compmedoids seurat clusters
compmedoids2 <- function (objectassay, part) 
  
{
  m <- c()
  for (i in sort(unique(as.numeric(part)))) {
    f <- names(part)[as.numeric(part) == i]
    if (length(f) == 1) {
      m <- append(m, f)
    }
    else {
      
      
      g <- apply(as.matrix(objectassay@data[, as.numeric(part) == 
                                              i]) - as.vector(pam(t(objectassay@data[, as.numeric(part) == 
                                                                                       i]), 1)$medoids), 2, sum) == 0
      m <- append(m, names(part)[as.numeric(part) == i][g])
      
    }
  }
  m
}

## plotmap for Seurat 

plotmap2 <- function (object, final = TRUE, tp = 1, fr = FALSE, um = FALSE, 
                     cex = 0.5, seurat=F, medoids, reduction="umap") 
{ 
  
  if (seurat) {
    if (reduction=="umap") {d <- object@reductions$umap@cell.embeddings}
    if (reduction=="pca") { d <- object@reductions$pca@cell.embeddings}
    if (reduction=="tsne") {d <- object@reductions$tsne@cell.embeddings}
    part <- object$seurat_clusters
  }
  else {
    if (length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 
        0) 
      stop("run comptsne/compfr/compumap before plotlabelsmap")
    if (final & length(object@cpart) == 0) 
      stop("run findoutliers before plotmap")
    if (!final & length(object@cluster$kpart) == 0) 
      stop("run clustexp before plotmap")
    if (!is.numeric(tp) | (is.numeric(tp) & tp > 1 | tp < 0)) 
      stop("tp has to be a number between 0 and 1 (transparency)")
    if (!is.logical(fr)) 
      stop("fr has to be TRUE or FALSE")
    if (!is.logical(um)) 
      stop("um has to be TRUE or FALSE")
    if (fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0) {
      if (dim(object@fr)[1] != 0) {
        fr <- TRUE
      }
      else if (dim(object@umap)[1] != 0) {
        um <- TRUE
      }
    }
    part <- if (final) 
      object@cpart
    else object@cluster$kpart
    if (fr) {
      d <- object@fr
    }
    else if (um) {
      d <- object@umap
    }
    else {
      d <- object@tsne
    }
  }
  if (seurat) { fcol <- rainbow(length(unique(as.numeric(part))))}
  else {
    if (is.null(object@fcol)) { fcol <- rainbow(length(unique(as.numeric(part))))}
    else { fcol <- object@fcol}
  }
  row.names(d) <- names(part)
  plot(d, xlab = "", ylab = "", cex = 0, axes = FALSE)
  for (i in 1:max(as.numeric(part))) {
    if (sum(as.numeric(part) == i) > 0) 
      points(d[as.numeric(part) == i, 1], d[as.numeric(part) == i, 2], col = adjustcolor(fcol[i], 
                                                                                         tp), pch = 20, cex = cex)
  }
  for (i in 1:max(as.numeric(part))) {
    if (sum(as.numeric(part) == i) > 0) 
      points(d[medoids[i], 1], d[medoids[i], 
                                 2], col = adjustcolor(fcol[i], tp), pch = 20, 
             cex = 4)
    if (sum(as.numeric(part) == i) > 0) 
      points(d[medoids[i], 1], d[medoids[i], 
                                 2], col = adjustcolor("white", tp), pch = 20, 
             cex = 3)
    if (sum(as.numeric(part) == i) > 0) 
      text(d[medoids[i], 1], d[medoids[i], 
                               2], i, col = adjustcolor("black", tp), cex = 0.75, 
           font = 4)
  }
}

## plotsymbolmap for Seurat

plotsymbolsmap2 <- function (object, types, subset = NULL, samples_col = NULL, cex = 0.5, 
                             fr = FALSE, um = FALSE, leg = TRUE, map = TRUE, seurat=T, reduction="umap") 
{
  if ( seurat == F)  {
    if (length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 
        0) 
      stop("run comptsne/compfr/compumap before plotlabelsmap")
    if (!is.logical(fr)) 
      stop("fr has to be TRUE or FALSE")
    if (!is.logical(um)) 
      stop("um has to be TRUE or FALSE")
    if (fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0) {
      if (dim(object@fr)[1] != 0) {
        fr <- TRUE
      }
      else if (dim(object@umap)[1] != 0) {
        um <- TRUE
      }
    }
  }
  if (is.null(subset)) 
    subset <- unique(types)
  h <- sort(unique(types)) %in% subset
  if (!is.null(subset)) {
    fp <- rep(FALSE, length(types))
    fp[types %in% subset] <- TRUE
  }
  if (is.null(samples_col)) {
    samples_col <- rainbow(length(unique(types[fp])))
  }
  else {
    samples_col <- samples_col[h]
  }
  if (seurat) {
    if (reduction=="umap") {d <- object@reductions$umap@cell.embeddings}
    if (reduction=="pca") { d <- object@reductions$pca@cell.embeddings}
    if (reduction=="tsne") {d <- object@reductions$pca@cell.embeddings}
  }
  else {
    if (fr) {
      d <- object@fr
    }
    else if (um) {
      d <- object@umap
    }
    else {
      d <- object@tsne
    }
  }
  if (map) {
    plot(d, xlab = "", ylab = "", axes = FALSE, cex = cex, 
         pch = 20, col = "grey")
    for (i in 1:length(unique(types[fp]))) {
      f <- types == sort(unique(types[fp]))[i]
      points(d[f, 1], d[f, 2], col = samples_col[i], pch = 20, 
             cex = cex)
    }
  }
  else {
    plot(d, xlab = "", ylab = "", axes = FALSE, cex = 0, 
         pch = 20, col = "grey", xlim = c(min(d[, 1]), max(d[, 
                                                             1])), ylim = c(min(d[, 2]), max(d[, 2])))
  }
  if (leg) 
    legend("topleft", legend = sort(unique(types[fp])), col = samples_col, 
           pch = 20, cex = 0.75, bty = "n")
}

### cell cycle scoring for RaceID

Scoring <- function(object, features, name=NULL, nbin = 24, ctrl= 100, n=sc@genes, st=NULL, nd=NULL, null=NULL ) {
  if ( is.null(features)) {
    stop("Input feature list is missing")
  }
  set.seed(1337)
  fdata <- object@ndata * min(object@counts)
  features.old <- features
  features <- lapply(seq_along(features), function(x, y, i) {
    z <- intersect(x[[i]], n)
    if ( length(z) == 0) {
      stop(paste("no feature of list element named ", y[[i]], " is element of n", sep = ""))
    } 
    return(z)
  }, y = names(features), x=features)
  names(features) <- names(features.old)
  feature.length <- length(features)
  universe <- n
  mean.universe <- apply(fdata[n,], 1, mean )
  mean.universe <- mean.universe[order(mean.universe)]
  mean.universe.binned <- ggplot2::cut_number(mean.universe + rnorm(length(mean.universe))/1e+30, n = nbin, labels=F, right=F)
  names(mean.universe.binned) <- names(mean.universe)
  
  ctrl.use <- vector(mode = "list", length = feature.length)
  for (i in 1:feature.length) {
    features.use <- features[[i]]   
    for (j in 1:length(x = features.use)) { 
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(sample(mean.universe.binned[which( mean.universe.binned == 
                                                                                   as.numeric(mean.universe.binned[features.use[j]]))], size = ctrl, replace = FALSE))) ### samples for each feature gene 100 control genes from the same bin as feature gene  
    }
  }
  ctrl.use <- lapply(ctrl.use, unique)
  ctrl.scores <- matrix(numeric(length = 1L), length(ctrl.use), ncol(fdata))
  
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- apply(fdata[features.use, ], 2, mean)    #### means across s.control or m.control genes for each cell
  }  
  features.scores <- matrix(numeric(length = 1L), feature.length, ncol(fdata))
  for (i in 1:feature.length) {
    features.use <- features[[i]]
    data.use <- fdata[features.use, , drop = FALSE]
    features.scores[i, ] <- apply(data.use, 2, mean) ### mean of feature genes across cells 
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(features.scores.use) <- paste0(names(features), "_factor")
  features.scores.use <- as.data.frame(t(features.scores.use))
  rownames(features.scores.use) <- colnames(fdata)
  
  if (length(features) == 1) { 
    assignments <- apply( features.scores.use, 1, function(scores) {
      if ( scores  < 0 ) { return(nd)}
      else { return(st)}
    })
  } 
  
  if ( length(features) > 1){
    if ( is.null(null)) {
      stop("set null")
    }
    assignments <- apply(features.scores.use, 1, FUN = function(scores) {
      if (all(scores < 0)) {
        return(null)
      }
      else {
        if (length(which(scores == max(scores))) > 1) {
          return("Undecided")
        }
        else {
          return(c(st, nd)[which(scores == max(scores))])
        }
      }
    })
  }
  features.scores.use <- cbind(features.scores.use, Phase=assignments)
  
  return(features.scores.use)
  
}
