##Code_FigureS3F-H
###R version 3.6.1 ### RaceID3 version v0.1.6

### custom functions to load
source("~/Desktop/Functions_for_code_upload/Functions_to_load.R")

# load sce R file from Walker et al. 
sc_experiment <- readRDS("~/Desktop/RawData_FigureReproduce/GSE131038_ilc_R_sce_object_Figure7A.RData")
# load cluster partition of Walker et al.
cluster <- read.delim("~/Desktop/RawData_FigureReproduce/GSE131038_ilc_R_sce_object_clusters_Figure7A.tsv", header = F)

# filter cells with expression of non ILCs in cluster 1 and 2 and filter genes 
raw_data_cl1_cl2 <- data.frame(sc_experiment@assays$`.->data`$counts)[,cluster$V2 %in% c(1,2)]
raw <- raw_data_cl1_cl2
gene.df <- matrix(ncol=2, nrow=nrow(raw))
genes_all <- rownames(raw)
genes_all2 <- sub("\\..+", "", genes_all )

gene.df <- bitr(genes_all2, fromType = "ENSEMBL", toType = c("SYMBOL"),OrgDb = org.Mm.eg.db, drop=F)
nas <- which(is.na(gene.df$SYMBOL)) ## Enrez ID correspond to Gms don't have a symbol
gene.df <- gene.df[-nas,] ### genes without Gms 
genes_ensemblID <- genes_all[genes_all2 %in% gene.df$ENSEMBL]

give_entrez <- function(gene, data) {
  x <- gene.df[gene.df$SYMBOL == gene,]
  y <- x$ENSEMBL[!is.na(x$ENSEMBL)]
  z <- rownames(data)[grep(y, rownames(data))]
  return(z)
}

non_ILC_genes <- c("Cd302", "Sox13", "Cd19", "Spi1", "Gata1")
out <- c()
for ( i in non_ILC_genes) {
  out <- append(out, colnames(raw)[which(raw[give_entrez(i, raw),] > 0)])
}
out <- unique(out)
cells <- as.character(cluster$V1[cluster$V2 %in% c(1,2)])
cells <- cells[!cells %in% out ]

# subset normalized and logtransformed count table and match gene symbols
walker_data <- sc_experiment@assays$`.->data`$logcounts
walker_data <- walker_data[genes_ensemblID,cells]

xx <- sapply(genes_ensemblID, function(x){
  zz <- sub("\\..+", "", x)
  symbol <- gene.df[gene.df$ENSEMBL == zz,"SYMBOL"]
  return(symbol)
})
xx <- as.character(xx)
rownames(walker_data) <- xx

# combine counts of EnsemblIDs which correspond to the same gene symbol
non_unique_genes <- names(table(xx)[which(table(xx) > 1)])

for ( i in 1:length(non_unique_genes)) {
  cat(paste(dim(walker_data)[1], ",", sep = ""))
  indixes <- which(rownames(walker_data) == non_unique_genes[i])
  if ( length(indixes) == 2){
    vrow <- var(rowSums(walker_data[indixes,]))
    
    if (vrow < 0.5) { 
      cat("mean,")
      new <- apply(walker_data[indixes,], 2, mean)
      walker_data[indixes[1],] <- new }
    else {
      cat("sum,")
      new <- apply(walker_data[indixes,], 2, sum)
      walker_data[indixes[1],] <- new }
    
    
  }  
  else {
    new <- apply(walker_data[indixes,], 2, sum)
    walker_data[indixes[1],] <- new }
  
  walker_data <- walker_data[-indixes[2:length(indixes)],]
  cat(paste(dim(walker_data)[1], ",", sep = ""))
}

# combine with Il18r1+Icos+ BM data 
sc_bm <- readRDS("~/Desktop/RaceIDobjects_mainFigures/sc_BM.RData")
prdata <- as.data.frame(as.matrix(sc_bm@expdata[,colnames(sc_bm@ndata)]))
inter <- intersect(rownames(prdata), rownames(walker_data))
combined <- cbind(prdata[inter, ], walker_data[inter,] )

# RaceID parameter and CGenes initialization 
ccor <- 0.65
FGenes <- NULL
mintotal <- 1500
minexpr <- round(5*mintotal/3000,0)

Gm            <- rownames(prdata)[grepl("Gm\\d+", rownames(prdata))] 
RP            <- rownames(prdata)[grepl("RP", rownames(prdata))] 
Hsp           <- rownames(prdata)[grepl("Hsp", rownames(prdata))] 
Igl            <- rownames(prdata)[grepl("Igl", rownames(prdata))] 
Igh           <- rownames(prdata)[grepl("Igh", rownames(prdata))] 
A4300  <- rownames(prdata)[grepl("A4300", rownames(prdata))]
CGenes        <- c( Gm, RP,  Hsp, "Scgb1a1", "Jchain", Igl, Igh, "Igkc", "Malat1", "Xist", A4300, "Mid1", "Kcnq1ot1")

# Run RaceID for object initialisation and gene/cell filtering 

sc <- SCseq(combined)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")

### Seurat standart integration 
celseq <- CreateSeuratObject(counts = prdata[inter,], project = "celseq", min.cells = 1, min.features = 200)
smartseq <- CreateSeuratObject(counts=walker_data[inter,], project = "smartseq", min.cells = 1, min.features = 200)
combined.list <- list(celseq=celseq, smartseq=smartseq)

# normalize and select variable features
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- NormalizeData(combined.list[[i]], verbose = FALSE)
  combined.list[[i]] <- FindVariableFeatures(combined.list[[i]], selection.method = "vst",  nfeatures = 3000, verbose = FALSE)
}

# filter out genes from variable features which are not element of CGenes  
combined.list$celseq@assays$RNA@var.features <- combined.list$celseq@assays$RNA@var.features[combined.list$celseq@assays$RNA@var.features %in% sc@genes]
combined.list$smartseq@assays$RNA@var.features <- combined.list$smartseq@assays$RNA@var.features[combined.list$smartseq@assays$RNA@var.features %in% sc@genes]

# anchoring and integration
restinglung.anchors <- FindIntegrationAnchors(object.list = combined.list, dims = 1:30)
restinglung.integrated <- IntegrateData(anchorset = restinglung.anchors, dims = 1:30)
dimred <- data.frame(as.matrix(restinglung.integrated@assays$integrated@data))

### Run RaceID

sc@dimRed$x <- dimred
sc <- compdist(sc,metric="pearson",FSelect=F,knn=NULL) ### no feature selection as background model might not be applicable due to dataset integration
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=5,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids") ##initially 4 - for previous analysis
sc <- comptsne(sc,rseed=15555,perplexity=30)
sc <- compumap(sc)
sc@medoids <- compmedoids(sc, sc@cluster$kpart)

# assign cluster partition from initial cluster analysis (Figure 2g for Il18r1+Icos+ or Walker et. al cl1=ILCPs or cl2= ILC2Ps)
bm <- colnames(sc@ndata)[grep("BM", colnames(sc@ndata))]
bm_part <- sc_bm@cpart 
bm_part_2 <- paste("Il18r1+Icos+cl", bm_part, sep = "")
names(bm_part_2) <- names(sc_bm@cpart)

mac <- colnames(sc@ndata)[grep("BM", colnames(sc@ndata), invert = T)]
zz <- sapply(mac, function(x) {
  z <- cluster[cluster$V1 == x,"V2"]
})
zz2 <- paste("Walker.et.al.cl", zz, sep = "")
names(zz2) <- names(zz)

symbol <- c(bm_part_2, zz2)
symbol2 <- sapply(colnames(sc@ndata), function(x){
  z <- as.character(symbol[names(symbol) == x])
})

library(RColorBrewer)
set1 <- brewer.pal(9, "Set1")
dark <- brewer.pal(6, "Dark2")
set1[6] <- dark[6]
set2 <- brewer.pal(8, "Set2")
colloc <- c(set1, set2[1])

# plot UMAP of cluster identity 
plotmap(sc, final = F, um = T, cex = 1.25)

# plot UMAP of sample idendity as cluster partition of initial clustering
plotsymbolsmap(sc, types = as.character(symbol2), cex = 1.25, um = T, samples_col = colloc)

# cluster composition
get_cluster_composition(sc, cluster_size = 10, norm = F, symbols = as.character(symbol2) , color = colloc, leg = F, order = c(5,1,2,4,3)) ## old order  order = c(4,1,3,2)
get_cluster_composition(sc, cluster_size = 10, norm = F, symbols = as.character(symbol2) , color = colloc, map = F)

#fraction dotplot of candidate cells integrated clusters
cluster <- c(5,1,2,4,3)
genes2 <- c( "Tbx21","Rorc", "Gata3", "Il1rl1", "Il17rb","Il2ra", "Il18r1", "Icos", "Thy1", "Il7r", "Zbtb16", "Tcf7", "Cd7", "Igf1r","Itgae", "Cd3e",  "S100a6","S1pr1","Vim","Klf2","Cd69","Mki67", "Id2","Kit","Arg1","Klrg1","Cxcl2", "Csf2", "Il2", "Calca","Il5","Il13","Pdcd1", "Ctla4","Il1r2","H2-Ab1", "Areg", "Klf4", "Rgs2")
fracdotplot(sc, genes2, cluster = cluster, limup = 1, limdown = -1)
