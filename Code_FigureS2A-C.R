### Code figure S2 A-C 
###R version 3.6.1 ### RaceID3 version v0.1.6

### custom functions to load
source("~/Desktop/Functions_for_code_upload/Functions_to_load.R")

### loading of uninfected data 
setwd("~/Desktop/ILC_Paper/")
wt_lung <- readRDS("~/Desktop/RawData_FigureReproduce/lung_full_blood.RData")
wt_lung <- wt_lung[,grep("wt", colnames(wt_lung))]

#loading of cells worked with in uninfected dataset 
sc_wt <- readRDS("~/Desktop/RaceIDobjects_mainFigures/sc_wtlung.RData")
# only take cells worked with in Figure1
wt_lung <- wt_lung[, colnames(wt_lung) %in% colnames(sc_wt@ndata)]
wt_lung <- Matrix(as.matrix(wt_lung), sparse = T)

### loading of lung and skin dataset of Ricardo-Gonzalez et al. 
x <- readMM("~/Desktop/RawData_FigureReproduce/GSM3303971_Adult_Lung_matrix.mtx")
x <- Matrix(as.matrix(x),sparse=TRUE)
colnames(x) <- paste("RicGonzLung", 1:ncol(x), sep="_")
genesx <- read.delim("~/Desktop/RawData_FigureReproduce/GSM3303971_Adult_Lung_genes.tsv", header = F)

y <- readMM("~/Desktop/RawData_FigureReproduce/GSM3303972_Adult_skin_matrix.mtx")
y <- Matrix(as.matrix(y),sparse=TRUE)
colnames(y) <- paste("RicGonzSkin", 1:ncol(y), sep="_")
genesy <- read.delim("~/Desktop/RawData_FigureReproduce/GSM3303972_Adult_skin_genes.tsv", header=F)

# test whether genes are the same for lung and skin dataset from Ricardo-Gonzalez et al. and combine datasets
identical(genesx, genesy)

#combine datasets
locks <- cbind(x, y)
rownames(locks) <- genesx$V2

#combine counts of EnsemblIDs which correspond to the same gene symbol
non_unique_genes <- names(table(rownames(locks))[which(table(rownames(locks)) > 1)])

for ( i in 1:length(non_unique_genes)) {
  cat(paste(dim(locks)[1], ",", sep = ""))
  indixes <- which(rownames(locks) == non_unique_genes[i])
  if ( length(indixes) == 2){
    vrow <- var(rowSums(locks[indixes,]))
    
    if (vrow < 0.5) { 
      cat("mean,")
      new <- apply(locks[indixes,], 2, mean)
      locks[indixes[1],] <- new }
    else {
      cat("sum,")
      new <- apply(locks[indixes,], 2, sum)
      locks[indixes[1],] <- new }
      }  
  else {
    new <- apply(locks[indixes,], 2, sum)
    locks[indixes[1],] <- new }
  
  locks <- locks[-indixes[2:length(indixes)],]
  cat(paste(dim(locks)[1], ",", sep = ""))
}

## take only the common genes between this study and Ricardo-Gonzalez et al.
inter <- intersect(rownames(wt_lung), rownames(locks))#[-1]
celseq <- wt_lung[inter,]
colnames(celseq) <- paste("Cellseq", colnames(celseq), sep = "_")
tenX <- locks[inter,]

## pre processing 
mintotal <- 2000 ### lower threshhold to include more cells from Ricardo-Gonzalez et al.
minfeat <- 200
getfeatures <- function(data) { apply(data > 0, 2, sum) } 

celseq <- celseq[,apply(celseq, 2, sum) >= 2000]
celseq <- celseq[,getfeatures(celseq) >= 200]

tenX <- tenX[,apply(tenX, 2, sum) >= 2000]
tenX <- tenX[,getfeatures(tenX) >= 200]

prdata <- cbind(celseq, tenX)
minexpr <- round(5*mintotal/3000,0)


### generic set of genes removed from clustering stored in CGenes 
Gm            <- rownames(prdata)[grepl("Gm\\d+", rownames(prdata))] 
RP            <- rownames(prdata)[grepl("RP", rownames(prdata))] 
Hsp           <- rownames(prdata)[grepl("Hsp", rownames(prdata))] 
Igl            <- rownames(prdata)[grepl("Igl", rownames(prdata))] 
Igh           <- rownames(prdata)[grepl("Igh", rownames(prdata))] 
A4300  <- rownames(prdata)[grepl("A4300", rownames(prdata))]
CGenes        <- c( Gm, RP,  Hsp, "Scgb1a1", "Jchain", Igl, Igh, "Igkc", "Malat1", "Xist", A4300, "Mid1", "Kcnq1ot1")

### retrieve cell cycle genes from Buettner et al. and correlating genes 
cell_cycle_genes_scLVM_paper <- as.character(read.delim("~/Desktop/RawData_FigureReproduce/Cell_cycle_genes.txt", header = F)[,1])
gene.df <- matrix(ncol=2, nrow=length(cell_cycle_genes_scLVM_paper))
gene.df <- bitr(cell_cycle_genes_scLVM_paper, fromType = "ENSEMBL", toType = "SYMBOL",OrgDb = org.Mm.eg.db, drop=T) ### to entrez gene list
genes_wo_cell_cycle <- unique(gene.df$SYMBOL)

# filter genes as in RaceID before deriving correlating genes
filt_gene <- as.logical(apply(prdata, 1, max) >= 5)
filt_gene <- rownames(prdata)[filt_gene]

# correlating filtered genes to cell cycle genes on library size normalized counts and combine with CGenes
prdata2 <- t(t(prdata)/colSums(prdata))  
ccor <- cor(t(as.matrix(prdata2[genes_wo_cell_cycle[genes_wo_cell_cycle %in% rownames(prdata2)],])), t(as.matrix(prdata2[filt_gene[!filt_gene %in% c(CGenes,genes_wo_cell_cycle)],])))
cc_cor <- as.logical(apply(ccor, 2, max, na.rm=T) >= 0.4)
cc_cor_genes <- colnames(ccor)[cc_cor]
CGenes <- c(CGenes,  cc_cor_genes, genes_wo_cell_cycle)

#Run RaceID filtering for final filtered gene set 
ccor <- 0.65
FGenes <- NULL
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")


### seurat worflow, dataset integration based on sc transform 

# seurat object initialisation and sc transformation
celseq <- CreateSeuratObject(counts = celseq, project = "celseq", min.cells = 1, min.features = 200)
tenX <- CreateSeuratObject(counts = tenX, project = "tenX", min.cells = 1, min.features = 200)
combined.list <- list(celseq=celseq, tenX=tenX)
for (i in 1:length(combined.list)) {
  combined.list[[i]]  <- SCTransform(combined.list[[i]], verbose = FALSE)
}

# select anchoring features and filter out genes which are not element of CGenes or proliferating genes 
restinglung.features2 <- SelectIntegrationFeatures(object.list = combined.list, nfeatures = 3000)
restinglung.features2 <- restinglung.features2[restinglung.features2 %in% sc@genes]

#find anchors and perform integration 
options(future.globals.maxSize = 2000 * 1024^2)
combined.list <- PrepSCTIntegration(object.list = combined.list, anchor.features = restinglung.features2, 
                                    verbose = FALSE)
restinglung.anchors2 <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", 
                                               anchor.features = restinglung.features2, verbose = FALSE)
restinglung.integrated2 <- IntegrateData(anchorset = restinglung.anchors2, normalization.method = "SCT", 
                                         verbose = FALSE)
#perform clustering 
restinglung.integrated2 <- RunPCA(restinglung.integrated2, npcs = 30, verbose = FALSE)
restinglung.integrated2 <- FindNeighbors(restinglung.integrated2, dims = 1:30)
restinglung.integrated2 <- FindClusters(restinglung.integrated2, resolution = 1.2)
restinglung.integrated2 <- RunUMAP(restinglung.integrated2, dims = 1:30)

#plot UMAP of cluster identity
medoids <- compmedoids2(restinglung.integrated2, restinglung.integrated2$seurat_clusters)
plotmap2(restinglung.integrated2, seurat = T, medoids = medoids, um = T)

#plot UMAP of sample identity
restinglung.integrated2@meta.data$symbols <- sub("\\_.+", "", colnames(restinglung.integrated2@assays$integrated@data))
plotsymbolsmap2(restinglung.integrated2, types = sub("\\_.+", "", names(restinglung.integrated2$symbols)), leg=F)
plotsymbolsmap2(restinglung.integrated2, types = sub("\\_.+", "", names(restinglung.integrated2$symbols)), map=F)

#cluster composition 
order_comp <- c(1, 2, 20, 5, 13, 7, 15, 6, 10, 16, 9, 14, 17, 19, 12, 4, 3, 8, 11, 18, 21 )
get_cluster_composition(restinglung.integrated2, norm = F, seurat = T, cluster_size = 10, leg = F, order = order_comp)
get_cluster_composition(restinglung.integrated2, norm = F, seurat = T, cluster_size = 10, map = F, order = order_comp)

#fraction dotplot of candidate cells in skin cells of integrated cluster 
restinglung.integrated2 <- getRC_object(restinglung.integrated2) #normalize rawcounts based on libary size
cluster <- c(14, 17, 19, 12, 4, 3, 8, 11)
genes2 <- c( "Tbx21","Rorc", "Gata3", "Il1rl1", "Il17rb","Il2ra", "Il18r1", "Icos", "Thy1", "Il7r", "Zbtb16", "Tcf7", "Cd7", "Igf1r","Itgae", "Cd3e",  "S100a6","S1pr1","Vim","Klf2","Cd69","Mki67", "Id2","Kit","Arg1","Klrg1","Cxcl2", "Csf2", "Il2", "Calca","Il5","Il13","Pdcd1", "Ctla4","Il1r2","H2-Ab1", "Areg", "Klf4", "Rgs2")
fracdotplot(restinglung.integrated2, genes2, cluster = cluster, sampclus = "RicGonzSkin", limup = 1, limdown = -1, mintotal = 2000, seurat = T)

#plot UMAP indicating log2 normalized expression
plotexpmap2(restinglung.integrated2, "Zbtb16", logsc = T, seurat = T, mintotal = 2000)


