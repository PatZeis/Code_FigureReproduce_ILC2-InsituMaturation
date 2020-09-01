## Code figure S1 h-k 
###R version 3.6.1 ### RaceID3 version v0.1.6

### custom functions to load
source("~/Desktop/Functions_for_code_upload/Functions_to_load.R")

### retrieve count matrix and cell cycle gene list from Buettner el. 
prdata <- readRDS("~/Desktop/RawData_FigureReproduce/lung_full_blood.RData")
prdata <- prdata[,grep("wt", colnames(prdata))]

### preprocessing 
mintotal <- 3000
minexpr <- round(5*mintotal/3000,0)

prdata <- prdata[grep("^(ERCC|mt)",row.names(prdata),invert=TRUE),] ### remove counts for mitochondrial genes and ERCC spike ins
cs <- apply(prdata,2,sum)  
prdata <- prdata[,cs >= mintotal]  ## only include cells with at least mintotal transcripts
f <- t(prdata["Kcnq1ot1",])/cs < .02  ## filter cells 
prdata <- prdata[,f]

### generic set of genes removed from clustering stored in CGenes 
Gm            <- rownames(prdata)[grepl("Gm\\d+", rownames(prdata))] 
RP            <- rownames(prdata)[grepl("RP", rownames(prdata))] 
Hsp           <- rownames(prdata)[grepl("Hsp", rownames(prdata))]
Igl            <- rownames(prdata)[grepl("Igl", rownames(prdata))] 
Igh           <- rownames(prdata)[grepl("Igh", rownames(prdata))] 
A4300  <- rownames(prdata)[grepl("A4300", rownames(prdata))]
CGenes        <- c( Gm, RP,  Hsp, "Scgb1a1", "Jchain", Igl, Igh, "Igkc", "Malat1", "Xist", A4300, "Mid1", "Kcnq1ot1")

### retrieve cell cycle genes from Buettner et al. and correlating genes 
cell_cycle_genes_scLVM_paper <- as.character(read.delim("~/Desktop/Cell_cycle/Cell_cycle_genes.txt", header = F)[,1])
gene.df <- matrix(ncol=2, nrow=length(cell_cycle_genes_scLVM_paper))
gene.df <- bitr(cell_cycle_genes_scLVM_paper, fromType = "ENSEMBL", toType = "SYMBOL",OrgDb = org.Mm.eg.db, drop=T) ### matching EntrezID with gene Symbol
genes_wo_cell_cycle <- unique(gene.df$SYMBOL) ### cell cycle genes as symbols

# filter genes as in RaceID before deriving correlating genes
filt_gene <- as.logical(apply(prdata, 1, max) >= 5)
filt_gene <- rownames(prdata)[filt_gene]

# correlating filtered genes to cell cycle genes on library size normalized counts and combine with CGenes
prdata2 <- t(t(prdata)/colSums(prdata)) 
ccor <- cor(t(as.matrix(prdata2[genes_wo_cell_cycle[genes_wo_cell_cycle %in% rownames(prdata2)],])), t(as.matrix(prdata2[filt_gene[!filt_gene %in% c(CGenes,genes_wo_cell_cycle)],])))  
cc_cor <- as.logical(apply(ccor, 2, max, na.rm=T) >= 0.4) 
cc_cor_genes <- colnames(ccor)[cc_cor] 
CGenes <- c(CGenes,  cc_cor_genes, genes_wo_cell_cycle)

### further RaceID parameters
outminc <- round(5*mintotal/3000,0)
ccor <- 0.65
FGenes <- NULL 

### Run RaceID functions 
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=0.65, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=200)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

## filterout ILC1s/NKcells and rerun RaceID functions
cl6 <- names(sc@cpart)[sc@cpart == 6]
prdata <- prdata[, !colnames(prdata) %in% cl6]
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=0.65, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=8,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=200)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)
sc <- compumap(sc)

## StemID lineage tree inference 
ltree     <- list()
ltr <- Ltree(sc)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr=20,nmode=T,fr=FALSE) 
ltr <- projback(ltr,pdishuf=2000,fast=F,rseed=17000)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.01)



### generation of plots 

#fractiondotplot
genes2 <- c( "Tbx21","Rorc", "Gata3", "Il1rl1", "Il17rb","Il2ra", "Il18r1", "Icos", "Thy1", "Il7r", "Zbtb16", "Tcf7", "Cd7", "Igf1r","Itgae", "Cd3e",  "S100a6","S1pr1","Vim","Klf2","Cd69","Mki67", "Id2", "Pdcd1","Kit","Arg1","Klrg1", "Cxcl2", "Csf2", "Il2", "Calca","Il5","Il13", "Areg", "Klf4", "Rgs2")
cluster <- c(1, 7, 9, 2, 5, 4, 3, 6, 12, 8)
fracdotplot(sc, genes2, cluster = cluster, limup = 1, limdown = -1)

#StemID lineagree tree
plotgraph(ltr,showCells=FALSE,tp=.3,scthr=0)

#pseudo temporal expression profile
# get pseudo temporal expression data and order normalised expression matrix based on this order
n <- cellsfromtree(ltr, c(1,7,9,2,5,8)) ### cxcl2 branch
ndata <- as.matrix(sc@ndata * min(sc@counts)) + 0.1
ndata <- data.frame(ndata)
ndata <- ndata[sc@genes,n$f]

fcol <- sc@fcol
y    <- sc@cpart[n$f]

genes_to_plot_up <- c("Gata3", "Il1rl1", "Arg1", "Klrg1", "Bcl11b", "Cxcl2")
genes_to_plot_down <- c("Il18r1", "Tcf7", "Zbtb16", "Cd7", "Ikzf2", "Maf")

#Plot pseudo-temporal gene expression profiles
plotexpression(ndata, y, genes_to_plot_up, n$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), legendpos = "topleft")
plotexpression(ndata, y, genes_to_plot_down, n$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), legendpos = "topright")

### cell cycle comparison 
# pre-requesites for cell cycle scoring 

first_letter_uper <- function(x) {
  y <- paste(toupper(str_split(x, pattern = "")[[1]][1]), paste(str_split(x, pattern = "")[[1]][-1], collapse = ""), sep = "")
  return(y)
}

s.genes <- cc.genes.updated.2019$s.genes
s.genes <- as.character(sapply(tolower(s.genes), function(x) first_letter_uper(x)))

g2m.genes <- cc.genes.updated.2019$g2m.genes
g2m.genes <- as.character(sapply(tolower(g2m.genes), function(x) first_letter_uper(x)))

features <- list(cell_cycle=c(s.genes, g2m.genes))

## scoring of partition with and without cell cycle removal 

# loading of RaceID3 clustered uninfected data without cell cycle removal 
sc_uninfected <- readRDS("~/Desktop/RaceIDobjects_mainFigures/sc_wtlung.RData")

# universe for sampling control genes 
uni_genes <- sc_uninfected@genes

#scoring 
cycl_non_cycl <- Scoring(sc_uninfected, features, st = "cycling", nd="noncycling", n=uni_genes)

cycl_non_cycl_cc_removal <- Scoring(sc, features, st = "cycling", nd="noncycling", n=uni_genes)

# visualisation of fraction cycling and non-cycling cells per cluster 

get_cluster_composition(sc_uninfected, cluster_size = 20, final = T, leg = T, order = c(1  ,5 ,11 , 9 , 8 , 4 ,10 , 12,3 , 6 ,  7  ,2) , symbols = cycl_non_cycl$Phase, norm = F)
get_cluster_composition(sc, cluster_size = 20, final = T, leg = T, order = c(1, 7, 9, 2, 5, 4, 3, 6, 12, 8), symbols = cycl_non_cycl_cc_removal$Phase, norm = F)
