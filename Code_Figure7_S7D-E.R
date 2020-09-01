###Code_Figure7_S7D-E
###R version 3.5.1 ### RaceID3 version v0.1.3 

###custom functions to load 
source("~/Desktop/Functions_for_code_upload/Functions_to_load.R")

### loading and preprocessing 
d15chim2 <- readRDS("~/Desktop/RawData_FigureReproduce/chimera2.RData")
d15chim1 <- readRDS("~/Desktop/RawData_FigureReproduce/chimera1.RData")
prdata <- cbind(d15chim1, d15chim2)
prdata <- prdata[grep("^(ERCC|mt)",row.names(prdata),invert=TRUE),]
cs <- apply(prdata,2,sum)
prdata <- prdata[,cs > 500]
cs <- cs[cs > 500]
f <- t(prdata["Kcnq1ot1",])/cs < .02
prdata <- prdata[,f]


### generic set of genes removed from clustering stored in CGenes 
Gm            <- rownames(prdata)[grepl("Gm\\d+", rownames(prdata))] 
RP            <- rownames(prdata)[grepl("RP", rownames(prdata))] 
Hsp           <- rownames(prdata)[grepl("Hsp", rownames(prdata))] 
Igl            <- rownames(prdata)[grepl("Igl", rownames(prdata))] 
Igh           <- rownames(prdata)[grepl("Igh", rownames(prdata))] 
A4300  <- rownames(prdata)[grepl("A4300", rownames(prdata))]
CGenes        <- c( Gm, RP,  Hsp, "Scgb1a1", "Jchain", Igl, Igh, "Igkc", "Malat1", "Xist", A4300, "Mid1", "Kcnq1ot1")

### RaceID parameters
mintotal <- 3000
minexpr <- round(5*mintotal/3000,0)
outminc <- round(5*mintotal/3000,0)
ccor <- 0.65
FGenes <- NULL 

### Run RaceID
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=0.65, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=200)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

###removal of contaminant and low quality cells
out <- unique(c(names(sc@cpart)[sc@cpart %in% c(8, 20)], names(sc@cluster$kpart)[sc@cluster$kpart %in% c(7)]))
prdata <- prdata[, !colnames(prdata) %in% out] 

sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=0.65, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=200)
sc <- findoutliers(sc,probthr=0.0001,outminc=outminc,outlg=2,outdistquant=0.95)

#StemID lineage
ltree     <- list()
ltr <- Ltree(sc)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr=9,nmode=T,fr=FALSE) ##knn removed
ltr <- projback(ltr,pdishuf=2000,fast=F,rseed=17000)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.01)

### generating plots
#tNSE map cluster identity
plotmap(sc)
#tSNE map donor or host cells
getype2 <- function(ma,y, y1){
  z <- sum(grepl(y, ma))
  ma[grep( y, ma)] <- paste(y1, 1:z, sep="_")
  return(ma)
}
host_don <- getype2(names(sc@cpart), "host", "Host")
host_don <- getype2(host_don, "don", "Donor")
cpart <- sc@cpart
names(cpart) <- host_don
marker_col <- brewer.pal(3, "Set1")
coloc <- marker_col[1:2]
plotsymbolsmap(sc, types = sub("\\_.+", "", names(cpart)), samples_col = coloc, cex = 1)
#StemID lineage tree
plotgraph(ltr,showCells=FALSE,showTsne=TRUE,tp=.3,scthr=0)
#pie-charts
types <- sub("\\_.+","", names(cpart))
types_num <- as.numeric(table(types))
types <- names(table(types))
a <- as.numeric(names(table(sc@cpart)))
for ( i in a){
  x <- table(sub("\\_.+", "", names(cpart[cpart == i])))
  lbls <- names(x)
  y <- which(types %in% lbls)  
  slices <- as.numeric(x)/types_num[y] 
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(pct, " % ", lbls,sep=" ")
  pie(slices, labels=lbls,cex.main=1.5, cex = 1.25,col=coloc[y], cex.main=2,main=paste("Cluster ", i, " distribution ", " N=", length(cpart[cpart == i]), sep=" "))
}
#fracdotplot
genes2 <- c( "Tbx21","Rorc", "Gata3", "Il1rl1", "Il17rb","Il2ra", "Il18r1", "Icos", "Thy1", "Il7r", "Zbtb16", "Tcf7", "Cd7", "Igf1r","Itgae", "Cd3e",  "S100a6","S1pr1","Vim","Klf2","Cd69","Mki67", "Id2","Kit","Arg1","Klrg1","Cxcl2", "Csf2", "Il2", "Calca","Il5","Il13","Pdcd1", "Ctla4","Il1r2","H2-Ab1", "Areg", "Klf4", "Rgs2")
cluster2 <- c(16, 17,  6, 10, 11,13, 3, 2, 4, 5, 12, 8,7,14, 15,1,19, 9)
#zscore
fracdotplot(sc, genes = genes2,cluster = cluster2, limup = 1,limdown =  -1)
#log2mean
fracdotplot(sc, genes = genes2,cluster = cluster2, limup = 4,limdown =  -4, zsc = F)

### MA plots
ndata <- data.frame(as.matrix(sc@ndata * min(sc@counts)) + 0.1)
colnames(ndata) <- names(cpart)
#cluster4
clust <- names(cpart)[cpart == 4]
clust_host <- clust[grep("Host", clust)]
clust_don <- clust[grep("Don", clust)]
xd <- diffexpnb(ndata, clust_host, clust_don, vfit = sc@background$vfit)
xd$res <- xd$res[rownames(xd$res) %in% sc@genes,]
plotdiffgenesnb(xd, padj = F, Aname = "Host", Bname = "Donor")
#cluster7
clust <- names(cpart)[cpart == 7]
clust_host <- clust[grep("Host", clust)]
clust_don <- clust[grep("Don", clust)]
xd <- diffexpnb(ndata, clust_host, clust_don, vfit = sc@background$vfit)
xd$res <- xd$res[rownames(xd$res) %in% sc@genes,]
plotdiffgenesnb(xd, padj = F, Aname = "Host", Bname = "Donor")

###enrichment 
colnames(sc@ndata) <- names(cpart)
names(sc@cpart) <- names(cpart)
get_enrichment(sc, c("Donor", "Host"), clustsize = 9)

#enrichment cluster 16+17
pop <- length(sc@cpart)
size_cluster <- length(sc@cpart[sc@cpart %in% c(16,17)])
#don
size_samp <- sum(grepl("Donor", names(sc@cpart)))
size_samp_clust <- sum(grepl("Donor", names(sc@cpart[sc@cpart %in% c(16,17)])))
fisher <- fisher.test(matrix(c(pop - size_samp, size_samp,size_cluster - size_samp_clust, size_samp_clust),ncol=2), alternative = "greater")
don_pval <- fisher$p.value
#host
size_samp2 <- sum(grepl("Host", names(sc@cpart)))
size_samp_clust2 <- sum(grepl("Host", names(sc@cpart[sc@cpart %in% c(16,17)])))
fisher2 <- fisher.test(matrix(c(pop - size_samp2, size_samp2,size_cluster - size_samp_clust2, size_samp_clust2),ncol=2), alternative = "greater")
host_pval <- fisher2$p.value

#expected frequency 
rel_size <- types_num/pop
rel_size_abs <- rel_size *100

#enrichment analysis
x <- table(as.character(sapply(names(sc@cpart[sc@cpart %in% c(16,17)]), function(x) str_split(x, pattern="[_]")[[1]][1])))
slices <- as.numeric(x)/types_num
pct <- round(slices/sum(slices)*100)
enrichment <- pct/rel_size_abs[y]
names(enrichment) <- lbls
barplot(enrichment,  names.arg=c(round(-log10(don_pval), 2), round(-log10(host_pval),2)), col = marker_col, ylim = range(0:2))
