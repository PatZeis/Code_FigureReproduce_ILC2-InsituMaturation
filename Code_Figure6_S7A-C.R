###Code_Figure6_S7A-C
###R version 3.5.1 ### RaceID3 version v0.1.3 

### custom functions to load
source("~/Desktop/Functions_for_code_upload/Functions_to_load.R")

### loading and preprocessing 
para1 <- readRDS("~/Desktop/RawData_FigureReproduce/Parabiosis1.RData")
para2 <- readRDS("~/Desktop/RawData_FigureReproduce/Parabiosis2.RData")

colnames(para1) <- paste("Sort1", colnames(para1), sep="")
colnames(para2) <- paste("Sort2", colnames(para2), sep="")
prdata <- cbind(para1,para2)
prdata <- prdata[grep("^(ERCC|mt)",row.names(prdata),invert=TRUE),]
cs <- apply(prdata,2,sum)
prdata <- prdata[,cs > 500]
cs <- cs[cs > 500]
f <- t(prdata["Kcnq1ot1",])/cs < .02
prdata <- prdata[,f]

### generic set of genes removed from clustering stored in CGenes + cell cycle and low quality associated genes
Gm            <- rownames(prdata)[grepl("Gm\\d+", rownames(prdata))] 
RP            <- rownames(prdata)[grepl("RP", rownames(prdata))] 
Hsp           <- rownames(prdata)[grepl("Hsp", rownames(prdata))] 
Igl            <- rownames(prdata)[grepl("Igl", rownames(prdata))] 
Igh           <- rownames(prdata)[grepl("Igh", rownames(prdata))] 
A4300  <- rownames(prdata)[grepl("A4300", rownames(prdata))]
CGenes        <- c( Gm, RP,  Hsp, "Scgb1a1", "Jchain", Igl, Igh, "Igkc", "Malat1", "Xist", A4300, "Mid1", "Kcnq1ot1", "Lars2", "Mki67", "Pcna") 

### RaceID parameters
mintotal <- 2000
minexpr <- round(5*mintotal/3000,0)
outminc <- round(5*mintotal/3000,0)
ccor <- 0.4 
FGenes <- NULL

### Run RaceID
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=200)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

###removal of contaminant clusters
out <- unique(c(names(sc@cpart)[sc@cpart %in% c(4,5,7,8)], names(sc@cluster$kpart)[sc@cluster$kpart %in% c(2,3)]))
prdata <- prdata[, !colnames(prdata) %in% out] 
###re-run RaceID
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=200)
sc <- findoutliers(sc,probthr=0.0001,outminc=outminc,outlg=2,outdistquant=0.95)

### StemID lineage inference 
ltree     <- list()
ltr <- Ltree(sc)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr=20,nmode=T,fr=FALSE) ##knn removed
ltr <- projback(ltr,pdishuf=2000,fast=F,rseed=17000)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.01)

### renaming 
getype2 <- function(ma,y, y1){
  z <- sum(grepl(y, ma))
  ma[grep( y, ma)] <- paste(y1, 1:z, sep="_")
  return(ma)
}
host_don <- getype2(names(sc@cpart), "host", "Host")
host_don <- getype2(host_don, "don", "Donor")
cpart <- sc@cpart
names(cpart) <- host_don

### generation of plots
#tNSE map cluster identity
plotmap(sc)
#StemID lineage tree
plotgraph(ltr,showCells=FALSE,showTsne=TRUE,tp=.3,scthr=0)

###pie-charts
marker_col <- brewer.pal(3, "Set1")
coloc <- marker_col[1:2]
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

### MA plots

ndata <- data.frame(as.matrix(sc@ndata * min(sc@counts)) + 0.1)
colnames(ndata) <- names(cpart)
#cluster2
clust <- names(cpart)[cpart == 2]
clust_host <- clust[grep("Host", clust)]
clust_don <- clust[grep("Don", clust)]
xd <- diffexpnb(ndata, clust_host, clust_don, vfit = sc@background$vfit)
xd$res <- xd$res[rownames(xd$res) %in% sc@genes,]
plotdiffgenesnb(xd, padj = F, Aname = "Host", Bname = "Donor")
#cluster9
clust <- names(cpart)[cpart == 9]
clust_host <- clust[grep("Host", clust)]
clust_don <- clust[grep("Don", clust)]
xd <- diffexpnb(ndata, clust_host, clust_don, vfit = sc@background$vfit)
xd$res <- xd$res[rownames(xd$res) %in% sc@genes,]
plotdiffgenesnb(xd, padj = F, Aname = "Host", Bname = "Donor")

# tSNE map of marker gene 
plotexpmap(sc, "Gata3", logsc = T)
#fractiondotplot
cluster <- c(11, 4, 17, 3, 5, 14, 10, 18, 8, 9, 2, 1, 6, 7, 12, 13, 15)
genes2 <- c( "Tbx21","Rorc", "Gata3", "Il1rl1", "Il17rb","Il2ra", "Il18r1", "Icos", "Thy1", "Il7r", "Zbtb16", "Tcf7", "Cd7", "Igf1r","Itgae", "Cd3e",  "S100a6","S1pr1","Vim","Klf2","Cd69","Mki67", "Id2","Kit","Arg1","Klrg1","Cxcl2", "Csf2", "Il2", "Calca","Il5","Il13","Pdcd1", "Ctla4","Il1r2","H2-Ab1", "Areg", "Klf4", "Rgs2")
fracdotplot(sc, genes = genes2,cluster = cluster, limup = 1,limdown =  -1, zsc=T)
fracdotplot(sc, genes = genes2,cluster = cluster, limup = 4,limdown =  -4, zsc=F)

#cluster mapping 
sc_time <- readRDS("~/Desktop/RaceIDobjects_mainFigures/sc_timecourse.RData")
medoids_parabiosis <-  as.matrix(sc@ndata[,sc@medoids])

common_features <- Reduce(intersect, list(sc@cluster$features, sc_time@cluster$features))
common_genes <- Reduce(intersect, list(sc@genes, sc_time@genes))

QP <- function(k,m,norm=TRUE){
  library(quadprog)
  Dmat <- t(m) %*% m
  dvec <- t(k) %*% m
  if ( norm ){
    Amat <- cbind(rep(1,ncol(m)), diag(ncol(m)))
    bvec <- c(1,rep(0,ncol(m)))
    qp <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  }else{
    Amat <- diag(ncol(m))
    bvec <- rep(0,ncol(m))
    qp <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0, factorized=FALSE)
  }
  w <- qp$solution
  return( list(w=w, fit=m %*% w, residual= sum((m %*% qp$solution - k)**2), qp=qp)) 
}
comp <- apply(as.matrix(sc_time@ndata[common_features,]),2, function(x) {QP(as.vector(x), medoids_parabiosis[common_features,])})
comp_w1 <- lapply(comp, function(x) { x$w } )
df_w1 <- as.data.frame(comp_w1)
rownames(df_w1) <- paste("Weights parabiois cluster", rownames(df_w1))
scqp  <- sc_time
scqp@ndata <- df_w1

#plot weights 
plotexptsne(scqp, "Weights parabiois cluster 4")
plotexptsne(scqp, "Weights parabiois cluster 11")
            
### enrichment scores
colnames(sc@ndata) <- names(cpart)
names(sc@cpart) <- names(cpart)
get_enrichment(sc, c("Donor", "Host"))
