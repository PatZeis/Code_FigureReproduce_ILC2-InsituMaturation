### Code_FigureS4H-I 
###R version 3.6.1 ### RaceID3 version v0.1.6

### custom functions to load
source("~/Desktop/Functions_for_code_upload/Functions_to_load.R")

### loading and preprocessing 
prdata <- readRDS("~/Desktop/RawData_FigureReproduce/culture1.RData")
prdata2 <- readRDS("~/Desktop/RawData_FigureReproduce/culture2.RData")
prdata3 <- readRDS("~/Desktop/RawData_FigureReproduce/culture3.RData")
prdata <- cbind(prdata, prdata2, prdata3)
prdata <- prdata[,grep("LDL1|Progenit|BM|LinPos", colnames(prdata), invert = T)]

getype2 <- function(ma,y, y1){
  z <- sum(grepl(y, colnames(ma)))
  colnames(ma)[grep( y, colnames(ma))] <- paste(y1, 1:z, sep="_")
  return(ma)
}
prdata <- getype2(prdata, "LOP9IL33rnegIl18|LOP9IL33rnegIl18P5P6", "ST2negIl18pos")
prdata <- getype2(prdata, "LOP9IL33rnegIL25rIL33posIL25pos", "ST2negIl17rbnegIL33posIL25pos")
prdata <- getype2(prdata, "LOP9IL33rneg|LOP9IL33rnegP5P6", "ST2negIl18neg")
prdata <- getype2(prdata, "LungOP9Il18negP8", "ST2negIl17rbneg")
prdata <- getype2(prdata, "lungLOP9IL18neg", "PanIl18neg") # Pan=Il18r1+Icos+
prdata <- getype2(prdata, "lungLOP9IL18pos", "PanIl18pos") # Pan=Il18r1+Icos+

prdata <- prdata[grep("^(ERCC|mt)",row.names(prdata),invert=TRUE),]
cs <- apply(prdata,2,sum)
prdata <- prdata[,cs > 500]
cs <- cs[cs > 500]
f <- t(prdata["Kcnq1ot1",])/cs < .02
prdata <- prdata[,f]



### genes removed from clustering stored in CGenes 
Gm            <- rownames(prdata)[grepl("Gm\\d+", rownames(prdata))] 
RP            <- rownames(prdata)[grepl("RP", rownames(prdata))] 
Hsp           <- rownames(prdata)[grepl("Hsp", rownames(prdata))]
Igl            <- rownames(prdata)[grepl("Igl", rownames(prdata))] 
Igh           <- rownames(prdata)[grepl("Igh", rownames(prdata))] 
A4300  <- rownames(prdata)[grepl("A4300", rownames(prdata))]
CGenes        <- c( Gm, RP,  Hsp, "Scgb1a1", "Jchain", Igl, Igh, "Igkc", "Malat1", "Xist", A4300, "Mid1", "Kcnq1ot1", "Gzmb", "Ucp2", "Calca") 

### RaceID parameters
mintotal <- 3000
minexpr <- round(5*mintotal/3000,0)
outminc <- round(5*mintotal/3000,0)
ccor <- 0.65
FGenes <- NULL 

### Run RaceID functions
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=0.65, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=30)
sc <- compumap(sc)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

### remove Kcnq1ot1 cluster and rerun RaceID
out <- names(sc@cpart)[sc@cpart == 7] 
prdata <- prdata[, !colnames(prdata) %in% out]
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=0.65, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=30)
sc <- compumap(sc)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

### generation of plots

# plot UMAP for status of Il18 stimulation  
il18 <- colnames(sc@ndata)
il18[grep("Il18pos", il18, invert = T)] <- "Il18neg"
il18[grep("Il18pos", il18)] <- "Il18pos"
marker_col <- brewer.pal(3, "Set1")
marker_col <- marker_col[c(2,1)]
plotsymbolsmap(sc, types = il18, um = T, cex = 1.25, leg = T, samples_col = marker_col)

# plot UMAP for sample information: Pan=Il18r1+Icos+
col <- rainbow(length(unique(sub("\\_.+", "", colnames(sc@ndata)))))
col <- col[c(2, 3,1, 4, 5, 6)]
plotsymbolsmap(sc, types = sub("\\_.+", "", colnames(sc@ndata)), um = T, cex = 1.25, leg = T, samples_col = col)

# fractiondotplot
genes <- c("Gata3", "Rorc", "Tbx21", "Il1rl1", "Il17rb", "Il18r1", "Zbtb16", "Tcf7", "Mki67", "Id2", "Klrg1", "Il13", "Il5", "Areg", "Cxcl2", "Furin", "Calca", "Cxcr3", "Gzmc")
samples <- unique(sub("\\_.+", "", colnames(sc@ndata)))[c(1, 2, 4, 5, 3, 6)]
fracdotplot(sc, genes, samples=samples, limup = 1, limdown = -1)
