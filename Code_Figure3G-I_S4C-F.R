### Code_Figure3G-I_S4E-F
###R version 3.5.1 ### RaceID3 version v0.1.3 

### custom functions to load
source("~/Desktop/Functions_for_code_upload/Functions_to_load.R")

### loading and preprocessing 
prdata <- readRDS("~/Desktop/RawData_FigureReproduce/culture1.RData")
prdata2 <- readRDS("~/Desktop/RawData_FigureReproduce/culture2.RData")
prdata3 <- readRDS("~/Desktop/RawData_FigureReproduce/culture3.RData")
prdata <- cbind(prdata, prdata2, prdata3)
prdata <- prdata[,grep("LDL1|BM|LinPos|IL18pos|Il18P|Il18\\_\\d+|IL25r", colnames(prdata), invert = T)]

### renaming
getype2 <- function(ma,y, y1){
  z <- sum(grepl(y, colnames(ma)))
  colnames(ma)[grep( y, colnames(ma))] <- paste(y1, 1:z, sep="_")
  return(ma)
}
prdata <- getype2(prdata, "lungLOP9IL18neg", "cPan")
prdata <- getype2(prdata, "LOP9IL33rneg|LOP9IL33rnegP5P6", "aST2neg")
prdata <- getype2(prdata, "LungOP9Il18negP8", "bST2Il17rbDN")
prdata <- getype2(prdata, "Progenit", "dInput")


### data loading and preprocessing
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
mintotal <- 1500
minexpr <- round(5*mintotal/3000,0)
outminc <- round(5*mintotal/3000,0)
ccor <- 0.65
FGenes <- NULL 

### pre-Run for investigating batch genes 
prdata2 <- prdata[,grep("Pan|ST2neg", colnames(prdata))]

Pan <- colnames(prdata2)[grep("Pan", colnames(prdata2))]
ST2neg <- colnames(prdata2)[grep("ST2neg", colnames(prdata2))]
#DN2 <- colnames(prdata)[grep("ST2Il17rbDN", colnames(prdata))]
sc <- SCseq(prdata2)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=list(Pan, ST2neg), knn=10, FGenes=FGenes, CGenes=CGenes, ccor=0.65, bmode="RaceID")

### Run RaceID
CGenes <- c( Gm, RP,  Hsp, "Scgb1a1", "Jchain", Igl, Igh, "Igkc", "Malat1", "Xist", A4300, "Mid1", "Kcnq1ot1", "Gzmb", "Ucp2", "Calca")
mintotal <- 3000
minexpr <- round(5*mintotal/3000,0)
outminc <- round(5*mintotal/3000,0)

sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=30)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

cl10 <- names(sc@cpart)[sc@cpart == 10] ##knq1ot1 cells
prdata <- prdata[,!colnames(prdata) %in% cl10]

sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=30)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

### generation of plots
#tSNE map cluster identity
plotmap(sc, final = F)
#tSNE map of sample identity
marker_col <- brewer.pal(5, "Set1")
marker_col <- marker_col[c(1,2, 5, 3)]
plotsymbolsmap(sc, types = sub("\\_.+", "", colnames(sc@ndata)), samples_col = marker_col, cex = 1.25, leg = F)
plotsymbolsmap(sc, types = sub("\\_.+", "", colnames(sc@ndata)), samples_col = marker_col, cex = 1.25, map = F)

#cluster and sample composition 
get_cluster_composition(sc, cluster_size = 10, norm = T, symbols = NULL , color = marker_col, leg = F, order = c(5,6,7,1,2,4,3,8))
get_cluster_composition(sc, cluster_size = 10, norm = T, symbols = NULL , color = marker_col, map = F, order = c(5,6,7,1,2,4,3,8))

get_cluster_composition(sc, cluster_size = 10, norm = F, symbols = NULL , color = sc@fcol, leg = F, sample = T, order = c("cPan", "aST2neg", "bST2Il17rbDN","dInput" ))
get_cluster_composition(sc, cluster_size = 10, norm = F, symbols = NULL , color = sc@fcol, map = F, sample = T, order = c("cPan", "aST2neg", "bST2Il17rbDN","dInput" ))

#fracdotplot cluster
genes <- c("Gata3", "Rorc", "Tbx21", "Il1rl1", "Il17rb", "Il18r1", "Zbtb16", "Tcf7", "Mki67", "Id2", "Klrg1", "Il13", "Il5", "Areg", "Cxcl2", "Calca", "Cxcr3", "Gzmc")
cluster <- c(5, 6, 7, 1, 2, 4, 3, 8)
fracdotplot(sc, genes, cluster = cluster, limup = 4, limdown = -4, zsc = F)
#fracdotplot sample S4E
fracdotplot(sc, genes, samples = c("cPan", "aST2neg", "bST2Il17rbDN", "dInput"), limup = 4, limdown = -4, zsc = F)


###combined analysis of uninfected lung (Figure 1C,F), and sorted Il18r1+Icos+ST2-Il17rb- ILCs used as input for the in vitro cultures 
prdata <- readRDS("~/Desktop/RawData_FigureReproduce/lung_full_blood.RData")
prdata <- prdata[,grep("wt", colnames(prdata))]
prdata2 <- readRDS("~/Desktop/RawData_FigureReproduce/culture2.RData")
prdata2 <- prdata2[, grep("Progenit", colnames(prdata2))]
prdata <- cbind(prdata, prdata2)

#pre-processing
prdata <- prdata[grep("^(ERCC|mt)",row.names(prdata),invert=TRUE),]
cs <- apply(prdata,2,sum)
prdata <- prdata[,cs > 500]
cs <- cs[cs > 500]
f <- t(prdata["Kcnq1ot1",])/cs < .02
prdata <- prdata[,f]

# generic set of genes removed from clustering stored in CGenes 
Gm            <- rownames(prdata)[grepl("Gm\\d+", rownames(prdata))] ### also FGenes also functional genes
RP            <- rownames(prdata)[grepl("RP", rownames(prdata))] ### also FGenes same as Rik
Hsp           <- rownames(prdata)[grepl("Hsp", rownames(prdata))] ## okay
Igl            <- rownames(prdata)[grepl("Igl", rownames(prdata))] ## okay
Igh           <- rownames(prdata)[grepl("Igh", rownames(prdata))] ## okay
A4300  <- rownames(prdata)[grepl("A4300", rownames(prdata))]
CGenes        <- c( Gm, RP,  Hsp, "Scgb1a1", "Jchain", Igl, Igh, "Igkc", "Malat1", "Xist", A4300, "Mid1", "Kcnq1ot1")

#RaceID parameter
mintotal <- 3000
minexpr <- round(5*mintotal/3000,0)
outminc <- round(5*mintotal/3000,0)
ccor <- 0.65
FGenes <- NULL 
#Run RaceID
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=200)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

getype2 <- function(ma,y, y1){
  z <- sum(grepl(y, ma))
  ma[grep( y, ma)] <- paste(y1, 1:z, sep="_")
  return(ma)
}
naming <- getype2(colnames(sc@ndata), "wt", "uninfected lung Figure1")
naming <- getype2(naming, "Pro", "ProgenitorInput")
marker_col <- brewer.pal(3, "Set1")[1:2]
plotsymbolsmap(sc, types = sub("\\_.+", "", naming), samples_col = marker_col, cex = 1)

##histograms Figure S4D
#load uninfected lung partition 
sc_wt_part <- readRDS("~/Desktop/RaceIDobjects_mainFigures/sc_wtlung.RData")@cpart

progenitor_cluster <- names(sc_wt_part)[sc_wt_part %in% c(1,5,11)]
non_progenitor_cluster <- names(sc_wt_part)[!names(sc_wt_part) %in% progenitor_cluster]
input <- colnames(sc@ndata)[grep("Progenit", colnames(sc@ndata))]
#scaled normalized data 
ndata <- as.matrix(sc@ndata * min(sc@counts)) + 0.1
ndata <- data.frame(ndata)

# Il1rl1 transcript histogram in the different populations 
d1 <- hist(as.numeric(ndata["Il1rl1", colnames(ndata) %in% non_progenitor_cluster]), breaks = 30)
d2 <- hist(as.numeric(ndata["Il1rl1",colnames(ndata) %in% progenitor_cluster]), breaks = d1$breaks, add=TRUE)
d3 <- hist(as.numeric(ndata["Il1rl1",input]), breaks = d1$breaks, add=T)
d1$counts <- d1$counts/sum(d1$counts)
d2$counts <- d2$counts/sum(d2$counts)
d3$counts <- d3$counts/sum(d3$counts)
plot( d1, col=rgb(0,0,1,1/4), xlim = range(0:25), ylim = range(0:1), xlab = "No. Il1rl1 transcripts", main = "Histogram Il1rl1 transcripts mature ILC2")
plot( d2, col=rgb(1,0,0,1/4), xlim = range(0:25), ylim = range(0:1), xlab = "No. Il1rl1 transcripts", main = "Histogram Il1rl1 transcripts immature ILCs")
plot( d3, col=rgb(0,1,0,1/4), xlim = range(0:25), ylim = range(0:1), xlab = "No. Il1rl1 transcripts", main = "Histogram Il1rl1 transcripts Input")

# Il17rb transcript histograms in the different populations 
d1 <- hist(as.numeric(ndata["Il17rb", colnames(ndata) %in% non_progenitor_cluster]), breaks = 30)
d2 <- hist(as.numeric(ndata["Il17rb",colnames(ndata) %in% progenitor_cluster]), breaks = d1$breaks, add=TRUE)
d3 <- hist(as.numeric(ndata["Il17rb",input]), breaks = d1$breaks, add=T)
d1$counts <- d1$counts/sum(d1$counts)
d2$counts <- d2$counts/sum(d2$counts)
d3$counts <- d3$counts/sum(d3$counts)
plot( d1, col=rgb(0,0,1,1/4), xlim=range(0:6), ylim = range(0:1), xlab = "No. Il17rb transcripts", main = "Histogram Il17rb transcripts mature ILC2s")
plot( d2, col=rgb(1,0,0,1/4), xlim=range(0:6), ylim = range(0:1), xlab = "No. Il17rb transcripts", main = "Histogram Il17rb transcripts immature ILCs")
plot( d3, col=rgb(0,1,0,1/4), xlim=range(0:6), ylim = range(0:1), xlab = "No. Il17rb transcripts", main = "Histogram Il17rb transcripts Il18r1+Icos+ST2-Il17rb- Input")

### Code for BM culture S4F 

### data laoding and pre-processing and renaming
prdata <- readRDS("~/Desktop/RawData_FigureReproduce/BMculture1.RData")
prdata <- prdata[, grep("IL18pos|BMDL1", colnames(prdata), invert = T)]
prdata2 <- readRDS("~/Desktop/RawData_FigureReproduce/culture2.RData")
prdata2 <- prdata2[, grep("BMOP9Il18negP7P5|BMOP9Il18negP7P5P6|BMOP9Il18negP8\\_\\d\\_", colnames(prdata2))]

prdata <- cbind(prdata, prdata2)
prdata <- getype2(prdata, "BMOP9Il18negP8", "ST2Il17rbDN")
prdata <- getype2(prdata, "BMOP9Il18negP7P5|BMOP9Il18negP7P5P6", "ST2neg")
prdata <- getype2(prdata, "BMOP9IL18neg", "Pan")

prdata <- prdata[grep("^(ERCC|mt)",row.names(prdata),invert=TRUE),]
cs <- apply(prdata,2,sum)
prdata <- prdata[,cs > 500]
cs <- cs[cs > 500]
f <- t(prdata["Kcnq1ot1",])/cs < .02
prdata <- prdata[,f]

### consistent set of genes for removed from clustering stored in CGenes 
Gm            <- rownames(prdata)[grepl("Gm\\d+", rownames(prdata))] 
RP            <- rownames(prdata)[grepl("RP", rownames(prdata))] 
Hsp           <- rownames(prdata)[grepl("Hsp", rownames(prdata))] 
Igl            <- rownames(prdata)[grepl("Igl", rownames(prdata))] 
Igh           <- rownames(prdata)[grepl("Igh", rownames(prdata))] 
A4300  <- rownames(prdata)[grepl("A4300", rownames(prdata))]
CGenes <- c( Gm, RP,  Hsp, "Scgb1a1", "Jchain", Igl, Igh, "Igkc", "Malat1", "Xist", A4300, "Mid1", "Kcnq1ot1", "Gzmb", "Ucp2", "Calca")

### RaceID paramenter
mintotal <- 3000
minexpr <- round(5*mintotal/3000,0)
outminc <- round(5*mintotal/3000,0)
ccor <- 0.65
FGenes <- NULL 
### Run RaceID
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=30)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

###fracdotplot
genes <- c("Gata3", "Rorc", "Tbx21", "Il1rl1", "Il17rb", "Il18r1", "Zbtb16", "Tcf7", "Mki67", "Id2", "Il13", "Il5", "Areg", "Cxcl2", "Calca", "Cxcr3", "Gzmc", "Cxcr5", "Fcer1g", "Klrb1b", "Bace2")
fracdotplot(sc, genes, samples = unique(sub("\\_.+", "", colnames(sc@ndata))), limup = 4, limdown = -4, zsc = F)
