### Code_Figure2G-H
###R version 3.5.1 ### RaceID3 version v0.1.3 

### custom functions to load
source("~/Desktop/Functions_for_code_upload/Functions_to_load.R")

### data loading and preprocessing
prdata <- readRDS("~/Desktop/RawData_FigureReproduce/Il18r1_Icos_BM.RData")
prdata <- prdata[grep("^(ERCC|mt)",row.names(prdata),invert=TRUE),]
cs <- apply(prdata,2,sum)
prdata <- prdata[,cs > 500]
cs <- cs[cs > 500]
f <- t(prdata["Kcnq1ot1",])/cs < .02
prdata <- prdata[,f]

### generic set of genes removed from clustering stored in CGenes 
Gm            <- rownames(prdata)[grepl("Gm\\d+", rownames(prdata))] ### also FGenes also functional genes
RP            <- rownames(prdata)[grepl("RP", rownames(prdata))] ### also FGenes same as Rik
Hsp           <- rownames(prdata)[grepl("Hsp", rownames(prdata))] ## okay
Igl            <- rownames(prdata)[grepl("Igl", rownames(prdata))] ## okay
Igh           <- rownames(prdata)[grepl("Igh", rownames(prdata))] ## okay
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
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=6,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=30)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

### generation of plots
#tSNE map of cluster
plotmap(sc)
#fractiondotplot
genes2 <- c( "Tbx21","Rorc", "Gata3", "Il1rl1", "Il17rb","Il2ra", "Il18r1", "Icos", "Thy1", "Il7r", "Zbtb16", "Tcf7", "Cd7", "Igf1r","Itgae", "Cd3e",  "S100a6","S1pr1","Vim","Klf2","Cd69","Mki67", "Id2","Kit","Arg1","Klrg1","Cxcl2", "Csf2", "Il2", "Calca","Il5","Il13","Pdcd1", "Ctla4","Il1r2","H2-Ab1", "Areg", "Klf4", "Rgs2")
cluster <- c(1,4,5,6,3,7,2)
fracdotplot(sc, genes2, cluster = cluster, limup = 1, limdown = -1)
