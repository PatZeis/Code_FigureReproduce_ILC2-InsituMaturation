###Code_Figure1 and S1A-F.
###R version 3.5.1 ### RaceID3 version v0.1.3 

### custom functions to load
source("~/Desktop/Functions_for_code_upload/Functions_to_load.R")

### data loading and preprocessing
prdata <- readRDS("~/Desktop/RawData_FigureReproduce/lung_full_blood.RData")
prdata <- prdata[,grep("wt", colnames(prdata))]
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
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=200)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

### remove cluster with NK cell signature and re-run RaceID
cl10 <- names(sc@cluster$kpart[sc@cluster$kpart == 10]) #### nk cell cluster
prdata <- prdata[,!colnames(prdata) %in% cl10]
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=200)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

### generation of plots
#tSNE map of cluster
plotmap(sc)
#tSNE map of candidate genes 
plotexpmap(sc, "Gata3", logsc = T)
#fractiondotplot
cluster <- c(1  ,5 ,11 , 9 , 8 , 4 ,10 , 12,3 , 6 ,  7  ,2)
genes2 <- c( "Tbx21","Rorc", "Gata3", "Il1rl1", "Il17rb","Il2ra", "Il18r1", "Icos", "Thy1", "Il7r", "Zbtb16", "Tcf7", "Cd7", "Igf1r","Itgae", "Cd3e",  "S100a6","S1pr1","Vim","Klf2","Cd69","Mki67", "Id2","Kit","Arg1","Klrg1","Cxcl2", "Csf2", "Il2", "Calca","Il5","Il13","Pdcd1", "Ctla4","Il1r2","H2-Ab1", "Areg", "Klf4", "Rgs2")
fracdotplot(sc, genes2, cluster = cluster, limup = 1, limdown = -1)

### StemID lineage inference 
ltree     <- list()
ltr <- Ltree(sc)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr=20,nmode=T,fr=FALSE) ##knn removed
ltr <- projback(ltr,pdishuf=2000,fast=F,rseed=17000)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.01)
#StemID lineage tree
plotgraph(ltr,showCells=FALSE,showTsne=TRUE,tp=.3,scthr=0)

### Derive pseudo-temporal order along the Cxcl2 trajectory(5>11>9>8>4>10), pseudo-temporal gene expression profiles and self organizing maps of modules of co-expressed genes along the trajectory  
n <- cellsfromtree(ltr, c(5,11,9,8,4,10)) ### cxcl2 branch
ndata <- data.frame(as.matrix(sc@ndata * min(sc@counts)) + 0.1)
ndata <- ndata[sc@genes,n$f]
s1d <- getsom(ndata,nb=500,alpha=.5)
ps  <- procsom(s1d,corthr=.75,minsom=10) 
fcol <- sc@fcol
y    <- sc@cpart[n$f]

#Plot pseudo-temporal gene expression profiles
genes_to_plot_up <- c("Gata3", "Il1rl1", "Arg1", "Klrg1", "Bcl11b", "Cxcl2")
genes_to_plot_down <- c("Il18r1", "Tcf7", "Zbtb16", "Cd7", "Ikzf2", "Maf")
plotexpression(ndata, y, genes_to_plot_up, n$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), legendpos = "topleft")
plotexpression(ndata, y, genes_to_plot_down, n$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), legendpos = "topright")

### Supplementary figures 

#plot3marker genes 
plot3marker(sc, cluster = 10, gene1 = "Cd69", gene2 = "Klf2", gene3 = "Cxcl2") 

#SOM and pathway analysis of modules
plotheatmap(ps$all.z, xpart=y, xcol=fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=FALSE)
go_analysis(ps$nodes, rownames(sc@ndata), node = 1)
go_analysis(ps$nodes, rownames(sc@ndata), node = 5)

#derive pseudo-temporal gene expression along the Areg trajectory 
n2 <- cellsfromtree(ltr, c(5,11,9,8,4,12))  ### areg branch
ndata2 <- data.frame(as.matrix(sc@ndata * min(sc@counts)) + 0.1)
ndata2 <- ndata2[sc@genes,n2$f]
fcol <- sc@fcol
y2    <- sc@cpart[n2$f]
genes_to_plot_up2 <- c("Gata3", "Il1rl1", "Arg1", "Klrg1", "Bcl11b", "Areg")
plotexpression(ndata2, y2, genes_to_plot_up2, n2$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), legendpos = "topleft")
plotexpression(ndata2, y2, genes_to_plot_down, n2$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), legendpos = "topright")
