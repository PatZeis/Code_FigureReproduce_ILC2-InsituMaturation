### Code_Figure4_Figure5_S5_S6
###R version 3.5.1 ### RaceID3 version v0.1.3 

### customm functions to load
source("~/Desktop/Functions_for_code_upload/Functions_to_load.R")

### loading and preprocessing 
prdata <- readRDS("~/Desktop/RawData_FigureReproduce/lung_full_blood.RData")
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
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=200)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

###removal of contaminant clusters
out <- unique(c(names(sc@cpart)[sc@cpart %in% c(10, 11, 17, 19, 23)], names(sc@cluster$kpart)[sc@cluster$kpart %in% c(3,5)]))
prdata <- prdata[,!colnames(prdata) %in% out]


### rerun RaceID
sc <- SCseq(prdata)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=200)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

### StemID lineage inference 
ltree     <- list()
ltr <- Ltree(sc)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr=20,nmode=T,fr=FALSE) 
ltr <- projback(ltr,pdishuf=2000,fast=F,rseed=17000)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.01)

### Derive pseudo-temporal order along the Cxcl2 trajectory(12>5>2>8>4>23), Il13 trajectory(12>5>2>8>6>16) and ivCD45 trajectory(15>2>8>6>16)
fcol <- sc@fcol
ndata <- data.frame(as.matrix(sc@ndata * min(sc@counts)) + 0.1)

samp <- c("wt", "d4", "d7", "d10", "iv45", "Blood", "d15")

n12 <- cellsfromtree(ltr, c(12,5, 2, 8, 6, 16))
ndata12 <- ndata[sc@genes,n12$f]
y12    <- sc@cpart[n12$f]
y12samp <- y12
for ( i in 1:length(samp)) { y12samp[grep(samp[i], names(y12samp))] <- i}


  
n14 <- cellsfromtree(ltr, c(12,5,2 ,8, 4, 23))
ndata14 <- ndata[sc@genes,n14$f]
y14  <- sc@cpart[n14$f]
y14samp <- y14
for ( i in 1:length(samp)) { y14samp[grep(samp[i], names(y14samp))] <- i}

n2 <- cellsfromtree(ltr, c(15, 2, 8, 6, 16))
ndata2 <- ndata[sc@genes,n2$f]
y2    <- sc@cpart[n2$f]
y2samp <- y2
for ( i in 1:length(samp)) { y2samp[grep(samp[i], names(y2samp))] <- i}
  
### generate plots
#tSNE map of cluster
plotmap(sc)
#tSNE map of candidate genes 
plotexpmap(sc, "Il18r1", logsc = T)
#fractiondotplot cluster
genes2 <- c( "Tbx21","Rorc", "Gata3", "Il1rl1", "Il17rb","Il2ra", "Il18r1", "Icos", "Thy1", "Il7r", "Zbtb16", "Tcf7", "Cd7", "Igf1r","Itgae", "Cd3e",  "S100a6","S1pr1","Vim","Klf2","Cd69","Mki67", "Id2","Kit","Arg1","Klrg1","Cxcl2", "Csf2", "Il2", "Calca","Il5","Il13","Pdcd1", "Ctla4","Il1r2","H2-Ab1", "Areg", "Klf4", "Rgs2")
cluster <- c(12, 1, 3, 2, 15, 5,  8, 4, 10, 23, 14, 6, 16, 17, 7, 19, 29, 9)
fracdotplot(sc, genes = genes2,cluster = cluster, limup = 1,limdown =  -1) ## z-core of the mean
fracdotplot(sc, genes = genes2,cluster = cluster, limup = 4,limdown =  -4, zsc=F) ## log2 mean

### sample renaming
getype2 <- function(ma,y, y1){
  z <- sum(grepl(y, colnames(ma)))
  colnames(ma)[grep( y, colnames(ma))] <- paste(y1, 1:z, sep="_")
  return(ma)
}

sc@ndata <- getype2(sc@ndata, "d15", "g.d15")
sc@ndata <- getype2(sc@ndata, "Blood", "f.d15Blood")
sc@ndata <- getype2(sc@ndata, "d4", "b.d4")
sc@ndata <- getype2(sc@ndata, "d7", "c.d7")
sc@ndata <- getype2(sc@ndata, "d10","d.d10")
sc@ndata <- getype2(sc@ndata, "iv45", "e.d10ivCD45")
sc@ndata <- getype2(sc@ndata, "wt","a.uninfected")

### continue generating plots 

#fractionplot samples
types <- unique(sub("\\_.+", "", colnames(sc@ndata)))
types <- types[order(types)]
fracdotplot(sc, genes = genes2, samples = types, limup = 1, limdown = -1)

#tSNE map of sample identity
set2 <- brewer.pal(7, "Set2")
set1 <- brewer.pal(7, "Set1")
set1 <- c(set1[7], set2[6],set1[5], set1[3], set1[4],  set1[1], set1[2]) 
plotsymbolsmap(sc, types = sub("\\_.+", "", colnames(sc@ndata)), samples_col = set1, cex = 1, leg = F, map = T)

# plot uninfected_lung_cluster in time-course data
cluster_partition_uninfected <- readRDS("~/Desktop/RaceIDobjects_mainFigures/sc_wtlung.RData")
cluster_partition_uninfected <- cluster_partition_uninfected@cpart
cl1.uninfected <- names(cluster_partition_uninfected)[cluster_partition_uninfected == 1]
getype2 <- function(ma,y, y1){
  z <-sum(  ma %in% y)
  ma[as.numeric(which(ma %in% y))] <- paste(y1, 1:z, sep="_")
  return(ma) }
symbol <- getype2(names(sc@cpart), cl1.uninfected, "cl1uninfected")
plotsymbolsmap(sc, sub("\\_.+", "", symbol), subset = "cl1uninfected")

#pie charts normalized sample contribution to clusters
types <- sub("\\_.+", "", colnames(sc@ndata))
types_num <- as.numeric(table(types))
types <- unique(types)
types <- types[order(types)]

coloc <- set1
cpart <- sc@cpart
names(cpart) <- colnames(sc@ndata)
a <- as.numeric(names(table(cpart)))
for ( i in a){
  x <- table(as.character(sapply(names(cpart[cpart==i]), function(x) str_split(x, pattern="[_]")[[1]][1])))
  lbls <- names(x)
  cat(paste(lbls, "\n", sep=""))
  y <- which(types %in% lbls)  
  slices <- as.numeric(x)/types_num[y] 
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(pct, " % ", lbls,sep=" ")
  pie(slices, labels=lbls,cex.main=1.5, cex = 1.25,col=coloc[y], cex.main=2,main=paste("Cluster ", i, " distribution ", " N=", length(cpart[cpart == i]), sep=" "))
}
#pie charts uninfected cluster distribution in time-course cluster
uninf_clust <- c(1,5,11,9,8)
clus_list <- list()
for ( i in 1:length(uninf_clust)) { 
  clus_list[[i]] <- names(cluster_partition_uninfected)[cluster_partition_uninfected == uninf_clust[i]]
  names(clus_list)[i] <- paste("clwt", uninf_clust[i], sep="" )}

for ( i in 1:length(clus_list)) { 
  x <- table(sc@cpart[names(sc@cpart) %in% clus_list[[i]]])
  x <- x[order(x, decreasing = T)]
  lbls <- names(x)
  slices <- as.numeric(x)
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(pct, "%", "cl", lbls,sep=" ")
  pie(slices, labels=lbls, cex.main=1.5, cex=1, col=sc@fcol[as.numeric(names(x))], main=paste("Sample ", names(clus_list)[i], " distribution ", " N=", sum(clus_list[[i]] %in% names(sc@cpart)), sep=" "))
}
# pie charts time-course sample distribution over cluster
sample <- unique(sub("\\_.+", "", colnames(sc@ndata)))
sample <- sample[order(sample)]
cluster <- as.numeric(names(table(sc@cpart)[table(sc@cpart) >= 20]))
cpart <- cpart[cpart %in% cluster]

for ( i in 1:length(sample)) { 
  x <- table(cpart[grep(sample[i], names(cpart))])
  x <- x[order(x, decreasing = T)]
  lbls <- names(x)
  slices <- as.numeric(x)
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(pct, "%", "cl", lbls,sep=" ")
  pie(slices, labels=lbls, cex.main=1.5, cex=1, col=sc@fcol[as.numeric(names(x))], main=paste("Sample ", sample[i], " distribution ", " N=", length(grep(sample[i], colnames(sc@ndata))), sep=" "))

}

#StemID lineage tree
plotgraph(ltr,showCells=FALSE,showTsne=TRUE,tp=.3,scthr=0)

##derive pseudo-temporal gene expression 
il13_up <- c("Gata3", "Il1rl1", "Arg1", "Klrg1", "Bcl11b", "Il13", "Neurl3", "Furin")
cxcl2_up <- c("Gata3", "Il1rl1", "Arg1", "Klrg1", "Bcl11b","Csf2", "Cxcl2")
downgenes <- c("Il18r1", "Tcf7", "Zbtb16","Cd7", "Ikzf2", "Maf")
migratory_up <- c("Gata3", "Il1rl1", "Bcl11b", "Furin", "Neurl3", "Il13")
migratory_down <- c("Vim", "Klrg1", "Cd48", "S1pr1", "Itgb7", "Itgb1", "Itga4", "Klf2")

#trajectory il13 (12>5>2>8>6>16)
plotexpression(ndata12, y12, il13_up, n12$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), samp = T, samppart = y12samp, sampcol = set1)
plotexpression(ndata12, y12, downgenes, n12$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), samp = T, samppart = y12samp, sampcol = set1, legendpos = "topright")

#trajectory cxcl2 (12>5>2>8>4>23
plotexpression(ndata14, y14, cxcl2_up, n14$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), samp = T, samppart = y14samp, sampcol = set1)
plotexpression(ndata14, y14, downgenes, n14$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), samp = T, samppart = y14samp, sampcol = set1)

#trajectory migratory (15>2>8>6>16)
plotexpression(ndata2, y2, migratory_up, n2$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), samp = T, samppart = y2samp, sampcol = set1)
plotexpression(ndata2, y2, migratory_down, n2$f, col = sc@fcol, cluster = F, alpha=.5, lwd = 2.5, leg = T, ylim = range(-0.5:10), samp = T, samppart = y2samp, sampcol = set1, legendpos = "topright")

# derive projection of cells of cluster12 
proj <- getproj(ltr,i=12)

#plot3marker genes for comparing expression of 3 genes between 2 clusters
colnames(sc@ndata) <- names(sc@cpart)
plot3marker(sc, cluster = c(6,12), gene1 = "Il18r1", gene2 = "Tcf7", gene3 = "Il1rl1")

### Blood
### filter for Blood cells
prdata3 <- prdata[, grep("Blood", colnames(prdata))]

### Run RaceID
sc <- SCseq(prdata3)
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=1, LBatch=NULL, knn=10, FGenes=FGenes, CGenes=CGenes, ccor=ccor, bmode="RaceID")
sc <- compdist(sc,metric="pearson",FSelect=T,knn=10)
sc <- clustexp(sc,sat=TRUE,samp=1000,cln=NULL,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
sc <- comptsne(sc,rseed=15555,perplexity=30)
sc <- findoutliers(sc,probthr=0.001,outminc=outminc,outlg=2,outdistquant=0.95)

### generating plots
#tSNE map cluster identity
plotmap(sc)
#fracdotplots
genes2 <- c( "Tbx21","Rorc", "Gata3", "Il1rl1", "Il17rb","Il2ra", "Il18r1", "Icos", "Thy1", "Il7r", "Zbtb16", "Tcf7", "Cd7", "Igf1r","Itgae", "Cd3e",  "S100a6","S1pr1","Vim","Klf2","Cd69","Mki67", "Id2","Kit","Arg1","Klrg1","Cxcl2", "Csf2", "Il2", "Calca","Il5","Il13","Pdcd1", "Ctla4","Il1r2","H2-Ab1", "Areg", "Klf4", "Rgs2")
cluster <- c(1:13)
fracdotplot(sc, genes2, cluster, limup = 1, limdown = -1)
