
data   <- list()
gene2iso <- read.csv("./wgEncodeGencodeBasicVM9_clean_genes2groups_2112016.tsv",sep="\t",header=FALSE) ### add path to file
ercc     <- read.csv("./ERCC_Controls_Analysis_length_ext.txt",sep="\t",header=TRUE) ### add path to file
for ( n in grep("nb.cell",names(ercc)) ){ercc[,n] <- 5/6 * ercc[,n]} # original CEL-seq
n4 <- 1:193
  
z <- list()
nl <- c(DL1IL8neg_1_13="", DL1IL18neg_2_14="", DL1IL18pos_1_15="", DL1IL18pos_2_16="", OP9IL18neg_1_9="", OP9IL18neg_2_10="", OP9IL18pos_1_11="", OP9IL18pos_2_12=""  );
for ( i in names(nl) ) nl[i] <- i
fl <- list()
for ( i in names(nl) ) fl[[i]] <- n4
files <- list.files(path = ".", pattern = ".coutt.csv$") ### path to decompressed processed data files (.coutt.csv) files
for ( sl in names(nl) ){
    cat("loading ",sl,"\n")
    x <- paste("_",sl,".coutt",".csv",sep="") 
    x <- read.csv(grep(x, files, value = T),sep="\t",header=TRUE)
    x <- x[,fl[[sl]]]
    x <- merge(data.frame(GENEID=c(as.vector(gene2iso[,1]),as.vector(ercc[,1])),GROUP=c(as.vector(gene2iso[,2]),as.vector(ercc[,1]))),x,by="GENEID",all.x=TRUE)[,-1]
    names(x)[1] <- "GENEID"
    x[is.na(x[,2]),-1] <- 0
    x <- x[order(x$GENEID),]
    z[[sl]]  <- x
    names(z[[sl]])  <- c("GENEID",paste(nl[[sl]],sub("X","",names(z[[sl]])[-1]),sep="_"))
    rownames(z[[sl]]) <- x$GENEID # new
}
for ( i in 1:length(z) ){ 
    cat("merging ",names(z)[[i]],"\n")
    y <- if ( i == 1 ) z[[i]] else cbind(y,z[[i]][rownames(y),-1])
}
row.names(y) <- y$GENEID
y <- y[,-1]
y <- y[,apply(y,2,sum)>500]
data <- y
colnames(data) <- paste("BM", colnames(data), sep = "")
saveRDS(data, "BMculture1.RData") 
  


