data   <- list()
gene2iso <- read.csv("./wgEncodeGencodeBasicVM9_clean_genes2groups_2112016.tsv",sep="\t",header=FALSE) ### add path to file
ercc     <- read.csv("./ERCC_Controls_Analysis_length_ext.txt",sep="\t",header=TRUE) ### add path to file
for ( n in grep("nb.cell",names(ercc)) ){ercc[,n] <- 5/6 * ercc[,n]} # original CEL-seq
n4 <- 1:193
z <- list()
nl <- c(Progenit_1="", Progenit_2="", Progenit_3="", Progenit_4="", LungOP9Il18negP8_1="", LungOP9Il18negP8_2="", LungOP9Il18negP8_3="", BMOP9Il18negP7P5="", BMOP9Il18negP7P5P6="", BMOP9Il18negP8_1="", BMOP9Il18negP8_2="");
for ( i in names(nl) ) nl[i] <- i
fl <- list()
for ( i in names(nl) ) fl[[i]] <- n4
for ( sl in names(nl) ){
    cat("loading ",sl,"\n")
    x <- read.csv(paste("./",sl,".coutt",".csv",sep=""),sep="\t",header=TRUE) ### path to processed data files (.coutt) files
    x <- x[,fl[[sl]]]
    x <- merge(data.frame(GENEID=c(as.vector(gene2iso[,1]),as.vector(ercc[,1])),GROUP=c(as.vector(gene2iso[,2]),as.vector(ercc[,1]))),x,by="GENEID",all.x=TRUE)[,-1]
    names(x)[1] <- "GENEID"
    x[is.na(x[,2]),-1] <- 0
    x <- x[order(x$GENEID),]
    z[[sl]]  <- x
    names(z[[sl]])  <- c("GENEID",paste(nl[[sl]],sub("X","",names(z[[sl]])[-1]),sep="_"))
    rownames(z[[sl]]) <- x$GENEID 
}
for ( i in 1:length(z) ){ 
    cat("merging ",names(z)[[i]],"\n")
    y <- if ( i == 1 ) z[[i]] else cbind(y,z[[i]][rownames(y),-1])
}
row.names(y) <- y$GENEID
y <- y[,-1]
y <- y[,apply(y,2,sum)>500] 

data <- y

         
bm_in_lung_number <- c(1, 17, 33, 49, 65, 81, 25, 41, 57, 73, 89)
bm_in_lung <- paste("LungOP9Il18negP8", 1, bm_in_lung_number, sep = "_")
bmP7P6 <- paste("BMOP9Il18negP7P5P6", 145:192, sep="_")
getype2 <- function(ma,y, y1){
    z <-sum(  colnames(ma) %in% y)
    colnames(ma)[as.numeric(which(colnames(ma) %in% y))] <- paste(y1, 1:z, sep="_")
    return(ma) }
data <- getype2(data, bm_in_lung, "BMOP9Il18negP8")
data <- getype2(data, bmP7P6, "BMOP9Il18negST2negLinPosOutput")
saveRDS(data, "culture2.RData")


