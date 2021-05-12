data   <- list()
gene2iso <- read.csv("./wgEncodeGencodeBasicVM9_clean_genes2groups_2112016.tsv",sep="\t",header=FALSE) ### add path to file
ercc     <- read.csv("./ERCC_Controls_Analysis_length_ext.txt",sep="\t",header=TRUE) ### add path to file
for ( n in grep("nb.cell",names(ercc)) ){ercc[,n] <- 5/6 * ercc[,n]} # original CEL-seq
n4 <- 1:193
z <- list()
nl <- c(BMILCPP8_1="", BMILCPP8_2="", BMILCPP8_3="", BMILCPP8_4="", BMILCPP9_1="",BMILCPP9_2="", BMILCPP9_3="", BMILCPP9_4="");
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
    rownames(z[[sl]]) <- x$GENEID 
}
for ( i in 1:length(z) ){ 
    cat("merging ",names(z)[[i]],"\n")
    y <- if ( i == 1 ) z[[i]] else cbind(y,z[[i]][rownames(y),-1])
}
row.names(y) <- y$GENEID
y <- y[,-1]
data <- y

double_P8 <- c(paste("BMILCPP8_1", 1:96, sep="_"), paste("BMILCPP8_2", 97:(97+31), sep="_"))
il18r1_P8 <- c(paste("BMILCPP8_2", 129:192, sep="_"), paste("BMILCPP8_3", 1:64, sep="_"))
icos_P8 <- c(paste("BMILCPP8_3", 65:96, sep="_"), paste("BMILCPP8_4", 97:192, sep="_"))
double_P9 <- c(paste("BMILCPP9_1", 1:96, sep="_"), paste("BMILCPP9_2", 97:(97+31), sep="_"))
il18r1_P9 <- c(paste("BMILCPP9_2", 129:192, sep="_"), paste("BMILCPP9_3", 1:64, sep="_"))
icos_P9 <- c(paste("BMILCPP9_3", 65:96, sep="_"), paste("BMILCPP9_4", 97:192, sep="_"))

getype2 <- function(ma,y, y1){
    z <-sum(  colnames(ma) %in% y)
    colnames(ma)[as.numeric(which(colnames(ma) %in% y))] <- paste(y1, 1:z, sep="_")
    return(ma)
}

data <- getype2(data, double_P8, "Il18r1IcosMouse1")
data <- getype2(data, il18r1_P8, "Il18r1Mouse1")
data <- getype2(data, icos_P8, "IcosMouse1")
data <- getype2(data, double_P9, "Il18r1IcosMouse2")
data <- getype2(data, il18r1_P9, "Il18r1Mouse2")
data <- getype2(data, icos_P9, "IcosMouse2")
data <- data[,grep("^Icos|Il18r1M", colnames(data), invert=T)]
data <- data[,grep("BMILCP", colnames(data), invert=T)]
      
saveRDS(data, "Il18r1_Icos_BM.RData") 
  


