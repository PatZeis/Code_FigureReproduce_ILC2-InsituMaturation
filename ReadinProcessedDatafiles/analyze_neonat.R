data   <- list()
gene2iso <- read.csv("./wgEncodeGencodeBasicVM9_clean_genes2groups_2112016.tsv",sep="\t",header=FALSE) ### add path to file
ercc     <- read.csv("./ERCC_Controls_Analysis_length_ext.txt",sep="\t",header=TRUE) ### add path to file
for ( n in grep("nb.cell",names(ercc)) ){ercc[,n] <- 5/6 * ercc[,n]} # original CEL-seq
n4 <- 1:193
z <- list()
nl <- c(NeoP8_1="", NeoP8_2="", NeoP8_3="", NeoP8_4="");
for ( i in names(nl) ) nl[i] <- i
fl <- list()
for ( i in names(nl) ) fl[[i]] <- n4
for ( sl in names(nl) ){
    cat("loading ",sl,"\n")
    x <- read.csv(paste("./",sl,".coutt",".csv",sep=""),sep="\t",header=TRUE) ### path to processed data files (.coutt.csv) files
    x <- x[,fl[[sl]]]
    x <- merge(data.frame(GENEID=c(as.vector(gene2iso[,1]),as.vector(ercc[,1])),GROUP=c(as.vector(gene2iso[,2]),as.vector(ercc[,1]))),x,by="GENEID",all.x=TRUE)[,-1]
    names(x)[1] <- "GENEID"
    x[is.na(x[,2]),-1] <- 0
    x <- x[order(x$GENEID),]
    z[[sl]]  <- x
    names(z[[sl]])  <- c("GENEID",paste(nl[[sl]],sub("X","",names(z[[sl]])[-1]),sep="_"))
    rownames(z[[sl]]) <- x$GENEID # new
}
for ( i in 1:length(z) ){ # new
    cat("merging ",names(z)[[i]],"\n")
    y <- if ( i == 1 ) z[[i]] else cbind(y,z[[i]][rownames(y),-1])
}

row.names(y) <- y$GENEID
y <- y[,-1]
data <- y

vec <- 0:23 
CD45_2 <- c(1,9,2,10,3,11, 4, 12) 
CD45_1 <- c(5,13,6,14,7,15,8,16)
 
givenumbers <- function(rows){
    vec <- 0:23
    numbers <- c()
    for (i in 1:length(rows)) {
        for (n in 1:length(vec)){
            numbers <- append(numbers, rows[i]+16*vec[n])
        }
    }
    return(numbers)
}
CD45_2_number <- givenumbers(CD45_2) 
CD45_1_number <- givenumbers(CD45_1)
      
givecells <- function(numbers, descri, a=c(1,2,3,4)) {
    cells <- c()
    for (i in 1:length(numbers)){
        if(numbers[i] %in% 1:96){
            cells <- append(cells, paste(descri, a[1], numbers[i], sep="_"))
        }
        if(numbers[i] %in% 97:192){
            cells <- append(cells, paste(descri, a[2],numbers[i], sep="_"))
        }
        if(numbers[i] %in% 193:288){
            cells <- append(cells, paste(descri, a[3], numbers[i]-192, sep="_"))
        }
        if(numbers[i] %in% 289:384){
            cells <- append(cells, paste(descri, a[4], numbers[i]-192, sep="_"))
        }
    }
    cells
} 

NeonatST2 <- givecells(CD45_2_number, "NeoP8")
NeonatST2neg <- givecells(CD45_1_number, "NeoP8")

n <- names(data)
names(n) <- n
n[NeonatST2] <- sub("NeoP8","NeoST2",NeonatST2)
n[NeonatST2neg] <- sub("NeoP8","NeoST2neg", NeonatST2neg)

names(data) <- as.character(n)
saveRDS(data, "neonat.RData")
  
