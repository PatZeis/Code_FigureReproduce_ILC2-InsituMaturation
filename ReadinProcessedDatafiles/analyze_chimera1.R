
data   <- list()
gene2iso <- read.csv("./wgEncodeGencodeBasicVM9_clean_genes2groups_2112016.tsv",sep="\t",header=FALSE) ### add path to file
ercc     <- read.csv("./ERCC_Controls_Analysis_length_ext.txt",sep="\t",header=TRUE) ### add path to file
for ( n in grep("nb.cell",names(ercc)) ){ercc[,n] <- 5/6 * ercc[,n]} # original CEL-seq
n4 <- 1:193
 z <- list()
nl <- c(d15chim1_1="",d15chim1_2="", d15chim1_3="", d15chim1_4="", d15chim2_1="", d15chim2_2="", d15chim2_3="", d15chim2_4="");
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

      
    
  
vec <- 0:23 
### top half of the plate = host
### bottom half of the plate = donor

host <- c(1,9,2,10,3,11, 4, 12)
donor <- c(5,13,6,14,7,15,8,16)
donor_2 <- c(5,13,6,14)
donor_3 <- c(7,15,8,16) 
 
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
chim1_host_number <- givenumbers(host)
chim1_donor_number <- givenumbers(donor)
chim2_host_number <- givenumbers(host)
chim2_donor_number <- givenumbers(donor_2)
chim3_host_number <- givenumbers(donor_3)
      
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

chim1_host <- givecells(chim1_host_number, "d15chim1")
chim1_donor <- givecells(chim1_donor_number, "d15chim1")
chim2_host <- givecells(chim2_host_number, "d15chim2")
chim2_donor <- givecells(chim2_donor_number, "d15chim2")
chim3_donor <- givecells(chim3_host_number, "d15chim2")
      
n <- names(data)
	
names(n) <- n
n[chim1_host] <- sub("d15chim1","d15C12",chim1_host)
n[chim1_donor] <- sub("d15chim1","d15C11",chim1_donor)
n[chim2_host] <- sub("d15chim2","d15C22",chim2_host)
n[chim2_donor] <- sub("d15chim2","d15C21",chim2_donor)
n[chim3_donor] <- sub("d15chim2","d15C31",chim3_donor)
names(data) <- as.character(n)
      
getype2 <- function(ma,y, y1){
    z <- sum(grepl(y, colnames(ma)))
    colnames(ma)[grep( y, colnames(ma))] <- paste(y1, 1:z, sep="_")
    return(ma)
}
data <- getype2(data, "d15C12", "std15CIhost")
data <- getype2(data, "d15C11", "std15CIdon")
data <- getype2(data, "d15C22", "std15CIIhost")
data <- getype2(data, "d15C21", "std15CIIdon")
data <- getype2(data, "d15C31", "std15CIIIdon")
data <- data[grep("d15chim", colnames(data), invert=T)]
saveRDS(data, "chimera1.RData")

  


