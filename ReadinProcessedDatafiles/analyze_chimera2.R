data   <- list()
gene2iso <- read.csv("./wgEncodeGencodeBasicVM9_clean_genes2groups_2112016.tsv",sep="\t",header=FALSE) ### add path to file
ercc     <- read.csv("./ERCC_Controls_Analysis_length_ext.txt",sep="\t",header=TRUE) ### add path to file
for ( n in grep("nb.cell",names(ercc)) ){ercc[,n] <- 5/6 * ercc[,n]} # original CEL-seq
n4 <- 1:193
z <- list()
nl <- c(Chim1_1="", Chim1_2="", Chim1_3="", Chim1_4="", Chim1_5="",Chim1_6="", Chim1_7="", Chim1_8="", Chim2_1="", Chim2_2="", Chim2_3="", Chim2_4="", Chim2_5="", Chim2_6="", Chim2_7="", Chim2_8="", ChieldD_1="", ChieldD_2="", ChieldD_3="", ChieldD_4="");
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
host_number <- givenumbers(host) 
donor_number <- givenumbers(donor)
     
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

Chim1host1 <- givecells(host_number , "Chim1")
Chim1host2 <- givecells(host_number , "Chim1", a=c(5,6,7,8))
Chim2host1 <- givecells(host_number , "Chim2")
Chim2host2 <- givecells(host_number , "Chim2", a=c(5,6,7,8))
Chim1don1 <- givecells(donor_number, "Chim1")
Chim1don2 <- givecells(donor_number, "Chim1", a=c(5,6,7,8))
Chim2don1 <- givecells(donor_number, "Chim2")
Chim2don2 <- givecells(donor_number, "Chim2", a=c(5,6,7,8))   
   
ChimdonS1 <- givecells(host_number, "ChieldD")
ChimdonS2 <- givecells(donor_number, "ChieldD")
     
n <- names(data)
     
names(n) <- n
n[Chim1host1] <- sub("Chim1","CI1host",Chim1host1)
n[Chim1host2] <- sub("Chim1","CI2host", Chim1host2)
n[Chim2host1] <- sub("Chim2","CII1host",Chim2host1)
n[Chim2host2] <- sub("Chim2","CII2host", Chim2host2)
n[Chim1don1] <- sub("Chim1","CI1don",Chim1don1)
n[Chim1don2] <- sub("Chim1","CI2don",Chim1don2)
n[Chim2don1] <- sub("Chim2","CII1don",Chim2don1)
n[Chim2don2] <- sub("Chim2", "CII2don", Chim2don2)
n[ChimdonS1] <- sub("ChieldD", "CI3don", ChimdonS1)
n[ChimdonS2] <- sub("ChieldD", "CII3don", ChimdonS2)
                        
names(data) <- as.character(n)
data <- data[,grep("Chim", colnames(data), invert = T)]
colnames(data) <- paste("d15", colnames(data), sep="")
saveRDS(data,"chimera2.RData") 
 


