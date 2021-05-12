data   <- list()
gene2iso <- read.csv("./wgEncodeGencodeBasicVM9_clean_genes2groups_2112016.tsv",sep="\t",header=FALSE)  ### add path to file
ercc     <- read.csv("./ERCC_Controls_Analysis_length_ext.txt",sep="\t",header=TRUE)  ### add path to file
for ( n in grep("nb.cell",names(ercc)) ){ercc[,n] <- 5/6 * ercc[,n]} # original CEL-seq
n4 <- 1:193
z <- list()
nl <- c(Blood_1="Blood_1", Blood_2="Blood_2", Blood_3="Blood_3", Blood_4="Blood_4",mergedInf_M3_P1_I5="d15m3_1",mergedInf_M3_P1_I6="d15m3_2", mergedInf_M3_P1_I7="d15m3_3", mergedInf_M3_P1_I8="d15m3_4", Inf_M3_P2_I5="d15m3_5", Inf_M3_P2_I8="d15m3_6", Inf_M3_P2_I9="d15m3_7", mergedInf_M4_P1_I9="d15m4_1", mergedInf_M4_P1_I10="d15m4_2", mergedInf_M4_P1_I11="d15m4_3", mergedInf_M4_P1_I12="d15m4_4", Inf_M4_P2_I1="d15m4_5", Inf_M4_P2_I2="d15m4_6", Inf_M4_P2_I3="d15m4_7", Inf_M4_P2_I4="d15m4_8", mergedWT_M1_P1_I1="wt1_1", mergedWT_M1_P1_I2="wt1_2", mergedWT_M1_P1_I3="wt1_3", mergedWT_M1_P1_I4="wt1_4", WT_M1_P2_I9="wt1_5", WT_M1_P2_I10="wt1_6", WT_M1_P2_I11="wt1_7", WT_M1_P2_I12="wt1_8", WT_M2_P1_I5="wt2_1", WT_M2_P1_I6="wt2_2", WT_M2_P1_I7="wt2_3", WT_M2_P1_I8="wt2_4", WT_M2_P2_I1="wt2_5", WT_M2_P2_I2="wt2_6", WT_M2_P2_I3="wt2_7", WT_M2_P2_I4="wt2_8", d7_1="d7_1", d7_2="d7_2", d7_3="d7_3", d7_4="d7_4", d10_1="d10_1", d10_2="d10_2", d10_3="d10_3", d10_4="d10_4", iv45_1="iv45_1", iv45_2="iv45_2", iv45_3="iv45_3", iv45_4="iv45_4", infLd7_1="d7m1_1", infLd7_2="d7m1_2", infLd7_3="d7m1_3", infLd7_4="d7m1_4", infLd7_5="d7m1_5", infLd7_6="d7m1_6", infLd7_7="d7m1_7", d4_1="d4_1", d4_2="d4_2", d4_3="d4_3", d4_4="d4_4",d4_5="d4_5", d4_6="d4_6", d4_7="d4_7", d4_8="d4_8");
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

d41I <- givecells(CD45_2_number, "d4")
d41II <- givecells(CD45_1_number, "d4")
d42I <- givecells(CD45_2_number, "d4", a=c(5,6,7,8))
d42II <- givecells(CD45_1_number, "d4", a=c(5,6,7,8))
    
n <- names(data)
names(n) <- n
n[d41I] <- sub("d4","d4m1",d41I)
n[d41II] <- sub("d4", "d4m2", d41II)
n[d42I] <- sub("d4","d4m1",d42I)
n[d42II] <- sub("d4", "d4m2", d42II)
names(data) <- as.character(n)
saveRDS(data, "lung_full_blood.RData") 


