data   <- list()
gene2iso <- read.csv("./wgEncodeGencodeBasicVM9_clean_genes2groups_2112016.tsv",sep="\t",header=FALSE)### add path to file
ercc     <- read.csv("./ERCC_Controls_Analysis_length_ext.txt",sep="\t",header=TRUE)### add path to file
for ( n in grep("nb.cell",names(ercc)) ){ercc[,n] <- 5/6 * ercc[,n]} # original CEL-seq
n4 <- 1:193
z <- list()
nl <- c(exp2_Para11_1="Para11_1",exp2_Para11_2="Para11_2", exp2_Para11_3="Para11_3", exp2_Para11_4="Para11_4", exp2_Para11_5="Para11_5", exp2_Para11_6="Para11_6",exp2_Para11_7="Para11_7", exp2_Para11_8="Para11_8", exp2_Para12_1="Para12_1", exp2_Para12_2="Para12_2", exp2_Para12_3="Para12_3", exp2_Para12_4="Para12_4", exp2_Para12_5="Para12_5", exp2_Para12_6="Para12_6", exp2_Para12_7="Para12_7", exp2_Para21_1="Para21_1", exp2_Para21_2="Para21_2", exp2_Para21_3="Para21_3", exp2_Para21_4="Para21_4", exp2_Para21_5="Para21_5", exp2_Para21_6="Para21_6", exp2_Para21_7="Para21_7", exp2_Para21_8="Para21_8",exp2_Para22_1="Para22_1", exp2_Para22_2="Para22_2", exp2_Para22_3="Para22_3", exp2_Para22_5="Para22_5", exp2_Para22_6="Para22_6", exp2_Para22_7="Para22_7", exp2_Para22_8="Para22_8");
fl <- list()
for ( i in names(nl) ) fl[[i]] <- n4
for ( sl in names(nl) ){
    cat("loading ",sl,"\n")
    x <- read.csv(paste("./",sl,".coutt",".csv",sep=""),sep="\t",header=TRUE)   ### path to processed data files (.coutt) files
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

PI1st1host <- givecells(CD45_2_number, "Para11")
PI1st2host <- givecells(CD45_2_number, "Para11", a=c(5,6,7,8))
PI2nd1host <- givecells(CD45_1_number, "Para12") ## order swapped by sorting operator donor top
PI2nd2host <- givecells(CD45_2_number, "Para12", a=c(5,6,7,8))
PI1st1don <- givecells(CD45_1_number, "Para11")
PI1st2don <- givecells(CD45_1_number, "Para11", a=c(5,6,7,8))
PI2nd1don <- givecells(CD45_2_number, "Para12") ## order swapped by sorting operator host down
PI2nd2don <- givecells(CD45_1_number, "Para12", a=c(5,6,7,8))      
PII1st1host <- givecells(CD45_2_number, "Para21")
PII1st2host <- givecells(CD45_2_number, "Para21", a=c(5,6,7,8))
PII2nd1host <- givecells(CD45_2_number, "Para22")
PII2nd2host <- givecells(CD45_2_number, "Para22", a=c(5,6,7,8))
PII1st1don <- givecells(CD45_1_number, "Para21")
PII1st2don <- givecells(CD45_1_number, "Para21", a=c(5,6,7,8))
PII2nd1don <- givecells(CD45_1_number, "Para22")
PII2nd2don <- givecells(CD45_1_number, "Para22", a=c(5,6,7,8))
      

PI2nd2host <- PI2nd2host[grep("_8_", PI2nd2host, invert=T)] ## missing sample
PI2nd2don <- PI2nd2don[grep("_8_", PI2nd2don, invert=T)]   ## missing sample
PII2nd1host <- PII2nd1host[grep("_4_", PII2nd1host, invert=T)] ## missing sample
PII2nd1don <- PII2nd1don[grep("_4_", PII2nd1don, invert=T)] ## missing sample
      
n <- names(data)
	
names(n) <- n
n[PI1st1host] <- sub("Para11","PI1st1host",PI1st1host)
n[PI1st2host] <- sub("Para11","PI1st2host",PI1st2host)
n[PI2nd1host] <- sub("Para12","PI2nd1host",PI2nd1host)
n[PI2nd2host] <- sub("Para12","PI2nd2host", PI2nd2host)
n[PI1st1don] <- sub("Para11","PI1st1don",PI1st1don)
n[PI1st2don] <- sub("Para11","PI1st2don",PI1st2don)
n[PI2nd1don] <- sub("Para12","PI2nd1don",PI2nd1don)
n[PI2nd2don] <- sub("Para12", "PI2nd2don", PI2nd2don)
n[PII1st1host] <- sub("Para21", "PII1st1host", PII1st1host)
n[PII1st2host] <- sub("Para21", "PII1st2host", PII1st2host)
n[PII2nd1host] <- sub("Para22", "PII2nd1host", PII2nd1host)
n[PII2nd2host] <- sub("Para22", "PII2nd2host", PII2nd2host)
n[PII1st1don] <- sub("Para21", "PII1st1don", PII1st1don)
n[PII1st2don] <- sub("Para21", "PII1st2don", PII1st2don)
n[PII2nd1don] <- sub("Para22", "PII2nd1don", PII2nd1don)
n[PII2nd2don] <- sub("Para22", "PII2nd2don", PII2nd2don)

names(data) <- as.character(n)
data <- data[,grep("Para", colnames(data), invert = T)]
saveRDS(data, "Parabiosis2.RData")
  

