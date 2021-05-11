if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("karyoploteR")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

setwd("/Users/stefanolise/Documents/ROSETREES/TR019/COVERAGE_PROFILE/")

covHC <- read.csv("TR019.100K.cov.gz",sep = "\t",header = T )
covHC <- subset(covHC,X.chr == 1 | X.chr == 2 | X.chr == 3 | X.chr == 4 | X.chr == 5 | X.chr == 6 |
                X.chr == 7 | X.chr == 8 | X.chr == 9 | X.chr == 10 | X.chr == 11 | X.chr == 12 |
                X.chr == 13 | X.chr == 14 | X.chr == 15 | X.chr == 16 | X.chr == 17 | X.chr == 18 |
                X.chr == 19 | X.chr == 20 | X.chr == 21 | X.chr == 22 | X.chr == "X")

iGL <- 8
covHC <- covHC[covHC[,iGL] !=0,]
binSize <- 100000
ymin=0
ymax=3

iCol <- 5
normFact <- median(covHC[,iGL])/median(covHC[,iCol])
#normFact <- 1/median(covHC[,iCol])

kp <- plotKaryotype(chromosomes=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'),plot.type = 4,genome="hs37d5",cex=0.5)
for(i in c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)){
  chr <- i
  normFactChr <- normFact
  chr_name <- sprintf("%d",chr)
  if (i == 23){
    chr <- 'X'
    chr_name <- sprintf("X")
    normFactChr <- normFact/2
  }
  idx <- which(covHC[,1] == chr)
  covRatio <- normFactChr*covHC[idx,iCol]/covHC[idx,iGL]
#  covRatio <- normFact*covHC[idx,iCol]
  pos <- covHC[idx,2] + binSize/2
  kpPoints(kp, chr=chr_name, y=covRatio,x=pos,cex=0.4,ymin=ymin,ymax=ymax,plot.type=4)
  a1 <- rep(0,times=length(pos))
  kpLines(kp,chr=chr_name,y=a1,x=pos,cex=0.2,ymin=ymin,ymax=ymax,plot.type=4, col="red")
  finalPos <- max(pos)
  y1 <- seq(ymin, ymax, by = 1)
  x1 <- rep(finalPos,times=length(y1))
  kpLines(kp,chr=chr_name,y=y1,x=x1,cex=0.2,ymin=ymin,ymax=ymax,plot.type=4, col="red")
}
kpAxis(kp, ymin=ymin, ymax=ymax,numticks=11,cex=0.5,side=2)
iCol
colnames(covHC)[iCol]