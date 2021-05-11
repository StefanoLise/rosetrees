if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("karyoploteR")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

setwd("/Users/stefanolise/Documents/ROSETREES/TR064/COVERAGE_PROFILE/")

cov_gl <- read.csv("TR064_GL.1M.seg.gz",sep = "\t",header = T )
cov_bl <- read.csv("TR064_BL.1M.seg.gz",sep = "\t",header = T )
cov_bl <- cov_bl[cov_gl$count !=0,]
cov_gl <- cov_gl[cov_gl$count !=0,]

binSize <- 1000000
ymin=0
ymax=3


normFact <- median(cov_gl$count)/median(cov_bl$count)
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
  idx <- which(cov_gl[,2] == chr)
  covRatio <- normFactChr*cov_bl[idx,5]/cov_gl[idx,5]
#  covRatio <- normFact*covHC[idx,iCol]
  pos <- cov_gl[idx,3] + binSize/2
  kpPoints(kp, chr=chr_name, y=covRatio,x=pos,cex=0.4,ymin=ymin,ymax=ymax,plot.type=4)
  a1 <- rep(0,times=length(pos))
  kpLines(kp,chr=chr_name,y=a1,x=pos,cex=0.2,ymin=ymin,ymax=ymax,plot.type=4, col="red")
  finalPos <- max(pos)
  y1 <- seq(ymin, ymax, by = 1)
  x1 <- rep(finalPos,times=length(y1))
  kpLines(kp,chr=chr_name,y=y1,x=x1,cex=0.2,ymin=ymin,ymax=ymax,plot.type=4, col="red")
}
kpAxis(kp, ymin=ymin, ymax=ymax,numticks=11,cex=0.5,side=2)