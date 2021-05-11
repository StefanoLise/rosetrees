setwd("/Users/stefanolise/Documents/ROSETREES/TR019/COVERAGE_PROFILE/")

covHC <- read.csv("TR064.100k.cov.gz",sep = "\t",header = T )
covHC <- subset(covHC,X.chr == 1 | X.chr == 2 | X.chr == 3 | X.chr == 4 | X.chr == 5 | X.chr == 6 |
                  X.chr == 7 | X.chr == 8 | X.chr == 9 | X.chr == 10 | X.chr == 11 | X.chr == 12 |
                  X.chr == 13 | X.chr == 14 | X.chr == 15 | X.chr == 16 | X.chr == 17 | X.chr == 18 |
                  X.chr == 19 | X.chr == 20 | X.chr == 21 | X.chr == 22 | X.chr == "X")

iGL <- 8
covHC <- covHC[covHC[,iGL] !=0,]
binSize <- 100000

iCol <- 5
chr <- 13
normFact <- median(covHC[,iGL])/median(covHC[,iCol])
normFactChr <- normFact
if (chr == 'X'){
  normFactChr <- normFact/2
}
idx <- which(covHC[,1] == chr)
covRatio <- normFactChr*covHC[idx,iCol]/covHC[idx,iGL]
pos <- covHC[idx,2] + binSize/2
plot(pos,covRatio,
     #     xlim=c(60e6,70e6),
     ylim=c(0.,2))
x1<-90e6
lines(c(x1,x1),c(0.,2.),lty=1,lwd=2)
y1<-0.70
lines(c(0e6,250e6),c(y1,y1),lty=1,lwd=2)


lines(c(210e6,210e6),c(0.,2.),lty=1,lwd=2)
