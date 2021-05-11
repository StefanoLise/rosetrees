setwd("/Users/stefanolise/Documents/ROSETREES/TR019/HAP_BAF")
hapBAF.GL <- read.csv("TR019_GL.HRC.SHAPEIT4.100snp.tsv",sep = "\t",header = T)
hapBAF.BL <- read.csv("TR019_BL.HRC.SHAPEIT4.400snp.tsv",sep = "\t",header = T)

hapBAF <- hapBAF.BL
hapBAF$bin_pos <- (hapBAF$bin_start + hapBAF$bin_end)/2

chr=22
idx <- which(hapBAF[,1] == chr & hapBAF$n_snp > 75)
plot(hapBAF$bin_pos[idx],hapBAF$af_l[idx],pch=16,col=2,
     xlab=paste("chr",chr," (in Mb)",sep=""),ylab="BAF",
#     xlim=c(80e6,130e6),
     ylim=c(0.2,0.8),axes =T,cex.lab=1.3
)
points(hapBAF$bin_pos[idx],hapBAF$af_r[idx],pch=16,col=4)

ticks=c(0,20e6,40e6,60e6,80e6,100e6,120e6,140e6)
labels=c(0,20,40,60,80,100,120,140)
axis(1,at=ticks,labels=labels)

x1<-8e6
lines(c(x1,x1),c(0.,1.),lty=1,lwd=2)
y1<-0.48
lines(c(0e6,250e6),c(y1,y1),lty=1,lwd=2)

