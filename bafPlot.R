#library(Hmisc)
pat_id <-'TR064'
work_dir <-  paste0("/Users/stefanolise/Documents/ROSETREES/",pat_id,"/HAP_BAF")
sample_gl <- paste0(pat_id,"_GL.HRC.SHAPEIT4.100snp")
sample_bl <- paste0(pat_id,"_BL.HRC.SHAPEIT4.400snp")
out_file <-  sample_bl
setwd(work_dir)
hapBAF.GL <- read.csv(paste0(sample_gl,".tsv"),sep = "\t",header = T)
hapBAF.BL <- read.csv(paste0(sample_bl,".tsv"),sep = "\t",header = T)

hapBAF <- hapBAF.BL
hapBAF$bin_pos <- (hapBAF$bin_start + hapBAF$bin_end)/(2)


pdf(paste0("PLOTS/",out_file,".pdf"))
for(i in c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)){
  chr <- i
#  out_file <- paste0('chr',chr)
#  print(out_file)
#  pdf(paste("PLOTS/",out_file,".pdf"))
  idx <- which(hapBAF[,1] == chr & hapBAF$n_snp > 75)
  plot(hapBAF$bin_pos[idx],hapBAF$af_l[idx],pch=16,col=2,
       xlab=paste("chr",chr," (in Mb)",sep=""),ylab="BAF",
#       xlim=c(80e6,130e6),
      ylim=c(0.,1),
      xaxt="n",
      cex.lab=1.3
  )
  points(hapBAF$bin_pos[idx],hapBAF$af_r[idx],pch=16,col=4)

  ticks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260)
  labels=c(0,'',20,'',40,'',60,'',80,'',100,'',120,'',140,'',160,'',180,'',200,'',220,'',240,'',260)
  if (chr < 4){
    ticks=c(0,20,40,60,80,100,120,140,160,180,200,220,240,260)
    labels=c(0,'',40,'',80,'',120,'',160,'',200,'',240,'')
  }
  ticks <- ticks*1e6
  axis(1,at=ticks,labels=labels)
#  minor.tick(nx=2, ny=1, tick.ratio=0.5)
#  dev.off()
}
dev.off()

x1<-8e6
lines(c(x1,x1),c(0.,1.),lty=1,lwd=2)
y1<-0.48
lines(c(0e6,250e6),c(y1,y1),lty=1,lwd=2)

