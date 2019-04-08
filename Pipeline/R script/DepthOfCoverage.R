args <- commandArgs(TRUE)
if(length(args)!=2){
  stop('must provide sample and depth input\nusage:Rscript ReadDepth.R <sample> <depth_file>')
}
sample<-args[1]
input<-args[2]
data<-read.table(input,header=T,sep="\t")
cov<-data$mean_coverage
max<-max(cov)
min<-min(cov)
mean<-mean(cov)
depth<-c()
for(i in 1:(3*mean)){
  depth[i]<-length(cov[cov>i])/length(cov)
}
png(paste0(sample,".png"),width=1024,height=768)
plot(depth,lwd=3,type="l",col="red",xlab="Base Read Depth",ylab="Cumulative Bases",main=sample)
#txt_x<- which(depth==100)[1]
txt1<-paste("Mean Depth",round(mean),sep="    ")
txt2<-paste(">=1X ",paste0(round(depth[1],4)*100,"%"),sep="           ")
txt3<-paste(">=100X ",paste0(round(depth[100],4)*100,"%"),sep="            ")
text(x=100,y=0.5,labels=paste(txt1,txt2,txt3,sep="\n"))
abline(v=1,col="firebrick",lwd=1,lty=2)
abline(v=100,col="blue",lwd=1,lty=2)
grid()

#pch=c(15,16,17) shape control
legend(x="topright",c("Depth","1X","100X"),col=c("red","firebrick","blue"),text.col=c("red","firebrick","blue"),lty=c(2,2,2))
dev.off()