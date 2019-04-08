library(GeneticRiskScore)
args <- commandArgs(T)
if(length(args) < 1){
    cat("Argument:tab format snp annotation file\n")
    quit('no')
}
name<-gsub(pattern=".txt",replacement="",x=args[1])
data<-read.table(args[1],header=F,sep="\t",stringsAsFactors=FALSE)
data<-data[data[,10]!='a',]
#获得所有疾病列表
disease.grs<-data.frame(disease=unique(data[,11]),geneticRiskScore=c(rep(0,length(unique(data[,11])))))
#计算每个疾病的风险值
for(i in 1:length(disease.grs$disease)){
    print(disease.grs$disease[i])
    disease.data<-data[data[,11]==disease.grs$disease[i],]
    grs<-geneticRiskScore(snp.effect_freq=as.double(disease.data[,10]),snp.effect_risk=as.double(disease.data[,9]),snp.genotype=disease.data[,13])
    disease.grs[i,c("geneticRiskScore")]<-grs
}

#输出结果
write.table(disease.grs,paste0(name,".grs.txt"),row.names=F,sep="\t",quote=F)
