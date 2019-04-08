library(GeneticRiskScore)
data<-read.table('risk_db.final.txt',sep="\t",stringsAsFactors = F,fileEncoding = "UTF-8")
disease<-unique(data[,11])
for(i in 1:length(disease)){
  i=27
  dd<-data[data[,11]==disease[i],]
  
  cat(disease[i],"\t")
  cat(nrow(dd),"\t")
  
  #风险位点等位基因频率
  snp.alt_freq<-dd[,10]
  
  #risk loci genotype frequency distribution
  snp.genotype_freq<-t(sapply(snp.alt_freq,genotype_frequency))
  
  #风险位点OR值
  snp.alt_risk<-dd[,9]
  
  #根据提供的风险位点,计算hom/het/wild 三种基因型的相对风险值[]
  snp.genotype_rr<-t(sapply(snp.alt_risk,genotype_relative_risk))
  #calculation population average risk
  #计算人群平均相对风险
  pop_avg<-c()
  for(i in 1:nrow(snp.genotype_rr)){
    pop_avg<-c(pop_avg,average_poplation_risk(snp.genotype_rr[i,],as.double(snp.alt_freq[i])))
  }
  #calculation genotype risk to population average risk
  snp.genotype_pop.rr<-snp.genotype_rr/pop_avg
  
  if(nrow(dd)>=2&nrow(dd)<=19){
    #get all genotype combination
    b<-snp.genotype_pop.rr[1,]
    len<-nrow(snp.genotype_pop.rr)
    for(i in 2:len){
      b<-as.vector(outer(b,c(snp.genotype_pop.rr[i,]),'*'))
      #b<-outer(b,c(snp.genotype_pop.rr[i,]),*)
    }
    cat(min(b),"\t",max(b),"\t",mean(b),"\t",sd(b),"\n")
  }else if(nrow(dd)>=20){
    rr<-matrix(nrow=nrow(snp.genotype_pop.rr),ncol=2)
    for(j in 1:nrow(snp.genotype_freq)){
      #sort genotype frequency
      s<-sort(snp.genotype_freq[j,],decreasing = T)
      index<-which(snp.genotype_freq[j,] %in% s[1:2])
      #get the max and the second
      if(s[1]>=0.8){
        rr[j,1]<-snp.genotype_pop.rr[j,index[1]]
        
      }else{
        
        rr[j,]<-snp.genotype_pop.rr[j,index]
      }
     
      
    }
    b<-rr[1,]
    len<-nrow(rr)
    for(i in 2:len){
      b<-as.vector(outer(b,c(rr[i,]),'*'))
      #b<-outer(b,c(snp.genotype_pop.rr[i,]),*)
    }
    cat(min(b),"\t",max(b),"\t",mean(b),"\t",sd(b),"\n")
  }
  
}


