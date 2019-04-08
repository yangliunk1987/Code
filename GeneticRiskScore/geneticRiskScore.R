geneticRiskScore <-
function(snp.effect_freq,snp.effect_risk,snp.genotype){
    #if(snp.effect_freq | snp.effect_risk | snp.genotype )
#风险位点等位基因频率
snp.alt_freq<-snp.effect_freq
#风险位点OR值
snp.alt_risk<-snp.effect_risk
#风险位点样本检测基因型
snp.test_genotype<-snp.genotype

#根据提供的风险位点,计算hom/het/wild 三种基因型的相对风险值[]
snp.genotype_rr<-t(sapply(snp.alt_risk,genotype_relative_risk))

#获得位点样本检测基因型的对应风险值
snp.test_genotype_rr<-get_genotype_rr(snp.genotype_rr,snp.test_genotype)

#calculation population average risk
#计算人群平均相对风险
pop_avg<-c()
for(i in 1:nrow(snp.genotype_rr)){
    pop_avg<-c(pop_avg,average_poplation_risk(snp.genotype_rr[i,],snp.alt_freq[i]))
}

#计算每个位点相对人群平均风险值
#relative to population relative risk
snp.pop_rr<-snp.test_genotype_rr/pop_avg

#个体总的相对风险值
disease.risk<-prod(snp.pop_rr)
#最终结果
return(disease.risk)

}
