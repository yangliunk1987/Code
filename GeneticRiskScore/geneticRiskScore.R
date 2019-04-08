geneticRiskScore <-
function(snp.effect_freq,snp.effect_risk,snp.genotype){
    #if(snp.effect_freq | snp.effect_risk | snp.genotype )
#����λ���λ����Ƶ��
snp.alt_freq<-snp.effect_freq
#����λ��ORֵ
snp.alt_risk<-snp.effect_risk
#����λ��������������
snp.test_genotype<-snp.genotype

#�����ṩ�ķ���λ��,����hom/het/wild ���ֻ����͵���Է���ֵ[]
snp.genotype_rr<-t(sapply(snp.alt_risk,genotype_relative_risk))

#���λ�������������͵Ķ�Ӧ����ֵ
snp.test_genotype_rr<-get_genotype_rr(snp.genotype_rr,snp.test_genotype)

#calculation population average risk
#������Ⱥƽ����Է���
pop_avg<-c()
for(i in 1:nrow(snp.genotype_rr)){
    pop_avg<-c(pop_avg,average_poplation_risk(snp.genotype_rr[i,],snp.alt_freq[i]))
}

#����ÿ��λ�������Ⱥƽ������ֵ
#relative to population relative risk
snp.pop_rr<-snp.test_genotype_rr/pop_avg

#�����ܵ���Է���ֵ
disease.risk<-prod(snp.pop_rr)
#���ս��
return(disease.risk)

}