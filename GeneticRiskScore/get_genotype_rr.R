get_genotype_rr <-
function(snp.genotype_rr,snp.test_genotype){
    if(nrow(snp.genotype_rr)!=length(snp.test_genotype)){
    stop("snp numbers  is not equal between [genotype relative risk] and [test gentoype] ")
    }
    c<-c()
    for(i in 1:length(snp.test_genotype)){
    if(snp.test_genotype[i]=="het"){
            c[i]<-snp.genotype_rr[i,2]
    }else if(snp.test_genotype[i]=="hom"){
    c[i]<-snp.genotype_rr[i,1]
    }else{
    c[i]<-snp.genotype_rr[i,3]
    }
    }
    return(c)
}
