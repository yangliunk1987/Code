average_poplation_risk <-
function(genotype_relative_risk,allele_frequency){
    avg_risk<-(allele_frequency^2*genotype_relative_risk[1])+(2*allele_frequency*(1-allele_frequency)*genotype_relative_risk[2])+((1-allele_frequency)^2*genotype_relative_risk[3])
    return(avg_risk)
}
