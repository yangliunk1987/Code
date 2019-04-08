genotype_frequency <-
function(alt_freq){
    genotype_freq<-c(alt_freq^2,2*alt_freq*(1-alt_freq),(1-alt_freq)^2)
    return(genotype_freq)
}
