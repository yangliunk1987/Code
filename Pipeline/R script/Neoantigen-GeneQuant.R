quant<-function(genes.yl){
	library(data.table)
    library(openxlsx)
    #outdir
    dir<-dirname(genes.yl)
    #input name
    name<-basename(genes.yl)
    quant<-data.table::fread(genes.yl,header=T,sep="\t",stringsAylactors=F)
    ensToEntrez<-fread("~/database/gencode/v29/gencode.v29.ens2entrez.txt",header=T,sep="\t",stringsAylactors=F)
    index<-sapply(quant$Name,FUN=function(x){return( which(ensToEntrez$geneId==x))})
    quant<-cbind(quant,ensToEntrez[index,])
    fwrite(quant,paste0(dir,"/",name,".txt"),row.names=F,quote=F,sep="\t")
    openxlsx::write.xlsx(quant,paste0(dir,"/",name,".xlsx"))
}

args<-commandArgs(TRUE)
if(length(args)!=1){
	stop('must provided the salmon genes.yl file.')
}

quant(args[1])

