workflow neoantigen_quantify{
	String name
	String f1
    String f2
    String outputDir

    call salmon{
    	input:
    	    name=name,
    	    f1=f1,
    	    f2=f2,
    	    outputDir=outputDir
    }

    call gene_quant{
    	input:
    	    name=name,
    	    outputDir=outputDir,
    	    genes_yl=salmon.gene
    }

}


task salmon{
	String name
	String f1
    String f2
    String outputDir

    command <<<
        /home/yangliu/biosoft/rnaseq/salmon-0.11.3-linux_x86_64/bin/salmon quant \
        --libType A \
        --seqBias \
        --gcBias \
        --threads 10 \
        --index /home/yangliu/database/gencode/v29/genecode_v29/ \
        --geneMap /home/yangliu/database/gencode/v29/gencode.v29.annotation.gtf \
        --mates1 ${f1} \
        --mates2 ${f2} \
        --output ${outputDir}/${name}
    >>>

    output{
    	String transcript="${outputDir}/${name}/quant.yl"
    	String gene="${outputDir}/${name}/quant.genes.yl"
    }
}


task gene_quant{
	String genes_yl
	String name
	String outputDir

	command <<<
	    Rscript /home/yangliu/biosoft/easybio/rscripts/Neoantigen-GeneQuant.R ${genes_yl}
	>>>

	output{
		String txt = "${outputDir}/${name}/quant.genes.yl.txt"
		String xlsx = "${outputDir}/${name}/quant.genes.yl.xlsx"
	}
}
