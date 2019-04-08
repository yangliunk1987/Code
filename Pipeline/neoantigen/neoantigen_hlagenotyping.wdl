workflow neoantigen_halgenotyping{
	String name
	String f1
    String f2
    String inputDir
    String outputDir

    call opitype{
        input:
            name=name,
            f1=f1,
            f2=f2,
            inputDir=inputDir,
            outputDir=outputDir
    }
}



task opitype{
	String name
	String f1
    String f2
    String inputDir
    String outputDir


    command <<<
        docker run \
        -v ${inputDir}:/data/ \
        -v ${outputDir}:/output/ \
        -t docker.io/fred2/optitype \
        --input /data/${f1} /data/${f2} \
        --dna \
        --outdir /output/ \
        --prefix ${name} \
        --verbose
    >>>

    output{
    	File coverage_plot = "${outputDir}/${name}_coverage_plot.pdf"
    	File result = "${outputDir}/${name}_result.tsv"
    }

}