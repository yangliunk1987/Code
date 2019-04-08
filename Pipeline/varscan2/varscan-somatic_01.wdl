workflow varscan2_somatic_01{
    String name
    String normal
    String tumor
    String inputDir
    String outputDir
    String ref
    String bed
    String accessBed
    String annotateGene

    parameter_meta{
        name: "sample name" 
        normal: "normal file prefix , [if file name is 'sample_B_R1.fq.gz' ,then normal is 'sample_B']"
        tumor: "tumor file prefix , [if file name is 'sample_T_R1.fq.gz' ,then tumor is 'sample_T']"
        inputDir: "directory contains raw fastq files"
        outputDir: "the result directory"
        ref: "reference fasta file,and the directory must contains bwa index files"
        bed: "the target region bed file"
        accessBed: "cnvkit access bed file"
        annotateGene: "cnvkit gene annotation file"
    }
    

    #创建目录
    call create_dir{
        input:
            outputDir=outputDir
    }

    #分别处理正常与组织样本[质控->比对->去除PCR重复]
    scatter (name in [normal,tumor]){
        call qc{
                input:
                    name=name,
                    inputDir=inputDir,
                    outputDir=outputDir
        }

        call alignment_bwa{
            input:
                name=name,
                outputDir=outputDir,
                ref=ref,
                read1=qc.outputFile[0],
                read2=qc.outputFile[1]
        }

        call markduplicates{
            input:
                name=name,
                outputDir=outputDir,
                ref=ref,
                sortedBam=alignment_bwa.sorted
        }
    }

    #完成前处理，开始进行突变检测
    #call mpileup_varscan2{
    #    input:
    #        name=name,
    #        normal=normal,
    #        tumor=tumor,
    #        ref=ref,
    #        bed=bed,
    #        outputDir=outputDir,
    #        rmdupBam=markduplicates.rmdup
    #}

    call generater_mpileup{
        input:
            name=name,
            normal=normal,
            tumor=tumor,
            ref=ref,
            bed=bed,
            outputDir=outputDir,
            rmdupBam=markduplicates.rmdup    

    }

    call mutation_calling{
        input:
            name=name,
            outputDir=outputDir,
            pileup=generater_mpileup.pileup

    }

    #annotate mutations
    call annovar{
        input:
            name=name,
            outputDir=outputDir,
            vcf= mutation_calling.vcf
    }

    #call copy number alteration with cnvkit
    call cnvkit{
        input:
            name=name,
            normal=normal,
            tumor=tumor,
            ref=ref,
            bed=bed,
            outputDir=outputDir,
            rmdupBam=markduplicates.rmdup,
            accessBed=accessBed,
            annotateGene=annotateGene
    }

    

}


#create project directory
task create_dir{
    String outputDir

    command <<<
        #创建目录
        if [ ! -d ${outputDir} ];then
            mkdir ${outputDir}
        fi

        #创建qc目录
        if [ ! -d ${outputDir}/qc ];then
            mkdir ${outputDir}/qc
        fi

        #创建alignment目录
        if [ ! -d ${outputDir}/alignment ];then
            mkdir ${outputDir}/alignment
        fi

        #创建mutation目录
        if [ ! -d ${outputDir}/mutation ];then
            mkdir ${outputDir}/mutation
        fi

        #创建cnv目录
        if [ ! -d ${outputDir}/cnvkit ];then
            mkdir ${outputDir}/cnvkit
        fi 
    >>>

}

#generator raw fastq to clean fastq
task qc{
    String name
    String inputDir
    String outputDir

    command <<<
        echo processing raw reads with fastp
        fastp -i ${inputDir}/${name}_R1.fq.gz -o ${outputDir}/qc/${name}_clean_R1.fq.gz \
        -I ${inputDir}/${name}_R2.fq.gz -O ${outputDir}/qc/${name}_clean_R2.fq.gz \
        -w 10 \
        --adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
        --adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
        -j ${outputDir}/qc/${name}.json \
        -h ${outputDir}/qc/${name}.html --report_title ${name}
    >>>

    output{
        Array[String] outputFile = [
        "${outputDir}/qc/${name}_clean_R1.fq.gz ",
        "${outputDir}/qc/${name}_clean_R2.fq.gz",
        "${outputDir}/qc/${name}.json",
        "${outputDir}/qc/${name}.html"
        ]
    }
}

#alignment clean fastq to reference
task alignment_bwa{
      String name
      String ref
      String outputDir
      String read1
      String read2
      command<<<
          bwa mem -M -t 10 ${ref} \
          ${outputDir}/qc/${name}_clean_R1.fq.gz ${outputDir}/qc/${name}_clean_R2.fq.gz | \
          samtools view -@ 2 -bh  -o - | samtools sort -@ 2  -o ${outputDir}/alignment/${name}.sorted.bam

          samtools index ${outputDir}/alignment/${name}.sorted.bam
      >>>

      output{
          String sorted = "${outputDir}/alignment/${name}.sorted.bam"
      }
}

#remove PCR duplicates 
task markduplicates{
    String name
    String ref
    String outputDir
    String sortedBam

    command<<<
        java -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx4G  -jar /home/biosoftware/install_pkg/picard-2.9.2/picard/build/libs/picard.jar MarkDuplicates \
        I=${outputDir}/alignment/${name}.sorted.bam \
        O=${outputDir}/alignment/${name}.rmdup.bam  \
        CREATE_INDEX=true \
        M=${outputDir}/alignment/${name}.rmdup.metrics.txt  \
        R=${ref}
    >>>

    output{
              String  rmdup = "${outputDir}/alignment/${name}.rmdup.bam"
    }
}


# generater mpileup file
task generater_mpileup{
    String name
    String normal
    String tumor
    String ref
    String bed
    String outputDir
    Array[String] rmdupBam

    command<<<
        samtools mpileup -Bq 20 -Q 20  -f ${ref} -l ${bed} \
        ${outputDir}/alignment/${normal}.rmdup.bam \
        ${outputDir}/alignment/${tumor}.rmdup.bam | gzip > ${outputDir}/alignment/${name}.pileup.gz
    >>>

    output{
        String pileup = "${outputDir}/alignment/${name}.pileup.gz"
    }

}

# call mutation with varscan2
task mutation_calling{
    String name
    String outputDir
    String pileup

    command<<<
        gzip ${outputDir}/alignment/${name}.pileup.gz | java -jar /home/biosoftware/varscan/VarScan.v2.4.3.jar somatic \
        --mpileup --output-snp ${outputDir}/mutation/${name}.snp.vcf \
        --output-indel ${outputDir}/mutation/${name}.indel.vcf \
        --min-var-freq 0.01 \
        --min-freq-for-hom 0.9 \
        --strand-filter 1 \
        --output-vcf 1 \
        --min-avg-qual 20  \
        --min-coverage-normal 10 \
        --min-coverage-tumor 30 --min-reads2 6 --min-strands2 2 --data-ratio 1

        java -jar /home/biosoftware/varscan/VarScan.v2.4.3.jar  processSomatic \
        ${outputDir}/mutation/${name}.snp.vcf \
        --min-tumor-freq 0.01 \
        --max-normal-freq 0.01 \
        --p-value 0.05

        java -jar /home/biosoftware/varscan/VarScan.v2.4.3.jar  processSomatic \
        ${outputDir}/mutation/${name}.indel.vcf \
        --min-tumor-freq 0.01 \
        --max-normal-freq 0.01 \
        --p-value 0.05
    >>>

    output{
        Array[String] vcf = [
        "${outputDir}/mutation/${name}.snp.Somatic.hc.vcf",
        "${outputDir}/mutation/${name}.indel.Somatic.hc.vcf",
        "${outputDir}/mutation/${name}.snp.Germline.hc.vcf",
        "${outputDir}/mutation/${name}.indel.Germline.hc.vcf",
        "${outputDir}/mutation/${name}.snp.LOH.hc.vcf",
        "${outputDir}/mutation/${name}.indel.LOH.hc.vcf"
        ]
    }

}

#convert bam to mpileup and call mutation with varscan2
task mpileup_varscan2{
    String name
    String normal
    String tumor
    String ref
    String bed
    String outputDir
    Array[String] rmdupBam
    command<<<
        samtools mpileup -Bq 20 -Q 20  -f ${ref} -l ${bed} \
        ${outputDir}/alignment/${normal}.rmdup.bam \
        ${outputDir}/alignment/${tumor}.rmdup.bam | java -jar /home/biosoftware/varscan/VarScan.v2.4.3.jar somatic \
        --mpileup --output-snp ${outputDir}/mutation/${name}.snp.vcf \
        --output-indel ${outputDir}/mutation/${name}.indel.vcf \
        --min-var-freq 0.01 \
        --min-freq-for-hom 0.9 \
        --strand-filter 1 \
        --output-vcf 1 \
        --min-avg-qual 20  \
        --min-coverage-normal 10 \
        --min-coverage-tumor 30 --min-reads2 6 --min-strands2 2 --data-ratio 1

        java -jar /home/biosoftware/varscan/VarScan.v2.4.3.jar  processSomatic \
        ${outputDir}/mutation/${name}.snp.vcf \
        --min-tumor-freq 0.01 \
        --max-normal-freq 0.01 \
        --p-value 0.05

        java -jar /home/biosoftware/varscan/VarScan.v2.4.3.jar  processSomatic \
        ${outputDir}/mutation/${name}.indel.vcf \
        --min-tumor-freq 0.01 \
        --max-normal-freq 0.01 \
        --p-value 0.05
    >>>

    output{
    	Array[String] vcf = [
    	"${outputDir}/mutation/${name}.snp.Somatic.hc.vcf",
    	"${outputDir}/mutation/${name}.indel.Somatic.hc.vcf",
    	"${outputDir}/mutation/${name}.snp.Germline.hc.vcf",
    	"${outputDir}/mutation/${name}.indel.Germline.hc.vcf",
    	"${outputDir}/mutation/${name}.snp.LOH.hc.vcf",
    	"${outputDir}/mutation/${name}.indel.LOH.hc.vcf"
    	]
    }
}


#annotation mutation with annovar
task annovar{
    String name
    String outputDir
    Array[String] vcf
    command<<<
        perl /home/database/annovar/table_annovar.pl \
        ${outputDir}/mutation/${name}.snp.Somatic.hc.vcf \
        /home/database/annovar/humandb/ \
        -buildver hg19 -nastring . -vcfinput -remove -otherinfo \
        -protocol refGene,cytoBand,esp6500siv2_all,1000g2015aug_all,snp138,cosmic81,gnomad_genome,exac03 \
        -operation g,r,f,f,f,f,f,f

        perl /home/database/annovar/table_annovar.pl \
        ${outputDir}/mutation/${name}.indel.Somatic.hc.vcf \
        /home/database/annovar/humandb/ \
        -buildver hg19 -nastring . -vcfinput -remove -otherinfo \
        -protocol refGene,cytoBand,esp6500siv2_all,1000g2015aug_all,snp138,cosmic81,gnomad_genome,exac03 \
        -operation g,r,f,f,f,f,f,f

        perl /home/database/annovar/table_annovar.pl \
        ${outputDir}/mutation/${name}.snp.Germline.hc.vcf \
        /home/database/annovar/humandb/ \
        -buildver hg19 -nastring . -vcfinput -remove -otherinfo \
        -protocol refGene,cytoBand,esp6500siv2_all,1000g2015aug_all,snp138,cosmic81,gnomad_genome,exac03 \
        -operation g,r,f,f,f,f,f,f

        perl /home/database/annovar/table_annovar.pl \
        ${outputDir}/mutation/${name}.indel.Germline.hc.vcf \
        /home/database/annovar/humandb/ \
        -buildver hg19 -nastring . -vcfinput -remove -otherinfo \
        -protocol refGene,cytoBand,esp6500siv2_all,1000g2015aug_all,snp138,cosmic81,gnomad_genome,exac03 \
        -operation g,r,f,f,f,f,f,f

        perl /home/database/annovar/table_annovar.pl \
        ${outputDir}/mutation/${name}.snp.LOH.hc.vcf \
        /home/database/annovar/humandb/ \
        -buildver hg19 -nastring . -vcfinput -remove -otherinfo \
        -protocol refGene,cytoBand,esp6500siv2_all,1000g2015aug_all,snp138,cosmic81,gnomad_genome,exac03 \
        -operation g,r,f,f,f,f,f,f

        perl /home/database/annovar/table_annovar.pl \
        ${outputDir}/mutation/${name}.indel.LOH.hc.vcf \
        /home/database/annovar/humandb/ \
        -buildver hg19 -nastring . -vcfinput -remove -otherinfo \
        -protocol refGene,cytoBand,esp6500siv2_all,1000g2015aug_all,snp138,cosmic81,gnomad_genome,exac03 \
        -operation g,r,f,f,f,f,f,f
    >>>
}


task updateVcfSequenceDictionary{
    String snpVcf
    String? indelVcf
    String sequenceDictionary
    String outputDir
    String name

    command <<<
        java -jar /home/biosoftware/install_pkg/picard-2.9.2/picard/build/libs/picard.jar \
        UpdateVcfSequenceDictionary \
        I=${outputDir}/mutation/${name}.indel.vcf \
        O=${outputDir}/mutation/ \
        SD=${sequenceDictionary}
    >>>
}

#call copy number alteration
task cnvkit{
    String name
    String ref
    String bed
    String tumor
    String normal
    String outputDir
    String accessBed
    String annotateGene
    #for task chain
    Array[String] rmdupBam

    command <<<
        echo run cnvkit batch to processing cnv calling
        cnvkit.py batch \
        ${outputDir}/alignment/${tumor}.rmdup.bam \
        --normal ${outputDir}/alignment/${normal}.rmdup.bam \
        --targets ${bed} \
        --fasta ${ref} \
        --access ${accessBed} \
        --output-reference ${outputDir}/cnvkit/${normal}_reference.cnn \
        --annotate ${annotateGene} \
        --drop-low-coverage  --scatter --diagram --output-dir ${outputDir}/cnvkit
        
        echo run cnvkit genemetrics to generater gene level cnv
        cnvkit.py genemetrics \
        ${outputDir}/cnvkit/${tumor}.rmdup.cnr \
        -s ${outputDir}/cnvkit/${tumor}.rmdup.cns \
        -t 0.4 \
        -m 5 \
        -o ${outputDir}/cnvkit/${tumor}.GeneMetrics.txt
    >>>
}