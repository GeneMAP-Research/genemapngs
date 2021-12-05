def getBamFastq() {
    if(params.dataFormat == "BAM") {
        return channel.fromFilePairs( params.bamDir + "*.bam", size: 1 )
                  .map { bamName, bamFile -> tuple(bamName, bamFile.first(), bamFile.first()) }
                  .ifEmpty{ error "\nERROR Oops... no BAM files were found in ${params.bamDir}\n" }
    }

    else if(params.dataFormat == "FASTQ") {
        return channel.fromFilePairs( params.bamDir + "*.f*q*", size: 2 )
                  .map { fastqName, fastqFile -> tuple(fastName, fastqFile.first()) }
                  .ifEmpty{ error "\nERROR Oops... no FASTQ files were found in ${params.bamDir}\n" }
    }

    else {  
        error "\nERROR: Please set the data format [BAM/FASTQ] in the config file: e.g. dataFormat = 'BAM' \n"
    }
}

process getFastqQualityReports() {
    format = "${params.dataFormat}"
    label 'fastqc'
    label 'bigMemory'
    input:
        tuple \
            val(dataName), \
            path(fileOne), \
            path(fileTwo)
    output:   
        val "${dataName}"
    script:   
    if( format == "FASTQ")
        """
        mkdir -p ${params.outputDir}/fastqc
        fastqc \
            -f fastq \
            -o ${params.outputDir}/fastqc/ \
            -t ${task.cpus} \
            ${fileOne} \
            ${fileTwo}
        """
    else if( format == "BAM" )
        """
        mkdir -p ${params.outputDir}/fastqc
        fastqc \
            -f bam \
            -o ${params.outputDir}/fastqc/ \
            -t ${task.cpus} \
            ${fileOne}
        """
}

/*
process getFastqQualityReports() {
    tag "processing ${fastqName}"
    label 'fastqc'
    label 'bigMemory'
    input:
        tuple \
            val(fastqName), \
            path(readOne), \
            path(readTwo)
    output:
        val "${fastqName}"
    script:
        """
        mkdir -p ${params.outputDir}/fastqc
        fastqc \
            -f fastq \
            -o ${params.outputDir}/fastqc/ \
            -t ${task.cpus} \
            ${readOne} \
            ${readTwo}
        """
}
*/

process getMultiQcFastqReports() {
    tag "Writing MULTIQC Report "
    label 'multiqc'
    label 'bigMemory'
    input:
        val(fastqName)
    script:
        """
        mkdir -p ${params.outputDir}/multiqc
        multiqc \
            ${params.outputDir}/fastqc/ \
            -o ${params.outputDir}/multiqc/
        """
}

