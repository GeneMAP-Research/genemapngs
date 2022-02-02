process getReadsQualityReports() {
    format = "${params.inputFileType}"
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

