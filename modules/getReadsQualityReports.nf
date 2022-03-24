process getBamQualityReports() {
    tag "processing ${bamName}"
    label 'fastqc'
    label 'fastqc_mem'
    input:
        tuple \
            val(bamName), \
            path(bamFile)
    output:
        publishDir path: "${params.outputDir}/fastqc/", mode: 'copy' 
        path "${bamName}*"
    script:   
        """
        fastqc \
            -f bam \
            -o . \
            -t ${task.cpus} \
            ${bamFile}
        """
}

process getFastqQualityReports() {
    tag "processing ${fastqName}"
    label 'fastqc'
    label 'fastqcMem'
    input:
        tuple \
            val(fastqName), \
            path(fastqReads)
    output:
        publishDir path: "${params.outputDir}/fastqc/", mode: 'copy'
        path "${fastqName}*"
    script:
        (readOne, readTwo) = fastqReads
        """
        fastqc \
            -f fastq \
            -o . \
            -t ${task.cpus} \
            ${readOne} \
            ${readTwo}
        """
}

process getMultiQcFastqReports() {
    tag "Writing MULTIQC Report"
    label 'multiqc'
    label 'multiqcMem'
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

