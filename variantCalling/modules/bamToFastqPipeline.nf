def getBamFiles() {
    return channel.fromFilePairs( params.bamDir + "*.bam", size: 1 )
                  .map { bamName, bamFile -> tuple(bamName, bamFile.first()) }
}

process sortBamByName() {
    tag "processing ${bamName}"
    label 'samtools'
    label 'bigMemory'
    input:
        tuple \
            val(bamName), \
            path(bamFile) 
    output:
        publishDir path: "${params.outputDir}/fastq/"
        tuple \
            val(bamName), \
            path("${bamName}.sortedByName.bam")
    script:
        """
        samtools \
            sort \
            --reference ${params.fastaRef} \
            -O BAM \
            --threads ${task.cpus} \
            -o "${bamName}.sortedByName.bam" \
            -n \
            ${bamFile}
        """
}

process convertBamToFastq() {
    tag "processing ${bamName}"
    label 'samtools'
    label 'mediumMemory'
    input:
        tuple \
            val(bamName), \
            path(bamFile)
    output:
        publishDir path: "${params.outputDir}/fastq/", mode: 'copy'
        tuple \
            val(bamName), \
            path("${bamName}_R1.fq.gz"), \
            path("${bamName}_R2.fq.gz")
    script:
        """
        samtools \
            fastq \
            -1 "${bamName}_R1.fq.gz" \
            -2 "${bamName}_R2.fq.gz" \
            -0 /dev/null \
            -s /dev/null \
            --threads ${task.cpus} \
            -n \
            ${bamFile}
        """
}
