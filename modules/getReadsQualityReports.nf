def getNpAdapter() {
    return channel.fromPath(projectDir + '/adapters/NexteraPE-PE.fa')
}
def getT3uAdapter() {
    return channel.fromPath(projectDir + '/adapters/TruSeq3-PE-2.fa')
}
def getT2pAdapter() {
    return channel.fromPath(projectDir + '/adapters/TruSeq2-PE.fa')
}
def getT3pAdapter() {
    return channel.fromPath(projectDir + '/adapters/TruSeq3-PE.fa')
}
def getT2sAdapter() {
    return channel.fromPath(projectDir + '/adapters/TruSeq2-SE.fa')
}
def getT3sAdapter() {
    return channel.fromPath(projectDir + '/adapters/TruSeq3-SE.fa')
}

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
    //////////////////debug true
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

process trimgalore() {
    tag "processing ${fastqName}"
    label 'trimgalore'
    label 'fastqcMem'
    input:
        tuple \
            val(fastqName), \
            path(fastqReads)
    output:
        publishDir path: "${params.outputDir}/trimmedreads/", mode: 'copy'
        path("*_val_R{1,2}.fastq.gz")
    script:
        (readOne, readTwo) = fastqReads
        """
        trim_galore \
            --paired \
            --clip_R1 ${params.headcrop} \
            --clip_R2 ${params.headcrop} \
            --three_prime_clip_R1 ${params.crop} \
            --three_prime_clip_R2 ${params.crop} \
            --basename ${fastqName} \
            --cores ${task.cpus} \
            -o . \
            ${readOne} \
            ${readTwo} 

        mv ${fastqName}_val_1.fq.gz ${fastqName}_val_R1.fastq.gz
        mv ${fastqName}_val_2.fq.gz ${fastqName}_val_R2.fastq.gz
        """
}

process trimmomatic() {
    tag "processing ${fastqName}"
    label 'trimatic'
    label 'fastqcMem'
    input:
        tuple \
            val(fastqName), \
            path(fastqReads), \
            path(adapter)
    output:
        publishDir path: "${params.outputDir}/trimmedReads/", mode: 'copy'
        path "${fastqName}_trimmed_R{1,2}_P.fq.gz"
    script:
        (readOne, readTwo) = fastqReads
        """
        trimmomatic PE \
            -threads ${task.cpus} \
            -phred33 \
            ${readOne} \
            ${readTwo} \
            ${fastqName}_trimmed_R1_P.fq.gz \
            ${fastqName}_trimmed_R1_U.fq.gz \
            ${fastqName}_trimmed_R2_P.fq.gz \
            ${fastqName}_trimmed_R2_U.fq.gz \
            ILLUMINACLIP:${adapter}:2:30:10 \
            SLIDINGWINDOW:4:15 \
            LEADING:20 \
            TRAILING:20 \
            CROP:${params.crop} \
            HEADCROP:${params.headcrop} \
            MINLEN:${params.min_length}
        """
}
