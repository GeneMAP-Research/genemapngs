/*----- SPLIT FASTQs INTO CHUNKS AND CONCAT SAM CHUNKS: not implemented -----
def getFastqChunks() {
    return channel
        .fromFilePairs( params.fastqDir + "*_R{1,2}*[fq,fastq]*", flat:true )
        .splitFastq(
            by: 1000000,
            pe: true,
            file: true,
            compress: true )
}

def concatenateAlignedChunks(alignedChunks) {
    return alignedChunks
        .collectFile(
            name: "concatenated.lgen",
            sort: true)
}
----------------------------------------------------------------------------*/

/*
def getInputs() {
    if(params.inputMode == "fastq")
        getFastq().vie()
    else if(params.inputMode == "bam")
        getBamFiles().view()
    
}

def getBamFiles() {
    return channel.fromFilePairs( params.bamDir + "*.bam", size: 1 )
                  .map { bamName, bamFile -> tuple(bamName, bamFile.first()) }
}
*/

def getFastq() {
    return channel.fromFilePairs( params.fastqDir + "*_R{1,2}*[fq,fastq]*", size: 2 )
                  .ifEmpty { error "\nERROR: Some fastq files could not be found!\n" }
                  .map { fqBase, fastq -> tuple(fqBase, fastq) }
}

process alignReadsToReference() {
    tag "processing ${fastqName}"
    label 'bwa_bgzip'
    label 'readAligner'
    input:
        tuple \
            val(fastqName), \
            path(reads)
    output:
        publishDir path: "${params.outputDir}/fastq/"
        tuple \
            val(fastqName), \
            path("${fastqName}.sam.gz")
    script:
        (readOne, readTwo) = reads
        """
        bwa \
            mem \
            -t ${task.cpus} \
            ${params.fastaRef} \
            ${readOne} \
            ${readTwo} \
            -R \"@RG\\tID:${fastqName}\\tSM:${fastqName}\\tPL:ILLUMINA\" | \
            bgzip -f -@ ${task.cpus} -c > "${fastqName}.sam.gz"
        """
}

process convertSamToBam() {
    tag "processing ${samName}"
    label 'samtools'
    label 'mediumMemory'
    input:
        tuple \
            val(samName), \
            path(samFile) 
    output:
        publishDir path: "${params.outputDir}/bam/"
        tuple \
            val(samName), \
            path("${samName}.bam")
    script:
        """
        samtools \
            view \
            -h \
            ${samFile} \
            -O BAM \
            --threads ${task.cpus} \
            -o "${samName}.bam"
        """
}

/*
*  process sortBam() {
*      tag "processing ${bamName}"
*      label 'samtools'
*      label 'bamSorter'
*      input:
*          tuple \
*              val(bamName), \
*              path(bamFile)
*      output:
*          publishDir path: "${params.outputDir}/bam/"
*          tuple \
*              val(bamName), \
*              path("${bamName}.sorted.bam")
*      script:
*          """
*          samtools \
*              sort \
*              --reference ${params.fastaRef} \
*              -O BAM \
*              --threads ${task.cpus} \
*              -o "${bamName}.sorted.bam" \
*              ${bamFile}
*          """
*  }
*
*    
*  process indexBam() {
*      tag "processing ${bamName}"
*      label 'samtools'
*      label 'bamIndexer'
*      input:
*          tuple \
*              val(bamName), \
*              path(bamFile)
*      output:
*          publishDir path: "${params.outputDir}/bam/"
*          tuple \
*              val(bamName), \
*              path("${bamFile}"), \
*              path("${bamFile}.bai")
*      script:
*          """
*          samtools \
*              index \
*              -@ ${task.cpus} \
*              ${bamFile}
*          """
*  }
*  
*  
*  process markDuplicates() {
*      tag "processing ${bamName}"
*      label 'gatk'
*      label 'duplicateMarker'
*      input:
*          tuple \
*              val(bamName), \
*              path(bamFile), \
*              path(bamIndex)
*      output:
*          publishDir path: "${params.outputDir}/bam/"
*          tuple \
*              val(bamName), \
*              path("${bamName}.dupsMarked.bam")
*      script:
*          """
*          gatk \
*              MarkDuplicates \
*              -I ${bamFile} \
*              -O "${bamName}.dupsMarked.bam" \
*              -M "${bamName}.bamMetrics.txt" \
*              --REMOVE_DUPLICATES
*          """
*  }
*/  

/*
*
*  Marked duplicate reads will be indexed 
*  by indexBam called as indexMarkedBam
*  in the workflow script
*  
*  The output will be a tuple suitable for
*  the next process:
*  - bamName
*  - bamFile
*  - bamIndex
*
******************************************/

/*
*  process recalibrateBaseQualityScores() {
*      tag "processing ${bamName}"
*      label 'gatk'
*      label 'baseRecalibrator'
*      input:
*          tuple \
*              val(bamName), \
*              path(bamFile), \
*              path(bamIndex)
*      output:
*          publishDir path: "${params.outputDir}/bam/"
*          tuple \
*              val(bamName), \
*              path("${bamName}.recal-table.txt")
*      script:
*          """
*          gatk \
*              BaseRecalibrator \
*              -I ${bamFile} \
*              -O "${bamName}.recal-table.txt" \
*              -R ${params.fastaRef} \
*              --known-sites ${params.dbsnp} \
*              --known-sites ${params.indelsRef} \
*              --known-sites ${params.hapmap} \
*              --known-sites ${params.omniRef} \
*              --known-sites ${params.snpRef}
*          """
*  }
*  
*  process applyBaseQualityRecalibrator() {
*      tag "processing ${bamName}"
*      label 'gatk'
*      label 'applyBqsr'
*      input:
*          tuple \
*              val(bamName), \
*              path(bamFile), \
*              path(bamIndex), \
*              path(recalTable) 
*      output:
*          publishDir path: "${params.outputDir}/bam/", mode: 'copy'
*          tuple \
*              val(bamName), \
*              path("${bamName}.bqsr.bam"), \
*              path("${bamName}.bqsr.bai")
*      script:
*          """
*          gatk \
*              ApplyBQSR \
*              -I ${bamFile} \
*              -O "${bamName}.bqsr.bam" \
*              -bqsr "${recalTable}"
*          """
*  }
*/

/*------------- GATK SPARK PIPELINES -----------*/
process buildBamIndex() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'bamIndexer'
    input:
        tuple \
            val(bamName), \
            path(bamFile)
    output:
        publishDir path: "${params.outputDir}/bam/"
        tuple \
            val(bamName), \
            path("${bamFile}"), \
            path("${bamName}.bai")
    script:
        """
        gatk \
            BuildBamIndex \
            -I ${bamFile} \
            -R ${params.fastaRef} \
            -O ${bamName}.bai
        """
}

process markDuplicates() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'duplicateMarker'
    input:
        tuple \
            val(bamName), \
            path(bamFile)
            //path(bamIndex)
    output:
        publishDir path: "${params.outputDir}/bam/"
        tuple \
            val(bamName), \
            path("${bamName}.dupsMarked.bam"), \
            path("${bamName}.dupsMarked.bam.bai")
    script:
        """
        gatk \
            MarkDuplicatesSpark \
            -I ${bamFile} \
            -O "${bamName}.dupsMarked.bam" \
            -M "${bamName}.bamMetrics.txt" \
            -OBI true \
            -- \
            --spark-master local[${task.cpus}]
        """
}

process fixBamTags() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'duplicateMarker'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        publishDir path: "${params.outputDir}/bam/"
        tuple \
            val(bamName), \
            path("${bamName}.fixed.bam"), \
            path("${bamName}.fixed.bai")
    script:
        """
        gatk \
            SetNmMdAndUqTags \
            -I ${bamFile} \
            -O ${bamName}.fixed.bam \
            -R ${params.fastaRef} \
            --CREATE_INDEX true
        """
}

process recalibrateBaseQualityScores() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'baseRecalibrator'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        publishDir path: "${params.outputDir}/bam/"
        tuple \
            val(bamName), \
            path("${bamName}.bqsr.bam"), \
            path("${bamName}.bqsr.bam.bai")
    script:
        """
        gatk \
            BQSRPipelineSpark \
            -I ${bamFile} \
            -O "${bamName}.bqsr.bam" \
            -R ${params.fastaRef} \
            --known-sites ${params.dbsnp} \
            --known-sites ${params.indelsRef} \
            --known-sites ${params.hapmap} \
            --known-sites ${params.omniRef} \
            --known-sites ${params.snpRef} \
            -- \
            --spark-master local[${task.cpus}]
        """
}

