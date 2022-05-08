/*----- SPLIT FASTQs INTO CHUNKS AND CONCAT SAM CHUNKS: not implemented -----
*  def getFastqChunks() {
*      return channel
*          .fromFilePairs( params.fastqDir + "*_R{1,2}*[fq,fastq]*", flat:true )
*          .splitFastq(
*              by: 1000000,
*              pe: true,
*              file: true,
*              compress: true )
*  }
*  
*  def concatenateAlignedChunks(alignedChunks) {
*      return alignedChunks
*          .collectFile(
*              name: "concatenated.lgen",
*              sort: true)
*  }
*/

def getInputBams() {
    return channel.fromFilePairs( params.inputDir + "*.bam", size: 1 )
                  .ifEmpty { error "\nERROR: Could not locate BAM files!\nPlease make sure you have specified the correct file type and that they exist in the input directory you specified" }
                  .map { bamName, bamFile -> tuple(bamName, bamFile.first()) }
}

def getInputFastqs() {
    return channel.fromFilePairs( params.inputDir + "*_R{1,2}*.[fq,fastq]*", size: 2 )
                  .ifEmpty { error "\nERROR: Some fastq files could not be found!\n" }
                  .map { fqBase, fastq -> tuple(fqBase, fastq) }
}

process sortBamByName() {
    tag "processing ${bamName}"
    label 'samtools'
    label 'bamSorter'
    input:
        tuple \
            val(bamName), \
            path(bamFile) 
    output:
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
    label 'bamSorter'
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

process alignReadsToReference() {
    tag "processing ${fastqName}"
    label 'bwa_bgzip'
    label 'readAligner'
    input:
        tuple \
            val(fastqName), \
            path(reads)
    output:
        tuple \
            val(fastqName), \
            path("${fastqName}.sam.gz")
    script:
        ( readOne, readTwo ) = reads
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

process buildRefHashTable() {
    tag "BULIDING HASH TABLE"
    label 'dragmap'
    label 'longRun'
    input:
        path fastaRef
    output:
        publishDir  path: "${params.referenceDir}", mode: 'copy'
        path "*"
    script:
        """
        dragen-os \
            --build-hash-table true \
            --ht-reference ${fastaRef}  \
            --output-file-prefix=${fastaRef.baseName}
        """
}

process dragenAligner() {
    tag "processing ${fastqName}"
    label 'dragmap'
    label 'readAligner'
    input:
        tuple \
            val(fastqName), \
            path(reads)
    output:
        tuple \
            val(fastqName), \
            path("${fastqName}.sam.gz")
    script:
        ( readOne, readTwo ) = reads
        """
        dragen-os \
            -r ${params.fastaRef} \
            -1 ${readOne} \
            -2 ${readTwo} | \
            bgzip -f -@ ${task.cpus} -c > ${fastqName}.sam.gz
        """
}

//--- TMAP ALIGNER FOR ION TORRENT ---

/*  process tmapAligner() {
*      tag "processing ${fastqName}"
*      label 'tmap'
*      label 'readAligner'
*      input:
*          tuple \
*              val(fastqName), \
*              path(reads)
*      output:
*          tuple \
*              val(fastqName), \
*              path("${fastqName}.sam.gz")
*      script:
*          ( readOne, readTwo ) = reads
*          """
*          tmap \
*              mapall \
*              -f ${params.fastaRef} \
*              -r ${readOne} \
*              -n ${task.cpus} \
*              -v stage1 map1 map2 map3 | \
*              bgzip -f -@ ${task.cpus} -c > "${fastqName}.sam.gz"
*          """
*  }
*/  

process convertSamToBam() {
    tag "processing ${samName}"
    label 'samtools'
    label 'samConverter'
    input:
        tuple \
            val(samName), \
            path(samFile) 
    output:
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

process sortBam() {
    tag "processing ${bamName}"
    label 'samtools'
    label 'bamSorter'
    input:
        tuple \
            val(bamName), \
            path(bamFile)
    output:
        tuple \
            val(bamName), \
            path("${bamName}.sorted.bam")
    script:
        """
        samtools \
            sort \
            --reference ${params.fastaRef} \
            -O BAM \
            --threads ${task.cpus} \
            -o "${bamName}.sorted.bam" \
            ${bamFile}
        """
}

process indexBam() {
    tag "processing ${bamName}"
    label 'samtools'
    label 'bamIndexer'
    input:
        tuple \
            val(bamName), \
            path(bamFile)
    output:
        tuple \
            val(bamName), \
            path("${bamFile}"), \
            path("${bamFile}.bai")
    script:
        bam = bamFile[0]
        """
        samtools \
            index \
            -@ ${task.cpus} \
            ${bamFile}
        """
}
 
process markDuplicates() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'duplicateMarker'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        tuple \
            val(bamName), \
            path("${bamName}.dupsMarked.bam")
    script:
        """
        gatk \
            --java-options "-XX:ParallelGCThreads=${task.cpus} -Xmx${task.memory.toGiga()}g" \
            MarkDuplicates \
            -I ${bamFile} \
            -O "${bamName}.dupsMarked.bam" \
            -M "${bamName}.bamMetrics.txt" 
        """
}
 
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
        tuple \
            val(bamName), \
            path("${bamName}.recal-table.txt")
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xmx${task.memory.toGiga()}g" \
            BaseRecalibrator \
            -I ${bamFile} \
            -O "${bamName}.recal-table.txt" \
            -R ${params.fastaRef} \
            --known-sites ${params.dbsnp} \
            --known-sites ${params.indelsRef} \
            --known-sites ${params.hapmap} \
            --known-sites ${params.omniRef} \
            --known-sites ${params.snpRef}
        """
}

process applyBaseQualityRecalibrator() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'applyBqsr'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex), \
            path(recalTable) 
    output:
        publishDir path: "${params.outputDir}/bam/", mode: 'copy'
        tuple \
            val(bamName), \
            path("${bamName}.bqsr.bam"), \
            path("${bamName}.bqsr.bam.bai")
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xmx${task.memory.toGiga()}g" \
            ApplyBQSR \
            -I ${bamFile} \
            -O "${bamName}.bqsr.bam" \
            -bqsr "${recalTable}"

        mv ${bamName}.bqsr.bai ${bamName}.bqsr.bam.bai
        """
}

/*------------- GATK SPARK PIPELINES -----------*/
process markDuplicatesSpark() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'duplicateMarkerSpark'
    input:
        tuple \
            val(bamName), \
            path(bamFile)
            //path(bamIndex)
    output:
        tuple \
            val(bamName), \
            path("${bamName}.dupsMarked.bam"), \
            path("${bamName}.dupsMarked.bam.bai")
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xmx${task.memory.toGiga()}g" \
            MarkDuplicatesSpark \
            -I ${bamFile} \
            -O "${bamName}.dupsMarked.bam" \
            -M "${bamName}.bamMetrics.txt" \
            -OBI true \
            -- \
            --conf 'spark.local.dir=${workDir}/temp'
            --conf 'spark.executor.cores=${task.cpus}'
            --spark-master local[${task.cpus}]
        """
}

process fixBamTags() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'fixBam'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
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

process recalibrateBaseQualityScoresSpark() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'baseRecalibratorSpark'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        publishDir path: "${params.outputDir}/bam/", mode: 'copy'
        tuple \
            val(bamName), \
            path("${bamName}.bqsr.bam"), \
            path("${bamName}.bqsr.bam.bai")
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xmx${task.memory.toGiga()}g" \
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

