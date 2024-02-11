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
    return channel.fromFilePairs( params.inputDir + "/*.bam", size: 1 )
                  .ifEmpty { error "\nERROR: Could not locate BAM files!\nPlease make sure you have specified the correct file type and that they exist in the input directory you specified" }
                  .map { bamName, bamFile -> tuple(bamName, bamFile.first()) }
}

def getInputCram() {
    return channel.fromFilePairs( params.inputDir + "/*.cram", size: 1 )
                  .ifEmpty { error "\nERROR: Could not locate BAM files!\nPlease make sure you have specified the correct file type and that they exist in the input directory you specified" }
                  .map { bamName, bamFile -> tuple(bamName, bamFile.first()) }
}

// "*_*{1,2}*.[fq,fastq]*"
def getInputFastqs() {
    println "\nWARN: Only file pairs with the extension '*_R{1,2}*.fastq.gz' will be read!\nIf there are other files without this extension, please rename them and run the workflow again.\n"
    return channel.fromFilePairs( params.inputDir + './*R{1,2}*.fastq.gz', size: 2)
                  .ifEmpty { error "\nERROR: Some fastq files could not be found!\n" }
                  .map { fqBase, fastq -> tuple(fqBase, fastq) }
}

def getSEInputFastqs() {
    return channel.fromFilePairs( params.inputDir + "/*.[fq,fastq]*", size: 1 )
                  .ifEmpty { error "\nERROR: Some fastq files could not be found!\n" }
                  .map { fqBase, fastq -> tuple(fqBase, fastq.first()) }
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

process bwaAligner() {
    tag "processing ${fastqName}"
    label 'bwa_bgzip'
    label 'readAligner'
    cache 'lenient'
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
            -r ${params.ref_dir} \
            -1 ${readOne} \
            -2 ${readTwo} \
            --RGID ${fastqName} \
            --RGSM ${fastqName} | \
            bgzip -f -@ ${task.cpus} -c > ${fastqName}.sam.gz
        """
}

//--- TMAP ALIGNER FOR ION TORRENT ---

process tmapAligner() {
    tag "processing ${fastqName}"
    label 'tmap'
    label 'tmap_mem'
    input:
        tuple \
            val(fastqName), \
            path(fastqFile)
    output:
        publishDir path: "${params.output_dir}/bam/"
        tuple \
            val(fastqName), \
            path("${fastqName}.bam")
    script:
        """
        tmap \
            mapall \
            -f ${params.fastaRef} \
            -o 1 \
            -R ID:${fastqName} \
            -R SM:${fastqName} \
            -R PG:TMAP \
            -R PL:IONTORRENT \
            --end-repair 1 \
            -J 25 \
            -r ${fastqFile} \
            -n ${task.cpus} \
            -v stage1 map1 map2 map3 \
            > "${fastqName}.bam"
        """
}

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

process fixAlignmentMate() {
    tag "processing ${bamName}"
    label 'samtools'
    label 'bamSorter'
    //publishDir \
    //    path: "${params.outputDir}/cram/", \
    //    mode: 'copy'
    input:
        tuple \
            val(bamName), \
            path(bamFile)
    output:
        tuple \
            val(bamName), \
            path("${bamName}.cram")
    script:
        """
        samtools \
            fixmate \
            -c \
            -m \
            -O CRAM \
            --reference ${params.fastaRef} \
            --threads ${task.cpus} \
            ${bamFile} \
            ${bamName}.cram
        """
}

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
            path("${samName}.cram")
    script:
        """
        samtools \
            view \
            -h \
            ${samFile} \
            -O CRAM \
            --threads ${task.cpus} \
            -o "${samName}.cram"
        """
}

process convertBamToCram() {
    tag "processing ${samBamName}"
    label 'samtools'
    label 'samConverter'
    publishDir \
        path: "${params.outputDir}/cram/", \
        mode: 'copy'    
    input:
        tuple \
            val(samBamName), \
            path(samBamFile)
    output:
        tuple \
            val(samBamName), \
            path("${samBamName}.cram"), \
            path("*.crai")
    script:
        """
        samtools \
            view \
            -h \
            -O CRAM \
            --reference ${params.fastaRef} \
            --threads ${task.cpus} \
            --write-index \
            ${samBamFile} \
            -o ${samBamName}.cram
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

process sortCram() {
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
            path("${bamName}.sorted.cram")
    script:
        """
        samtools \
            sort \
            --reference ${params.fastaRef} \
            -O CRAM \
            --threads ${task.cpus} \
            -o "${bamName}.sorted.cram" \
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
            path("${bamFile}.[bai,crai]")
    script:
        bam = bamFile[0]
        """
        samtools \
            index \
            -@ ${task.cpus} \
            ${bamFile}
        """
}

process indexAndCopyBam() {
    tag "processing ${bamName}"
    label 'samtools'
    label 'bamIndexer'
    publishDir \
        path: "${params.outputDir}/bam/", \
        mode: 'copy'
    input:
        tuple \
            val(bamName), \
            path(bamFile)
    output:
        tuple \
            val(bamName), \
            path("${bamFile}"), \
            path("*.bai")
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
    label 'samtools'
    label 'bamSorter'
    //publishDir \
    //    path: "${params.outputDir}/markedcram/", \
    //    mode: 'copy'
    input:
        tuple \
            val(bamName), \
            path(bamFile)
    output:
        tuple \
            val(bamName), \
            path("${bamName}.dupsMarked.cram")
    script:
        """
        samtools \
            markdup \
            --reference ${params.fastaRef} \
            -O CRAM \
            --threads ${task.cpus} \
            ${bamFile} \
            ${bamName}.dupsMarked.cram
        """
}
 
process markDuplicatesGatk() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'duplicateMarker'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        //publishDir path: "${params.outputDir}/markedbam/", mode: 'copy'
        tuple \
            val(bamName), \
            path("${bamName}.dupsMarked.bam"), \
            path("${bamName}.dupsMarked.bam.bai")
    script:
        """
        gatk \
            --java-options "-XX:ParallelGCThreads=${task.cpus} -Xmx${task.memory.toGiga()}g" \
            MarkDuplicates \
            -I ${bamFile} \
            --CREATE_INDEX true \
            -O "${bamName}.dupsMarked.bam" \
            -M "${bamName}.bamMetrics.txt"

        mv ${bamName}.dupsMarked.bai ${bamName}.dupsMarked.bam.bai
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
    publishDir \
        path: "${params.outputDir}/cram/", \
        mode: 'copy'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex), \
            path(recalTable) 
    output:
        tuple \
            val(bamName), \
            path("${bamName}.bqsr.cram"), \
            path("${bamName}.bqsr.cram.crai")
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xmx${task.memory.toGiga()}g" \
            ApplyBQSR \
            -I ${bamFile} \
            -O "${bamName}.bqsr.cram" \
            -bqsr "${recalTable}"

        mv ${bamName}.bqsr.crai ${bamName}.bqsr.bam.crai
        """
}

/*------------- GATK SPARK PIPELINES -----------*/
process markDuplicatesSpark() {
    beforeScript 'ulimit -c unlimited'
    cache 'lenient'
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
            --conf 'spark.local.dir=${workDir}/temp' \
            --conf 'spark.executor.cores=${task.cpus}' \
            --spark-master local[${task.cpus}]
        """
}

process fixBamTags() {
    tag "processing ${bamName}"
    beforeScript 'ulimit -c unlimited'
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
    beforeScript 'ulimit -c unlimited'
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

