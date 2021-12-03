#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getBamFiles;
    sortBamByName;
    convertBamToFastq;
} from "${projectDir}/modules/bamToFastqPipeline.nf"

workflow {
    println "\nBAM2FASTQ begins here\n"

    bam = getBamFiles()
    bamSortedByName = sortBamByName(bam)
    fastq = convertBamToFastq(bamSortedByName)
}

