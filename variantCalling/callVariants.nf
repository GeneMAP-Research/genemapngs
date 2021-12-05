#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getBamFiles;
    callVariants;
    //indexVcf;
    combineGvcfs;
    getPedFile;
    joinCallVariants
} from "${projectDir}/modules/variantCallingPipeline.nf"

workflow {
    println "\nVariant calling begins here\n"
    bam = getBamFiles()
    gvcf = callVariants(bam)
    gvcfList = gvcf.collect()
    combinedGvcf = combineGvcfs(gvcfList)
    ped = getPedFile(combinedGvcf)
    combinedGvcf
        .combine(ped)
        .set { join_call_input }
    vcf = joinCallVariants(join_call_input)
}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
