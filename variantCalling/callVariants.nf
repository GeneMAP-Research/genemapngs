#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getBamFile;
    getBamIndex;
    callVariants;
    //indexVcf;
    combineGvcfs;
    getPedFile;
    joinCallVariants
} from "${projectDir}/modules/variantCallingPipeline.nf"

include {
    indexBam;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nVariant calling begins here\n"
    bam = getBamFile()
    bamFileSet = indexBam(bam)
    gvcf = callVariants(bamFileSet)
    gvcfList = gvcf.collect()
    combinedGvcf = combineGvcfs(gvcfList)
    ped = getPedFile(combinedGvcf)
    combinedGvcf
        .combine(ped)
        .set { join_call_input }
    vcf = joinCallVariants(join_call_input) 
}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
