#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    //getGvcfFile;
    getPedFile;
    joinCallVariants
} from "${projectDir}/modules/variantCallingPipeline.nf"

workflow {
    println "\nGENOTYPE GVCF\n"
    combinedGvcf = getGvcfFile().toList().view()
    ped = getPedFile(combinedGvcf)
    combinedGvcf
        .combine(ped)
        .set { join_call_input }
    vcf = joinCallVariants(join_call_input) 
}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }

params.gvcf_dir = '/scratch/eshkev001/projects/scd-wgs/trypanogen/work/96/d957066cff06593e34aedd1a0fb86a/' 
params.gvcf_base = 'trypanogen'

def getGvcfFile() {
    return channel.fromPath(params.gvcf_dir + params.gvcf_base + '.*')
                  .ifEmpty { error "\nERROR: Could not locate a file! \n" }
}
