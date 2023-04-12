#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getGvcfFiles;
    getPedFile;
    joinCallVariants;
    glnexusJointCaller;
    convertBcfToVcf;
} from "${projectDir}/modules/variantCallingPipeline.nf"

workflow {
    println "\nGENOTYPE GVCF\n"
    gvcfList = getGvcfFiles().toList()

    if(params.jointCaller == "gatk") {
        combinedGvcf = combineGvcfs(gvcfList)
        ped = getPedFile(combinedGvcf)
        combinedGvcf
            .combine(ped)
            .set { join_call_input }
        vcf = joinCallVariants(join_call_input)
    } else if(params.jointCaller == 'glnexus') {
        bcf = glnexusJointCaller(gvcfList).view()
        vcf = convertBcfToVcf(bcf).view()
    } else {
        error: println "\nERROR: You must select a joint variant caller! Options are 'gatk' and 'glnexus'\n"
    }

//    ped = getPedFile(combinedGvcf)
//    combinedGvcf
//        .combine(ped)
//        .set { join_call_input }
//    vcf = joinCallVariants(join_call_input) 
}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
