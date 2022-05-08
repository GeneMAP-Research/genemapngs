#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getBamFileSet;
    callVariants;
    //indexVcf;
    combineGvcfs;
    getPedFile;
    joinCallVariants
    deepVariantCaller;
    glnexusJointCaller;
    convertBcfToVcf;
} from "${projectDir}/modules/variantCallingPipeline.nf"

include {
    indexBam;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nVariant calling begins here\n"

    bamFileSet = getBamFileSet()

    if(params.singleCaller == "gatk") {
        gvcf = callVariants(bamFileSet)
    } else if(params.singleCaller == 'deepvariant') {
        gvcf = deepVariantCaller(bamFileSet).view()
    } else { 
        error: println "\nERROR: You must select a single sample variant caller! Options are 'gatk' and 'deepvariant'\n" 
    }

    gvcfList = gvcf.collect()

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
}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
