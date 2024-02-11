#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getGvcfFiles;
    getPedFile;
    combineGvcfs;
    jointCallVariants;
    createGenomicsDbPerChromosome;
    callVariantsFromGenomicsDB;
    glnexusJointCaller;
    convertBcfToVcf;
} from "${projectDir}/modules/variantCallingPipeline.nf"

workflow {
    println "\nGENOTYPE GVCF\n"

    gvcfList = getGvcfFiles().toList()

    if(params.jointCaller == "gatk") {

        channel.from(1..22,'X','Y','M')
               .collect()
               .flatten()
               .map { chr -> "chr${chr}" }
               .combine(gvcfList.toList())
               .set { per_chrom_genomicsDB_input }
 
        genomicsDB = createGenomicsDbPerChromosome(per_chrom_genomicsDB_input).view()
        callVariantsFromGenomicsDB(genomicsDB).view()

/*
        combinedGvcf = combineGvcfs(gvcfList)
        ped = getPedFile(combinedGvcf)
        combinedGvcf
            .combine(ped)
            .set { join_call_input }
        vcf = jointCallVariants(join_call_input)
*/

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
