#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getGvcfFiles;
    getGenomicsdbWorkspaces;
    getPedFile;
    combineGvcfs;
    getVcfGenomicIntervals;
    getGenomicIntervalList;
    genotypeGvcfs;
    createGenomicsDbPerInterval;
    updateGenomicsDbPerInterval;
    callVariantsFromGenomicsDB;
    glnexusJointCaller;
    convertBcfToVcf;
} from "${projectDir}/modules/variantCallingPipeline.nf"

workflow {
    println "\nGENOTYPE GVCF\n"

    gvcfList = getGvcfFiles().toList()

    if(params.joint_caller == "gatk") {

        //channel.from(1..22,'X','Y','M')
        //       .collect()
        //       .flatten()
        //       .map { chr -> "chr${chr}" }
        //       .combine(gvcfList.toList())
        //       .set { per_chrom_genomicsDB_input }
        //
        //genomicsDB = createGenomicsDbPerChromosome(per_chrom_genomicsDB_input)


        if(params.interval == "NULL") {
            genomicInterval = getVcfGenomicIntervals(gvcfList).flatten()
        }
        else {
            genomicInterval = getGenomicIntervalList().flatten()

        }

        if(params.update == true) {
            workspace = getGenomicsdbWorkspaces().map { wrkspc -> tuple(wrkspc.simpleName, wrkspc) }
            genomicInterval
                .map { interval -> tuple(interval.simpleName + "_${params.output_prefix}-workspace", interval) }
                .join(workspace)
                .map {workspaceName, interval, workspace -> tuple(workspaceName, interval, workspace)}
                .set { workspace_interval }
            updateGenomicsDbPerInterval(workspace_interval, gvcfList)
        } else {
            genomicsDB = createGenomicsDbPerInterval(genomicInterval, gvcfList)
            //callVariantsFromGenomicsDB(genomicsDB).view()
        }

/*
        combinedGvcf = combineGvcfs(gvcfList)
        ped = getPedFile(combinedGvcf)
        combinedGvcf
            .combine(ped)
            .set { join_call_input }
        vcf = genotypeGvcfs(join_call_input)
*/

    } else if(params.joint_caller == 'glnexus') {
        bcf = glnexusJointCaller(gvcfList).view()
        vcf = convertBcfToVcf(bcf).view()
    } else {
        error: println "\nERROR: You must select a joint variant caller! Options are 'gatk' and 'glnexus'\n"
    }


//    ped = getPedFile(combinedGvcf)
//    combinedGvcf
//        .combine(ped)
//        .set { join_call_input }
//    vcf = genotypeGvcfs(join_call_input) 
}

workflow.onComplete { println "\nDone! Check results in ${params.output_dir}\n" }
