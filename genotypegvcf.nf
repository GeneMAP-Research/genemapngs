#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getGvcfFiles;
    getPedFile;
    combineGvcfs;
    getVcfGenomicIntervals;
    getGenomicIntervalList;
    getGenomicsdbWorkspaces;
    genotypeGvcfs;
    createGenomicsDbPerInterval;
    callVariantsFromGenomicsDB;
    collectIntervalsPerChromosome;
    concatPerIntervalVcfs;
    glnexusJointCaller;
    convertBcfToVcf;
} from "${projectDir}/modules/variantCallingPipeline.nf"

workflow {
    println "\nGENOTYPE GVCF\n"

    gvcfList = getGvcfFiles().toList()

    if(params.joint_caller == "gatk") {

        if(params.interval == "NULL") {
            genomicInterval = getVcfGenomicIntervals(gvcfList).flatten()
        }
        else {
            genomicInterval = getGenomicIntervalList().flatten().view()

        }

        //genomicsDB = createGenomicsDbPerChromosome(per_chrom_genomicsDB_input)
        //genomicsDB = createGenomicsDbPerInterval(genomicInterval, gvcfList)
        workspace = getGenomicsdbWorkspaces().map { wrkspc -> tuple(wrkspc.simpleName, wrkspc) }
        genomicInterval
            .map { interval -> tuple(interval.simpleName + "_${params.output_prefix}-workspace", interval) }
            .join(workspace)
            .map {workspaceName, interval, workspace -> tuple(workspaceName, interval, workspace)}
            .set { workspace_interval }
        vcfs = callVariantsFromGenomicsDB(workspace_interval).collect()
        vcfs_per_chrom_list = collectIntervalsPerChromosome(vcfs).flatten()
        concatPerIntervalVcfs(vcfs_per_chrom_list).view()

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
