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
        }

    } else {
        error: println "\nERROR: Joint caller must be 'gatk'\n"
    }

}

workflow.onComplete { println "\nDone! Check results in ${params.output_dir}\n" }
