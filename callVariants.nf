#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getBamFileSet;
    getCramFileSet;
    getAlignmentFileSet;
    getAlignmentGenomicIntervals;
    haplotypeCaller;
    haplotypeCallerWithIntervals;
    concatGvcfs;
    haplotypeCallerSpark;
    //combineGvcfsPerInterval;
    getPedFile;
    deepVariantCaller;
    glnexusJointCaller;
    dysguCallSvs;
    indexVcf;
    dysguMergeVcfs;
    getBamChuncks;
    mantaCallSvs;
    mergeMantaCandidateSvCalls;
    mergeMantaCandidateSmallIndelCalls;
    mergeMantaDiploidSvCalls;
    convertBcfToVcf;
    createGenomicsDb;
    createGenomicsDbPerInterval;
    callVariantsFromGenomicsDB;
    getGvcfFiles;
    getGenomicsdbWorkspaces;
    combineGvcfs;
    getVcfGenomicIntervals;
    getGenomicIntervalList;
    genotypeGvcfs;
    updateGenomicsDbPerInterval;
} from "${projectDir}/modules/variantCallingPipeline.nf"

include {
    indexAlignment;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nVariant calling begins here\n"

    bamFileSet = getAlignmentFileSet()

    /***************************************************************
    * ONE-STEP VARIANT CALLING FROM ALIGNMENT FILES: DYSGY & MANTA *
    *          SINGLE SAMPLE STRUCTURAL VARIANT CALLING            *
    ***************************************************************/
     
    if(params.single_caller.toUpperCase() == 'DYSGU') {
        dysguCallSvs(bamFileSet)
            .set { vcf }
        indexVcf(vcf)
            .collect()
            .set { vcfs }
        dysguMergeVcfs(vcfs)
    }
    else if(params.single_caller.toUpperCase() == 'MANTA') {
        bamFileSet
            .map { bamName, bamFile, bamIndex -> tuple(bamFile, bamIndex) }
            .collect()
            .set { manta_input }
        getBamChuncks(manta_input)
            .flatten().view()
            .set { bamlist }
        mantaCallSvs(bamlist)
            .collect().view()
            .set { manta_calls }
        mergeMantaDiploidSvCalls(manta_calls).view()
        mergeMantaCandidateSmallIndelCalls(manta_calls).view()
        mergeMantaCandidateSvCalls(manta_calls).view()
    }
    else {

    /*************************************************************************
    * TWO-STEP SMALL VARIANTS CALLING VIA GVCFs: GATK, DEEPVARIANT & GLNEXUS *
    *************************************************************************/

        if(params.single_caller.toUpperCase() == 'DEEPVARIANT') {
            gvcf = deepVariantCaller(bamFileSet)
        }
        else { // SINGLE CALLER DEFAULTS TO GATK //

            //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
            // haplotype caller per interval seems too cumbersome when  //
            // dealing with large samples sizes since each sample would //
            // have to be split into hundreds - thousands of intervals  //
            //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
 
            gvcf = haplotypeCaller(bamFileSet)
        }
    
        gvcfList = gvcf.collect().view()

        if(params.joint_caller.toUpperCase() == 'GLNEXUS')  {
            bcf = glnexusJointCaller(gvcfList)
            vcf = convertBcfToVcf(bcf)
        } 
        else {

            //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
            // If mode is not svarcall and not specifically jvarcall,   //
            // single sample variant calling leads directly to joint    //
            // calling                                                  //
            //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//


            if(params.mode != 'svarcall') {

                /**************************************************
                *          JOINT CALLER DEFAULTS TO GATK          *
                *    Using 'GenomicsDBImport' for efficiency      *
                **************************************************/

                if(params.interval == "NULL") {
                    genomicInterval = getVcfGenomicIntervals(gvcfList).flatten()
                }
                else {
                    genomicInterval = getGenomicIntervalList().flatten()
                }

                if(params.mode == 'jvarcall') {

                // genomicsdb workspaces must already exist //

                    workspace = getGenomicsdbWorkspaces().map { wrkspc -> tuple(wrkspc.simpleName, wrkspc) }
                    genomicInterval
                        .map { interval -> tuple(interval.simpleName + "_${params.output_prefix}-workspace", interval) }
                        .join(workspace)
                        .map {workspaceName, interval, workspace -> tuple(workspaceName, interval, workspace)}
                        .set { workspace_interval }
                    vcfs = callVariantsFromGenomicsDB(workspace_interval).collect()
                    vcfs_per_chrom_list = collectIntervalsPerChromosome(vcfs).flatten()
                    concatPerIntervalVcfs(vcfs_per_chrom_list).view()
                }
                else { // DEFAULT TO SINGLE AND JOINT CALLING

                    genomicsDB = createGenomicsDbPerInterval(genomicInterval, gvcfList)

                    workspace = genomicsDB.map { wrkspc -> tuple(wrkspc.simpleName, wrkspc) }
                    genomicInterval
                        .map { interval -> tuple(interval.simpleName + "_${params.output_prefix}-workspace", interval) }
                        .join(workspace)
                        .map {workspaceName, interval, workspace -> tuple(workspaceName, interval, workspace)}
                        .set { workspace_interval }
                    vcfs = callVariantsFromGenomicsDB(workspace_interval).collect()
                    vcfs_per_chrom_list = collectIntervalsPerChromosome(vcfs).flatten()
                    concatPerIntervalVcfs(vcfs_per_chrom_list).view()
                }

            }

            }
        }
    }
}

workflow.onComplete { println "\nDone! Check results in ${params.output_dir}\n" }
