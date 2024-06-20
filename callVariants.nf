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
    callVariantsFromExistingGenomicsDB;
    collectIntervalsPerChromosome;
    concatPerChromIntervalVcfs;
    getGvcfFiles;
    getGenomicsdbWorkspaces;
    combineGvcfs;
    getVcfGenomicIntervals;
    getGenomicIntervalList;
    getGenomicInterval;
    genotypeGvcfs;
    updateGenomicsDbPerInterval;
} from "${projectDir}/modules/variantCallingPipeline.nf"

include {
    indexAlignment;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nVariant calling begins here\n"

    if(params.mode == 'jvarcall') {

        //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        // YOU'RE EITHER STARING A NEW IMPORT OR UPDATING EXISTING GENOMICSDBs
        // OR CALLING VARIANTS FROM EXISTING GENOMICSDBs
        //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        if(params.imprt == true) { 

            //-=-=-=-=-=-=
            // NEW IMPORT
            //-=-=-=-=-=-=
 
           gvcfList = getGvcfFiles().toList()

           genomicInterval = getGenomicInterval(gvcfList)
           genomicsDB = createGenomicsDbPerInterval(genomicInterval, gvcfList)

        }
        else if(params.update == true) {

            //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            // UPDATING EXISTING GENOMICSDBs:
            // genomicsdb workspaces must already exist and gvcfs must be provided with interval list 
            //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            gvcfList = getGvcfFiles().toList()        
            genomicInterval = getGenomicInterval()
            workspace = getGenomicsdbWorkspaces().map { wrkspc -> tuple(wrkspc.simpleName, wrkspc) }
            genomicInterval
                .map { interval -> tuple(interval.simpleName + "_${params.output_prefix}-workspace", interval) }
                .join(workspace)
                .map {workspaceName, interval, workspace -> tuple(workspaceName, interval, workspace)}
                .set { workspace_interval }
            updateGenomicsDbPerInterval(workspace_interval, gvcfList)
        }
        else { 

            //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            // DEFAULT TO CALL VARIANTS FROM EXISTING GENOMICSDBs:
            // genomicsdb workspaces must already exist
            //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    
            workspace = getGenomicsdbWorkspaces().map { wrkspc -> tuple(wrkspc.simpleName, wrkspc) }
            vcfs = callVariantsFromExistingGenomicsDB(workspace).view().collect()
            vcfs_per_chrom_list = collectIntervalsPerChromosome(vcfs).flatten()
            concatPerChromIntervalVcfs(vcfs_per_chrom_list).view()

        }
    }
    else {

        bamFileSet = getAlignmentFileSet()

        //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        // ONE-STEP VARIANT CALLING FROM ALIGNMENT FILES: DYSGY & MANTA 
        // SINGLE SAMPLE STRUCTURAL VARIANT CALLING            
        //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         
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

        //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        // TWO-STEP SMALL VARIANT CALLING VIA GVCFs: GATK, DEEPVARIANT & GLNEXUS
        //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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


            if(params.mode == 'varcall') {
                if(params.joint_caller.toUpperCase() == 'GLNEXUS')  {
                    bcf = glnexusJointCaller(gvcfList)
                    vcf = convertBcfToVcf(bcf)
                }
                else {

                    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                    // JOINT CALLER DEFAULTS TO GATK
                    // Using 'GenomicsDBImport' for efficiency
                    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                    genomicInterval = getGenomicInterval()

                    genomicsDB = createGenomicsDbPerInterval(genomicInterval, gvcfList)

                    workspace = genomicsDB.map { wrkspc -> tuple(wrkspc.simpleName, wrkspc) }
                    genomicInterval
                        .map { interval -> tuple(interval.simpleName + "_${params.output_prefix}-workspace", interval) }
                        .join(workspace)
                        .map {workspaceName, interval, workspace -> tuple(workspaceName, interval, workspace)}
                        .set { workspace_interval }
                    vcfs = callVariantsFromGenomicsDB(workspace_interval).collect()
                    vcfs_per_chrom_list = collectIntervalsPerChromosome(vcfs).flatten()
                    concatPerChromIntervalVcfs(vcfs_per_chrom_list).view()
                }
            }

        }

    }

}

workflow.onComplete { println "\nDone! Check results in ${params.output_dir}\n" }
