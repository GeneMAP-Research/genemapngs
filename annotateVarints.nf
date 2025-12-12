#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getGenomicIntervalList;
    collectIntervalsPerChromosome;
    concatPerChromIntervalVcfs;
    concatPerChromosomeVcfs
} from "${projectDir}/modules/variantCallingPipeline.nf"

include {
    getVCF;
    indexVcf;
    splitVcfPerInterval;
    getChromFromIntervalList;
    validateVcfChunks;
    getChromFromVcf;
    annovarGRCh37;
    annovarGRCh38;
    minimalAnnovarGRCh37;
    minimalAnnovarGRCh38;
} from "${projectDir}/modules/annovarAnnotation.nf"

include {
    getVcf;
    getVcfIndex;
    splitMultiallelicSnvs;
    leftnormalizeSnvs;
} from "${projectDir}/modules/vcfQualityMetricsAndFiltering.nf"

workflow {
    println "\nANNOVAR annotation starts here\n"

    msg = "\nERROR: You must select a Build Version! Options are hg19 and hg38\n" 

    getVcf()
        .map { vcfname, vcf, index -> vcf }
        .set { vcf }

    if(params.left_norm == true) {
        //indexed_vcf = getVcfIndex(vcf)
        multisplit = splitMultiallelicSnvs(vcf)
        leftnorm = leftnormalizeSnvs(multisplit)
        //vcf = leftnorm.map { vcf, index -> tuple(vcf, index) }
        vcf = leftnorm
    }

    if(!(params.interval == "NULL")) {
        getChromFromVcf(vcf)
            .flatMap { chroms, vcf, index -> 
                chroms.collect { chrom ->
                    tuple("${chrom.simpleName}", vcf, index)
                }
            }
            .set { chrom_vcf }

        //indexVcf(chrom_vcf).view()
        //   .set { chrom_vcf_index }

        getGenomicIntervalList()
           .flatten()
           .set { interval_list }

        getChromFromIntervalList(interval_list)
           .map { chrom, interval -> tuple("${chrom.simpleName}", interval) }
           .groupTuple()
           .join(chrom_vcf)
           .flatMap { chr, intervals, vcf, index -> 
               intervals.collect { intval ->
                   tuple(chr, intval, vcf, index)
               }
           }.view()
           .set { split_input }

        splitVcfPerInterval(split_input)
            .collect()
            .set { vcfs }

        validateVcfChunks(vcfs)
            .collect()
            .flatten()
            .set { vcf }
    }

    if (params.minimal == true) {
       if (params.buildVersion == "hg19") {
           annotated = minimalAnnovarGRCh37(vcf).collect()
       }
       else if (params.buildVersion == "hg38") {
           annotated = minimalAnnovarGRCh38(vcf).collect()
       }
       else { error: println msg }
    }
    else {
       if (params.buildVersion == "hg19") {
           annotated = annovarGRCh37(vcf).collect()
       }
       else if (params.buildVersion == "hg38") {
           annotated = annovarGRCh38(vcf).collect()
       }
       else { error: println msg }
    }

    if(!(params.interval == "NULL")) {
        vcfs_per_chrom_list = collectIntervalsPerChromosome(annotated).flatten()
        vcfs_per_chrom = concatPerChromIntervalVcfs(vcfs_per_chrom_list).collect().view()
        //concatPerChromosomeVcfs(vcfs_per_chrom).view()
    }
}

