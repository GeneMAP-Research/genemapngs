def getVcf() {
    return channel.fromPath( params.vcf_dir + "*.vcf.gz" )
}

def getThousandGenomesReference() {
    return channel.fromFilePairs( params.ref_dir + "chr*1kg.phase3.v5a.vcf.{gz,gz.tbi}", size: 2 )
                  .ifEmpty { error "\nAn error occurred! Please check that the reference file and its index '.tbi' exist...\n" }
	   	  .map { chr, ref_file -> 
			 tuple( chr.replaceFirst(/chr/,""), ref_file.first(), ref_file.last())
	    	  }
}

process getVcfIndex() {
    tag "BCFTOOLS INDEX: ${input_vcf}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path input_vcf
    output:
        path "*.tbi"
    script:
        """
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${input_vcf}
        """
}

process vqsrSnp() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'vqsr'
    input:
        path(input_vcf)
        path(vcf_index)
    output:
        publishDir path: "${params.outputDir}/vqsr-tables/", mode: 'copy'
        tuple \
              path("${params.outPrefix}.snp.recal"), \
              path("${params.outPrefix}.snp.recal.idx"), \
              path("${params.outPrefix}.snp.tranches")
    script:
        """
        gatk \
           --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
           VariantRecalibrator \
           -R ${params.fastaRef} \
           -tranche 100.0 -tranche 99.95 -tranche 99.9 \
           -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
           -tranche 95.0 -tranche 94.0 \
           -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
           -V ${input_vcf} \
           --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
           --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.snpRef} \
           --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.omniRef} \
           --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp} \
           -an QD \
           -an MQ \
           -an MQRankSum \
           -an ReadPosRankSum \
           -an FS \
           -an SOR \
           -mode SNP \
           -O ${params.outPrefix}.snp.recal \
           --tranches-file ${params.outPrefix}.snp.tranches \
           --rscript-file ${params.outPrefix}.snp.plots.R
        """
}

process vqsrIndel() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'vqsr'
    input:
        path(input_vcf)
        path(vcf_index)
    output:
        publishDir path: "${params.outputDir}/vqsr-tables/", mode: 'copy'
        tuple \
              path("${params.outPrefix}.indel.recal"), \
              path("${params.outPrefix}.indel.recal.idx"), \
              path("${params.outPrefix}.indel.tranches")
    script:
        """
        gatk \
           --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
           VariantRecalibrator \
           -R ${params.fastaRef} \
           -V ${input_vcf} \
           -tranche 100.0 -tranche 99.95 -tranche 99.9 \
           -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
           -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
           -tranche 92.0 -tranche 91.0 -tranche 90.0 \
           --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
           --resource:1000G,known=false,training=true,truth=false,prior=12.0 ${params.indelsRef} \
           --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.omniRef} \
           --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp} \
           -an QD \
           -an DP \
           -an MQRankSum \
           -an ReadPosRankSum \
           -an FS \
           -an SOR \
           -mode INDEL \
           --max-gaussians 4 \
           -O ${params.outPrefix}.indel.recal \
           --tranches-file ${params.outPrefix}.indel.tranches \
           --rscript-file ${params.outPrefix}.indel.plots.R
        """
}

process applyVqsrSnp() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'applyVqsr'
    input:
        path(input_vcf)
        path(vcf_index)
        tuple \
            path(recal_table), \
            path(recal_index), \
            path(tranches)   
    output:
        publishDir path: "${params.outputDir}/vqsr/", mode: 'copy'
        path "${params.outPrefix}.snp.vqsr.vcf.gz"
    script:
        """
        gatk \
            --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            ApplyVQSR \
            -V ${input_vcf} \
            --recal-file ${recal_table} \
            --mode SNP \
            --truth-sensitivity-filter-level 99.9 \
            --create-output-variant-index true \
            --tranches-file ${tranches} \
            -R ${params.fastaRef} \
            -O ${params.outPrefix}.snp.vqsr.vcf.gz
        """
}

process applyVqsrIndel() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'applyVqsr'
    input:
        path(input_vcf)
        path(vcf_index)
        tuple \
            path(recal_table), \
            path(recal_index), \
            path(tranches)   
    output:
        publishDir path: "${params.outputDir}/vqsr/", mode: 'copy'
        path "${params.outPrefix}.indel.vqsr.vcf.gz"
    script:
        """
        gatk \
            --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            ApplyVQSR \
            -V ${input_vcf} \
            --recal-file ${recal_table} \
            --mode INDEL \
            --truth-sensitivity-filter-level 99.9 \
            --create-output-variant-index true \
            --tranches-file ${tranches} \
            -R ${params.fastaRef} \
            -O ${params.outPrefix}.indel.vqsr.vcf.gz
        """
}

process mergeVCFs() {
    tag "merging VCFs"
    label 'bcftools'
    label 'mediumMemory'
    input:
        tuple \
            path(snp_vcf), \
            path(indel_vcf)
    output:
        publishDir path: "${params.outputDir}/vqsr/", mode: 'copy'
        path"${params.outPrefix}.snp.indel.vqsr.vcf.{gz,gz.tbi}"
    script:
        """
        for i in $snp_vcf $indel_vcf; do bcftools index -ft --threads ${task.cpus} \$i; done

        bcftools \
            concat \
            -a -d all \
            --threads ${task.cpus} \
            -Oz \
            $snp_vcf \
            $indel_vcf | \
            tee ${params.outPrefix}.snp.indel.vqsr.vcf.gz | \
        bcftools index \
            --threads ${task.cpus} \
            -ft \
            --output "${params.outPrefix}.snp.indel.vqsr.vcf.gz.tbi"
        """
}

process filterGatkCalls() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    input:
        path input_vcf
    output:
        publishDir path: "${params.outputDir}/filtered/", mode: 'copy'
        path("${input_vcf.baseName}.filtered.vcf.gz")
    script:
        """
        bcftools \
            view \
            -i \'FILTER=="PASS"\' \
            --threads ${task.cpus} \
            ${input_vcf} | \
            bcftools \
                view \
                -i \'INFO/DP>=${params.minDP}\' \
                --threads ${task.cpus} | \
                bcftools \
                    view \
                    -i \'GQ>=${params.minGQ}\' \
                    --threads ${task.cpus} \
                    -Oz | \
                    tee "${input_vcf.baseName}.filtered.vcf.gz" | \
                bcftools index \
                    --threads ${task.cpus} \
                    -ft \
                    --output "${input_vcf.baseName}.filtered.vcf.gz.tbi"
        """
}

process filterGlnexusCalls() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    input:
        path(input_vcf)
        path(vcf_index)
    output:
        publishDir path: "${params.outputDir}/filtered/", mode: 'copy'
        tuple \
            path("${input_vcf.baseName}.filtered.vcf.gz"), \
            path("${input_vcf.baseName}.filtered.vcf.gz.tbi")
    script:
        """
        bcftools \
            view \
            -i \'FILTER!="MONOALLELIC"\' \
            --threads ${task.cpus} \
            ${input_vcf} | \
            bcftools \
                view \
                -i \'DP>=${params.minDP}\' \
                --threads ${task.cpus} | \
                bcftools \
                    view \
                    -i \'GQ>=${params.minGQ}\' \
                    --threads ${task.cpus} \
                    -Oz | \
                    tee "${input_vcf.baseName}.filtered.vcf.gz" | \
                bcftools index \
                    --threads ${task.cpus} \
                    -ft \
                    --output "${input_vcf.baseName}.filtered.vcf.gz.tbi"
        """
}

process splitMultiallelicSnvs() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index)
    output:
        tuple \
            path("${input_vcf.baseName}-tmp.vcf.gz"), \
            path("${input_vcf.baseName}-tmp.vcf.gz.tbi")
    script:
        """
        bcftools \
            norm \
            -m -any \
            --threads ${task.cpus} \
            -Oz \
            ${input_vcf} | \
            tee "${input_vcf.baseName}-tmp.vcf.gz" | \
        bcftools index \
            --threads ${task.cpus} \
            -ft \
            --output "${input_vcf.baseName}-tmp.vcf.gz.tbi"            
        """
}

process leftnormalizeSnvs() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index)
    output:
        publishDir path: "${params.outputDir}/filtered/", mode: 'copy'
        path "${input_vcf.baseName}-filtered-leftnorm.vcf.{gz,gz.tbi}"
    script:
        """
        bcftools \
            norm \
            -f ${params.fastaRef} \
            --threads ${task.cpus} \
            -Oz \
            ${input_vcf} | \
        bcftools \
            view \
            -c ${params.minAC} \
            --threads ${task.cpus} \
            -Oz | \
            tee "${input_vcf.baseName}-filtered-leftnorm.vcf.gz" | \
        bcftools index \
            --threads ${task.cpus} \
            -ft \
            --output "${input_vcf.baseName}-filtered-leftnorm.vcf.gz.tbi"
        """
}

process getCleanVcf() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index)
    output:
        publishDir path: "${params.outputDir}/filtered/", mode: 'copy'
        path "${input_vcf.baseName}-filtered-leftnorm-clean.vcf.{gz,gz.tbi}"
    script:
        """
        bcftools \
            view \
            --threads ${task.cpus} \
            -r \$(echo chr{1..22}, chrX | sed 's/[[:space:]]//g') \
            -Oz \
            "${input_vcf}" | 
            tee "${input_vcf.baseName}-filtered-leftnorm-clean.vcf.gz" | \
            bcftools index \
            --threads ${task.cpus} \
            -ft \
            --output "${input_vcf.baseName}-filtered-leftnorm-clean.vcf.gz.tbi"
        """
}

process getVcfStats() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    input:
        tuple \
            path(input_vcf), \
            path(vcf_index)
    output:
        publishDir path: "${params.outputDir}/vcfstats/", mode: 'copy'
        path "${input_vcf.baseName}.vcfstats.txt"
    script:
        """
        bcftools \
            stats \
            -F "${params.fastaRef}" \
            -s - \
            "${input_vcf}" > \
            "${input_vcf.baseName}.vcfstats.txt"
        """
}

process plotVcfStats() {
    tag "Processing ${vcfstat}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path vcfstat
    output:
        publishDir path: "${params.outputDir}_${vcfstat.baseName}", mode: 'copy'
        path "./*"
    script:
        """
        mkdir -p "${params.outputDir}${vcfstat.baseName}"
        plot-vcfstats \
           -p . \
           "${vcfstat}"
        """
}

