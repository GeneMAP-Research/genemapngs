def getBamFileSet() {
    return channel.fromFilePairs( params.outputDir + '/bam/' + "*.{bam,bam.bai}", size: 2, flat: true )
                  .ifEmpty { error "\nERROR: Could not locate a file! \n" }
                  .map { bamName, bamFile, bamIndex -> tuple(bamName, bamFile, bamIndex) }
}

process callVariants() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'variantCaller'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        publishDir path: "${params.outputDir}/gvcfs/"
        path "${bamName}.g.vcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            HaplotypeCaller \
            -I ${bamFile} \
            -R ${params.fastaRef} \
            -ERC GVCF \
            -OBI 'false' \
            -O "${bamName}.g.vcf.gz" 
        """
}

//            --java-options "-XX:ConcGCThreads=${task.cpus} -Xmx${task.memory.toGiga()}g" \

process callVariantsSpark() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'variantCaller'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        publishDir path: "${params.outputDir}/gvcfs/"
        path "${bamName}.g.vcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            HaplotypeCallerSpark \
            -I ${bamFile} \
            -R ${params.fastaRef} \
            --native-pair-hmm-threads ${task.cpus} \
            -ERC GVCF \
            -OBI 'false' \
            -O "${bamName}.g.vcf.gz" \
            -- \
            --spark-master local[${task.cpus}]
        """
}

process combineGvcfs() {
    tag "Combining GVCFs to ${params.outPrefix}.gvcf.gz"
    label 'gatk' 
    label 'variantCaller'
    input:
        path gvcfList
    output:
        publishDir path: "${params.outputDir}/vcf/"
        path "${params.outPrefix}.g.vcf.{gz,gz.tbi}"
    script:
        """
        for file in ${gvcfList}; do
            if [ \${file##*.} != "tbi" ]; then
                echo "-V \${file}";
            fi
        done > gvcf.list

        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            CombineGVCFs \
            -R ${params.fastaRef} \
            --arguments_file gvcf.list \
            --dbsnp ${params.dbsnp} \
            -O "${params.outPrefix}.g.vcf.gz"
        """
}

process getPedFile() {
    tag "Writing PED file"
    label 'bcftools'
    label 'smallMemory'
    echo true
    input:
        tuple \
            path(gvcf), \
            path(gvcfIndex)
    output:
        publishDir path: "${params.outputDir}/ped/"
        path "${params.outPrefix}.tmp.ped"
    script:
        """
        if [ ${params.ped} == NULL ]; then
           bcftools query -l ${gvcf} | \
               awk '{print \$1,\$1,"0","0","-9","-9"}' \
               > ${params.outPrefix}.ped
           ped="${params.outPrefix}.ped"
        else
           ped="${params.ped}"
        fi
        cp \${ped} ${params.outPrefix}.tmp.ped
        """
}

process joinCallVariants() {
    tag "Writing genotypes to ${params.outPrefix}.vcf.gz"
    label 'gatk'
    label 'variantCaller'
    input:
        tuple \
            path(gvcf), \
            path(gvcfIndex), \
            path(ped)
    output:
        publishDir path: "${params.outputDir}/vcf/", mode: 'copy'
        path "${params.outPrefix}.vcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            GenotypeGVCFs \
            -R ${params.fastaRef} \
            --dbsnp ${params.dbsnp} \
            -ped ${ped} \
            -V ${gvcf} \
            -O "${params.outPrefix}.vcf.gz"
        """
}

process deepVariantCaller() {
    tag "Writing genotypes to ${params.outPrefix}.vcf.gz"
    label 'deepvariant'
    label 'deepv_caller'
    beforeScript = 'module load python/anaconda-python-3.7'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        publishDir path: "${params.outputDir}/vcf/", mode: 'copy'
        path "${bamName}.g.vcf.{gz,gz.tbi}"
    script:
        """
        run_deepvariant \
            --model_type=WGS \
            --ref=${params.fastaRef} \
            --reads=${bamFile} \
            --output_gvcf=${bamName}.g.vcf.gz \
            --output_vcf=${bamName}.vcf.gz \
            --num_shards=${task.cpus}
        """
}

process glnexusJointCaller() {
    tag "Writing genotypes to ${params.outPrefix}.vcf.gz"
    label 'glnexus'
    label 'nexus_caller'
    input:
        path gvcfList
    output:
        publishDir path: "${params.outputDir}/vcf/"
        path "${params.outPrefix}_glnexus.bcf"
    script:
        """
        for file in ${gvcfList}; do
            if [ \${file##*.} != "tbi" ]; then
                echo "\${file}";
            fi
        done | sort -V > gvcf.list

        if [ -d "GLnexus.DB" ]; then rm -rf GLnexus.DB; fi

        glnexus_cli \
            --config DeepVariant \
            --list gvcf.list \
            --threads ${task.cpus} \
            > "${params.outPrefix}_glnexus.bcf"
        """
}

process convertBcfToVcf() {
    tag "Writing genotypes to ${bcf_file.baseName}.vcf.gz"
    label 'bcftools'
    label 'nexus_caller'
    input:
        path bcf_file
    output:
        publishDir path: "${params.outputDir}/vcf/", mode: 'copy'
        path "${bcf_file.baseName}.vcf.gz"
    script:
        """
        bcftools \
            view \
            --threads ${task.cpus} \
            -Oz \
            -o ${bcf_file.baseName}.vcf.gz \
            ${bcf_file}
        """
}
