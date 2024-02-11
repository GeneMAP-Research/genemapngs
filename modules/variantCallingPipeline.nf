def getBamFileSet() {
    return channel.fromFilePairs( params.outputDir + '/bam/' + "*.{bam,bam.bai}", size: 2, flat: true )
                  .ifEmpty { error "\nERROR: Could not locate a file! \n" }
                  .map { bamName, bamFile, bamIndex -> tuple(bamName, bamFile, bamIndex) }
}

def getCramFileSet() {
    return channel.fromFilePairs( params.outputDir + '/cram/' + "*.{cram,cram.crai}", size: 2, flat: true )
                  .ifEmpty { error "\nERROR: Could not locate a file! \n" }
                  .map { bamName, bamFile, bamIndex -> tuple(bamName, bamFile, bamIndex) }
}

def getGvcfFiles() {
    return channel.fromPath( params.gvcf_dir + "*.{gvcf,g.vcf}.{gz,gz.tbi}" )
                  .flatten()
}

process haplotypeCaller() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'variantCaller'
    cache 'lenient'
    publishDir \
        path: "${params.outputDir}/gvcfs/", \
        mode: 'copy'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
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

process haplotypeCallerSpark() {
    tag "processing ${bamName}"
    beforeScript 'ulimit -c unlimited'
    label 'gatk'
    label 'variantCaller'
    cache 'lenient'
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

process combineGvcfsPerChromosome() {
    tag "Combining GVCFs to ${params.outPrefix}.gvcf.gz"
    label 'gatk'
    label 'variantCaller'
    input:
        tuple \
            val(interval), \
            path(gvcfList)
    output:
        publishDir path: "${params.outputDir}/vcf/"
        tuple \
            val(interval), \
            path("${interval}_${params.outPrefix}.g.vcf.{gz,gz.tbi}")
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
            -L ${interval} \
            --dbsnp ${params.dbsnp} \
            -O "${interval}_${params.outPrefix}.g.vcf.{gz,gz.tbi}"
        """
}

process getGenomicIntervals() {
    tag "Combining GVCFs to ${params.outPrefix}.gvcf.gz"
    label 'gatk'
    label 'smallMemory'
    cache 'lenient'
    input:
        path gvcfList
    output:
        path "*"
    script:
        """
        vcf=\$(ls *.g.vcf.gz | head -1)

        zgrep 'contig' \$(readlink \${vcf}) | \
            cut -f1 -d',' | \
            cut -f3 -d'=' > intervals.list

        for interval in \$(cat intervals.list); do
            touch \$interval
        done
        """
}

process createGenomicsDb() {
    tag "Combining GVCFs to ${params.outPrefix}.gvcf.gz"
    label 'gatk'
    label 'variantCaller'
    cache 'lenient'
    input:
        path gvcfList
    output:
        publishDir path: "${params.outputDir}/vcf/"
        path "${params.outPrefix}-workspace"
    script:
        """
        for file in ${gvcfList}; do
            if [ \${file##*.} != "tbi" ]; then
                echo "-V \${file}";
            fi
        done > gvcf.list

        vcf=\$(ls *.g.vcf.gz | head -1)

        zgrep 'contig' \$(readlink \${vcf}) | \
            cut -f1 -d',' | \
            cut -f3 -d'=' > intervals.list

        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            GenomicsDBImport \
            -R ${params.fastaRef} \
            --arguments_file gvcf.list \
            -L intervals.list \
            --genomicsdb-workspace-path ${params.outPrefix}-workspace
        """
}

process createGenomicsDbPerChromosome() {
    tag "Combining GVCFs to ${interval}_${params.outPrefix}-workspace"
    label 'gatk'
    label 'variantCaller'
    input:
        tuple \
            val(interval), \
            path(gvcfList)
    output:
        publishDir path: "${params.outputDir}/vcf/"
        tuple \
            val(interval), \
            path("${interval}_${params.outPrefix}-workspace")
    script:
        """
        for file in ${gvcfList}; do
            if [ \${file##*.} != "tbi" ]; then
                echo "-V \${file}";
            fi
        done > gvcf.list

        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            GenomicsDBImport \
            -R ${params.fastaRef} \
            --arguments_file gvcf.list \
            -L ${interval} \
            --genomicsdb-workspace-path ${interval}_${params.outPrefix}-workspace
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

process genotypeGvcfs() {
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

process callVariantsFromGenomicsDB() {
    tag "Writing genotypes to ${interval}_${params.outPrefix}.vcf.gz"
    label 'gatk'
    label 'variantCaller'
    input:
        tuple \
            val(interval), \
            path(gvcf)
    output:
        publishDir path: "${params.outputDir}/vcf/", mode: 'copy'
        path "${interval}_${params.outPrefix}.vcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            GenotypeGVCFs \
            -R ${params.fastaRef} \
            --dbsnp ${params.dbsnp} \
            -V gendb://${gvcf} \
            -O "${interval}_${params.outPrefix}.vcf.gz"
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

process dysguCallSvs() {
    tag "Writing genotypes to ${bamName}.vcf.gz"
    label 'dysgu'
    label 'dysgu_caller'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        publishDir path: "${params.outputDir}/vcf/", mode: 'copy'
        path "${bamName}.vcf.gz"
    script:
        """
        dysgu \
            run \
            -p ${task.cpus} \
            ${params.fastaRef} \
            -x temp \
            -v2 \
            ${bamFile} | \
            bgzip -c > ${bamName}.vcf.gz \
        """
}

process indexVcf() {
    tag "processing ${vcf}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path(vcf)
    output:
        tuple \
            path(vcf), \
            path("${vcf}.tbi")
    script:
        """
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${vcf}
        """
}

process dysguMergeVcfs() {
    tag "Writing genotypes to ${params.outPrefix}_dysgu_sv.vcf.gz"
    label 'dysgu'
    label 'dysgu_caller'
    input:
        path(vcfs)
    output:
        publishDir path: "${params.outputDir}/vcf/", mode: 'copy'
        path "${params.outPrefix}_dysgu_sv.vcf.gz"
        """
        dysgu \
          merge \
          *.vcf.gz \
          -p ${task.cpus} | \
          bgzip -c > ${params.outPrefix}_dysgu_sv.vcf.gz
        """
}

process getBamChuncks() {
    tag "Splitting list of BAM files into chunks of 5..."
    input:
        path bam
    output:
        path("${params.outPrefix}_bamchunk_aaaaaaa*.txt")
    script:
        """
        readlink \$(ls *.bam) > ${params.outPrefix}_bamlist.txt
        split \
            -l 5 \
            -a 8 \
            --additional-suffix .txt \
            ${params.outPrefix}_bamlist.txt \
            ${params.outPrefix}_bamchunk_
        """
}

process mantaCallSvs() {
    tag "Running MANTA..."
    label 'manta'
    label 'dysgu_caller'
    publishDir \
        path: "${params.outputDir}/vcf/", \
        mode: 'copy'
    input:
        path bamlist
    output:
        path("*.vcf*")
    script:
        """
        bamlist=${bamlist}
        cat ${bamlist} | \
            awk '{print "--bam",\$1,"\\\\"}' | \
            sed "1 i --referenceFasta ${params.fastaRef} \$(if [[ ${params.exome} == true ]]; then echo --exome; fi)" | \
            sed 's/--exome/\\n--exome \\\\/g' | \
            sort -g | \
            sed '1 i configManta.py \\\\' \
            > \${bamlist/.txt/.sh}

        chmod 755 \${bamlist/.txt/.sh}

        ./\${bamlist/.txt/.sh}

        ./MantaWorkflow/runWorkflow.py \
            -j ${task.cpus}

        # ensure that there are no VCF files in the workdir
        [ -e "*.vcf*" ] && rm *.vcf*

        # get VCF results from manta workflow dir to nextflow workdir
        cp ./MantaWorkflow/results/variants/* .

        # rename manta results to include chunk filename
        for result in *.vcf*; do
            mv \${result} \${bamlist/.txt/}_\${result}
        done
        """
}

process mergeMantaDiploidSvCalls() {
    tag "mergeing diploid SV files..."
    label 'bcftools'
    label 'dysgu_caller'
    publishDir \
        path: "${params.outputDir}/vcf/", \
        mode: 'copy'
    input:
        path manta_calls
    output:
        path("*.vcf.gz*")
    script:
        """
        bcftools \
            merge \
            -m all \
            --threads ${task.cpus} \
            -Oz \
            \$(ls *_diploidSV.vcf.gz | tr '\n' ' ') | \
            tee ${params.outPrefix}_manta_diploidSV.vcf.gz | \
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            -o ${params.outPrefix}_manta_diploidSV.vcf.gz.tbi
        """
}

process mergeMantaCandidateSmallIndelCalls() {
    tag "mergeing candidate small indels files..."
    label 'bcftools'
    label 'dysgu_caller'
    publishDir \
        path: "${params.outputDir}/vcf/", \
        mode: 'copy'
    input:
        path manta_calls
    output:
        path("*.vcf.gz*")
    script:
        """
        bcftools \
            merge \
            -m all \
            --threads ${task.cpus} \
            -Oz \
            \$(ls *_candidateSmallIndels.vcf.gz | tr '\n' ' ') | \
            tee ${params.outPrefix}_manta_candidateSmallIndels.vcf.gz | \
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            -o ${params.outPrefix}_manta_candidateSmallIndels.vcf.gz.tbi
        """
}

process mergeMantaCandidateSvCalls() {
    tag "mergeing candidate SV files..."
    label 'bcftools'
    label 'dysgu_caller'
    publishDir \
        path: "${params.outputDir}/vcf/", \
        mode: 'copy'
    input:
        path manta_calls
    output:
        path("*.vcf.gz*")
    script:
        """
        bcftools \
            merge \
            -m all \
            --threads ${task.cpus} \
            -Oz \
            \$(ls *_candidateSV.vcf.gz | tr '\n' ' ') | \
            tee ${params.outPrefix}_manta_CandidateSvs.vcf.gz | \
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            -o ${params.outPrefix}_manta_CandidateSvs.vcf.gz.tbi
        """
}
