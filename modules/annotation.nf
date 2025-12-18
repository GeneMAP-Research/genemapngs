def getVCF() {
    return channel.fromFilePairs( params.vcf_dir + "*.vcf.{gz,gz.tbi}", size: 2 )
                  .ifEmpty { println "\nERROR: No VCF file found!\n" }
                  .map { vcf, index -> tuple(vcf, index) }
}

process getChromFromVcf() {
    tag "processing ${vcf}..."
    label 'bcftools'
    label 'variantCaller'
    //publishDir \
    //    path: "${params.output_dir}/gvcfs/", \
    //    mode: 'copy'
    input:
        tuple \
            path(vcf), \
            path(index)
    output:
        tuple \
            path("*.txt"), \
            path(vcf), \
            path(index)
    script:
        """
        bcftools \
            index \
            --stats \
            ${vcf} | \
        awk '{print \$1}' > chroms.list

        for chrom in \$(cat chroms.list); do echo \${chrom} > \${chrom}.txt; done
        """
}

process indexVcf() {
    tag "BCFTOOLS INDEX: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    input:
        tuple \
            val(chrom), \
            path(input_vcf)
    output:
        tuple \
            val(chrom), \
            path("${input_vcf}"), \
            path("${input_vcf}.tbi")
    script:
        """
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${input_vcf}
        """
}

process getChromFromIntervalList() {
    tag "processing ${interval}..."
    input:
        path(interval)
    output:
        tuple \
            path("*.chrom"), \
            path("${interval}")
    script:
        """
        echo ${interval} > \$(awk '{print \$1}' ${interval}).chrom
        """
}

process splitVcfPerInterval() {
    tag "processing ${interval.baseName}..."
    label 'bcftools'
    label 'variantCaller'
    //publishDir \
    //    path: "${params.output_dir}/gvcfs/", \
    //    mode: 'copy'
    input:
        tuple \
            val(chrom), \
            path(interval), \
            path(vcf), \
            path(index)
    output:
        path("${interval.baseName}.vcf.gz")
    script:
        """
        bcftools \
            view \
            -r \$(cat ${interval} | sed 's/ /:/1; s/ /-/1') \
            --threads ${task.cpus} \
            -Oz \
            -o ${interval.baseName}.vcf.gz \
            ${vcf}
        """
}

process validateVcfChunks() {
    tag "Validating VCF chunks..."
    input:
        path(vcf)
    output:
        path("*.validated.vcf.gz")
    script:
        """
        for vcf in ${vcf}; do
            vcfcontent=\$(zgrep -v '^#' \${vcf} | head -1 | cut -f1)
            if [[ ! \${vcfcontent} == "" ]]; then
                ln -s \${vcf} \$(basename \${vcf/.vcf.gz/.validated.vcf.gz})
            fi
        done
        """
}

process annovarGRCh37() {
    tag "processing ${vcfFile}"
    label 'annovar'
    label 'annovarMem'
    if(params.interval == 'NULL') {
        publishDir \
            path: "${params.output_dir}/annovar/", \
            mode: 'copy'
    } else {
        publishDir \
            path: "${params.output_dir}/tmp/"
    }
    input:
        path(vcfFile)
    output:
        path("${vcfFile.baseName}*multianno.{vcf.gz,txt.gz}")
    script:
        """
        table_annovar.pl \
            ${vcfFile} \
            ${params.annovarDB} \
            -buildver hg19 \
            -out "${vcfFile.baseName}" \
            -remove \
            -protocol refGene,knownGene,gwasCatalog,genomicSuperDups,cytoBand,exac03,avsnp151,abraom,dbnsfp33a,regsnpintron,esp6500siv2_all,SAS.sites.2015_08,EUR.sites.2015_08,EAS.sites.2015_08,AMR.sites.2015_08,ALL.sites.2015_08,AFR.sites.2015_08,nci60,clinvar_20210501,refGeneWithVer,revel,mitimpact24,intervar_20180118,hrcr1,gme,gnomad211_exome,gene4denovo201907,dbscsnv11,dbnsfp31a_interpro,cosmic70 \
            -operation g,g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
            -nastring '.' \
            -vcfinput \
            -polish \
            --dot2underline \
            --maxgenethread ${task.cpus} \
            --thread ${task.cpus}

        for i in ${vcfFile.baseName}*.{vcf,txt}; do bgzip -@ ${task.cpus} \${i}; done
        """
}

process minimalAnnovarGRCh37() {
    tag "processing ${vcfFile}"
    label 'annovar'
    label 'annovarMinimal'
    if(params.interval == 'NULL') {
        publishDir \
            path: "${params.output_dir}/", \
            mode: 'copy'
    } else {
        publishDir \
            path: "${workDir}/tmp/"
    }
    input:
        path(vcfFile)
    output:
        path("${vcfFile.baseName}*multianno.{vcf.gz,txt.gz}")
    script:
        """
        table_annovar.pl \
            ${vcfFile} \
            ${params.annovarDB} \
            -buildver hg19 \
            -out "${vcfFile.baseName}" \
            -remove \
            -protocol refGene,knownGene,cytoBand,gwasCatalog,avsnp151,gtexEqtlTissueWholeBlood,gtexEqtlTissueSpleen,gtexEqtlTissueXformedlymphocytes,gtexEqtlTissueWholeBlood,gtexEqtlTissueLiver,wgEncodeAffyRnaChipFiltTransfragsK562CellTotal,wgEncodeBroadHmmK562HMM_v2 \
            -operation g,g,r,r,f,r,r,r,r,r,r,r \
            -nastring '.' \
            -vcfinput \
            -polish \
            --dot2underline \
            --maxgenethread ${task.cpus} \
            --thread ${task.cpus}

        for i in ${vcfFile.baseName}*.{vcf,txt}; do bgzip -@ ${task.cpus} \${i}; done
        """
}

process annovarGRCh38() {
    tag "processing ${vcfFile}"
    label 'annovar'
    label 'annovarMem'
    publishDir path: "${params.output_dir}", mode: 'copy'
    input:
        path(vcfFile)
    output:
        path "${vcfFile.baseName}*multianno.{vcf.gz,txt.gz}"
    script:
        """
        table_annovar.pl \
            ${vcfFile} \
            ${params.annovarDB} \
            -buildver hg38 \
            -out "${vcfFile.baseName}" \
            -remove \
            -protocol refGene,knownGene,ucscGenePfam,cytoBand,keggPathway,dbnsfp47a,dbnsfp47a_interpro,dbscsnv11,intervar_20180118,cosmic70,exac03,gene4denovo201907,gnomad41_genome,gnomad41_exome,kaviar_20150923,ALL.sites.2015_08,gme,abraom,revel,avsnp151,clinvar_20240917,regsnpintron,gwasCatalog,GTEx_v8_eQTL,GTEx_v8_sQTL \
            -operation g,g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r \
            -nastring '.' \
            -vcfinput \
            -polish \
            --dot2underline \
            --maxgenethread ${task.cpus} \
            --thread ${task.cpus}

        for i in ${vcfFile.baseName}*.{vcf,txt}; do bgzip -@ ${task.cpus} \${i}; done
        """
}

process minimalAnnovarGRCh38() {
    tag "processing ${vcfFile}"
    label 'annovar'
    label 'annovarMinimal'
    if(params.interval == 'NULL') {
        publishDir \
            path: "${params.output_dir}/", \
            mode: 'copy'
    } else {
        publishDir \
            path: "${workDir}/tmp/"
    }
    input:
        path(vcfFile)
    output:
        path "${vcfFile.baseName}*multianno.{vcf.gz,txt.gz}"
    script:
        """
        table_annovar.pl \
            ${vcfFile} \
            ${params.annovarDB} \
            -buildver hg38 \
            -out "${vcfFile.baseName}" \
            -remove \
            -protocol refGene,knownGene,gwasCatalog,cytoBand,avsnp151 \
            -operation g,g,r,r,f \
            -nastring '.' \
            -vcfinput \
            -polish \
            --dot2underline \
            --maxgenethread ${task.cpus} \
            --thread ${task.cpus}

        for i in ${vcfFile.baseName}*.{vcf,txt}; do bgzip -@ ${task.cpus} \${i}; done
        """
}

