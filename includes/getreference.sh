#!/usr/bin/env bash

if [ $# -lt 2 ]; then
  echo -e "\nUsage: getreference.sh [reference directory] [build]\n"
  echo -e "build --> hg19, hg38, t2t\n"

else
  reference_dir=$1; build=$2
  home=$(pwd)

  mkdir -p ${reference_dir}/${build}/{refgenome,gatkbundles}
  cd ${reference_dir}/${build}

  if [[ $build == "hg38" ]]; then
    #hg38

    #wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.{gz,gz.tbi}

    echo hg38
    cd gatkbundles && gatkdir=$(pwd) && \
    wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.{gz,gz.tbi} && \
    wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.{gz,gz.tbi} && \
    wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.{gz,gz.tbi} && \
    wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.{gz,gz.tbi} && \
    wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.{vcf,vcf.idx} && \
    cd ../refgenome && refdir=$(pwd) && \
    wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.{dict,fasta,fastq.fai} && \
    wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.{alt,amb,ann,bwt,pac,sa} && \
    cd $home

    if [ $? -eq 0 ]; then
      cat <<EOF > configs/references/hg38.config
params {
    ref_dir = '${refdir}/'
    fastaRef = "\${params.ref_dir}Homo_sapiens_assembly38.fasta"
    gatk_dir = '${gatkdir}/'
    dbsnp = "\${params.gatk_dir}Homo_sapiens_assembly38.dbsnp138.vcf"
    indelsRef = "\${params.gatk_dir}Homo_sapiens_assembly38.known_indels.vcf.gz"
    hapmap = "\${params.gatk_dir}hapmap_3.3.hg38.vcf.gz"
    omniRef = "\${params.gatk_dir}1000G_omni2.5.hg38.vcf.gz"
    snpRef = "\${params.gatk_dir}1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    annovarDB = ''
    buildVersion = 'hg38'
}
EOF
    fi

  elif [[ $build == "hg19" ]]; then
    #hg19 ~ https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/

    #wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase3_v4_20130502.sites.vcf.{gz,gz.tbi}
    #wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.{gz,idx.gz}

    echo hg19
    cd refgenome && refdir=$(pwd) && \
    wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz && \
    wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.bwa_index.tar.gz && \
    gunzip hg19.p13.plusMT.no_alt_analysis_set.fa.gz && \
    tar zxvf hg19.p13.plusMT.no_alt_analysis_set.bwa_index.tar.gz && \
    cd ../gatkbundles && gatkdir=$(pwd) && \
    wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.{gz,idx.gz} && \
    wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.{gz,idx.gz} && \
    wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.{gz,idx.gz} && \
    wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.{gz,idx.gz} && \
    wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf.{gz,idx.gz} && \
    cd $home

    if [ $? -eq 0 ]; then
      cat <<EOF > configs/references/hg19.config
params {
    ref_dir = '${refdir}/'
    fastaRef = "\${params.ref_dir}/hg19.p13.plusMT.no_alt_analysis_set.fa"
    gatk_dir = '${gatkdir}/'
    dbsnp = "\${params.gatk_dir}/dbsnp_138.hg19.vcf.gz"
    indelsRef = "\${params.gatk_dir}/1000G_phase1.indels.hg19.sites.vcf.gz"
    hapmap = "\${params.gatk_dir}/hapmap_3.3.hg19.sites.vcf.gz"
    omniRef = "\${params.gatk_dir}/1000G_omni2.5.hg19.sites.vcf.gz"
    snpRef = "\${params.gatk_dir}/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz"
    annovarDB = ''
    buildVersion = 'hg19'
}
EOF
    fi
    cd $home
  elif [[ $build == "t2t" ]]; then
    #t2t
    echo t2t

    cd -
  fi
fi
