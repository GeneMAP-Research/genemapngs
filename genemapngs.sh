#!/usr/bin/env sh

#--- genemapgwas workflow wrapper ---#

function usage() {
   echo """
   ===================================================================
   GeneMAP-NGS ~ a wrapper for the nextflow-based genemapngs workflow
   ===================================================================

   Usage: genemapngs <workflow> <profile> [options] ...

           workflows:
           ---------
               test: Run test to see if workfow installed correctly.
                 qc: Check FASTQ or Alignment (BAM/CRAM) quality.
               trim: Trim adapters and poor quality bases from reads.
              align: Align/map reads to reference and post-alignment processing.
            varcall: Perform variant calling (both single and joint sample) in one run.
           svarcall: Perform only sinlge sample variant calling to generate gVCF files.
           jvarcall: Perform only joint (multi-sample) variant calling with pre-existing gVCF files.
          varfilter: Filter variant calls in VCF/BCF files.


           profiles: <executor>,<container>,<reference>
           ---------
          executors: local, slurm
         containers: singularity, aptainer, docker
          reference: hg19, hg38, t2t


           examples:
           ---------
         genemapngs qc slurm,singularity,hg38 [options]
         genemapngs align local,singularity,hg19 --bfile BEDFILE --out MYOUT --outdir MYPATH --pheno_file MYPHENO
   """
}


#####################################################################################################################
function checkprofile() {
   #profile was passed as only argument 
   #so it takes position $1 here
   if [[ "$1" == "" ]]; then
      echo "ERROR: please specify a profile to use!";
      usage;
      exit 1;
   elif [[ $1 == -* ]]; then
      echo "ERROR: please specify a valid profile!";
      usage;
      exit 1;
   else
      local profile="$1"
   fi
}

function check_ftype() {
   ftype=$1
   if [[ $ftype == -* ]]; then
      echo "ERROR: Invalid paramter value for option '--ftype'"
      exit 1;
   fi
}

function check_required_params() {
   for params_vals in $@; do
      #get each param and its value as an array
      param_val=( $(echo ${params_vals} | sed 's/,/ /g') )
      
      #slice the array to its consituent params and values
      param=${param_val[0]}
      val=${param_val[1]}

      #now check each param and its value
      if [[ $val == -* ]] || [[ $val == NULL ]]; then
         echo "ERROR: Invalid paramter value for option '--${param}'";
         break;
         exit 1;
      fi
   done
}

function check_optional_params() {
   for params_vals in $@; do
      #get each param and its value as an array
      param_val=( $(echo ${params_vals} | sed 's/,/ /g') )

      #slice the array to its consituent params and values
      param=${param_val[0]}
      val=${param_val[1]}

      #now check each param and its value
      if [[ $val == -* ]]; then
         echo "ERROR: Invalid paramter value for option '--${param}'";
         break;
         exit 1;
      fi
   done
}

function check_output_dir() {
   output_dir=$1
   if [[ $output_dir == -* ]]  || [[ $output_dir == NULL ]]; then
      #output_dir="${input_dir}/../" 
      if [ -d ${input_dir} ]; then
         output_dir="${input_dir}/../"
      elif [ -d ${gvcf_dir} ]; then
         output_dir="${gvcf_dir}/../"
      #elif [ -d ${genomicsdb_workspace_dir} ]; then
      #   output_dir="${genomicsdb_workspace_dir}/../"
      fi
      if [[ $output_dir == NULL* ]]; then
         echo "ERROR: Invalid paramter value for option '--output_dir'"
         exit 1;
      fi
   fi
}

function check_resources() {
   threads=$1
   njobs=$2
   if [[ $threads == -* ]]; then
      echo "ERROR: Invalid paramter value for option '--threads'"
      exit 1;
   fi
   if [[ $njobs == -* ]]; then
      echo "ERROR: Invalid paramter value for option '--njobs'"
      exit 2;
   fi
}

function setglobalparams() {
#- create the project nextflow config file
echo """includeConfig \"\${projectDir}/nextflow.config\"
includeConfig \"\${projectDir}/configs/profile-selector.config\"
includeConfig \"\${projectDir}/configs/resource-selector.config\"
"""
}

##################################################### USAGE #########################################################
function qcusage() {
   echo -e "\nUsage: genemapngs qc <profile> [options] ..."
   echo """
           options:
           --------
           --ftype              : Input file type; FASTQ, BAM, CRAM [default: FASTQ].
           --input_dir          : (required) Path to FASTQ/BAM/CRAM files.
           --output_dir         : (optional) Results will be saved to parent of input directory ['input_dir/../'].
           --threads            : number of computer cpus to use [default: 4].
           --njobs              : (optional) number of jobs to submit at once [default: 4]
           --help               : print this help message.
   """
}

function trimusage() {
   echo -e "\nUsage: genemapngs trim <profile> [options] ..."
   echo """
           options:
           --------
           --ftype              : Input file type; FASTQ, BAM, CRAM [default: FASTQ].
           --input_dir          : (required).
           --output_dir         : (optional) Results will be saved to parent of input directory ['input_dir/../'].
           --trimmer            : (optional) Options: trimmomatic, trimgalore [default: trimgalore].
           --adapter            : required if 'trimmomatic' selected [default: NP].
 
                         options:
                                  NP   --> NexteraPE-PE.fa
                                  T3U  --> TruSeq3-PE-2.fa [Illumina universal]
                                  T2P  --> TruSeq2-PE.fa
                                  T2S  --> TruSeq2-SE.fa
                                  T3P  --> TruSeq3-PE.fa
                                  T3S  --> TruSeq3-SE.fa

                              NB: 'trimgalore' will auto-detect adapters. Hence suitable to process
                                  samples from different sequencing companies in one batch.

           --min_length         : minimum read leangth to keep [default: 36].
           --headcrop           : number of bases to remove from the start of reads [default: 5].
           --crop               : number of bases to remove from the end of reads [default: 5].
           --threads            : number of computer cpus to use [default: 8].
           --njobs              : (optional) number of jobs to submit at once [default: 10]
           --help               : print this help message.
   """
}

function alignusage() {
   echo -e "\nUsage: genemapngs align <profile> [options] ..."
   echo """
           options:
           --------
           
           --ftype		: Input file type; FASTQ, BAM, CRAM [default: FASTQ].
           --input_dir          : (required) Path to FASTQ/BAM/CRAM files.
           --aligner		: Alignment tool (required); BWA, DRAGMAP [default: BWA].
           --output_dir		: (optional) [results will be saved to parent of input directory].
           --se			: If FASTQs are single-end (paired-end is assumed by default).
           --wgs		: If data is ehole-genome sequence (it runs whole exome - wes - by default)
           --dup_marker         : Duplicate marker tool (optional); sambamba, samtools [default: sambamba]
           --remove_dup         : Whether to remove duplicates (optional): true, false [default: false]
           --spark		: Whether to use GATK spark mode for post-alignment processing (it multi-threads). Is not used by default.
           --threads		: number of computer cpus to use  [default: 11].
           --njobs              : (optional) number of jobs to submit at once [default: 10]
           --help		: print this help message.
   """
}

function varcallusage() {
   echo -e "\nUsage: genemapngs varcall <profile> [options] ..."
   echo """
           options:
           --------

           --alignment_dir      : (required) Path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).
           --out                : Output prefix (optional) [default: my-ngs].
           --output_dir         : (optional) [results will be saved to parent of input directory]
           --scaller            : Single sample variant caller; gatk, deepvariant [default: gatk]
                                  For structural variant calling, use 'svarcall'
           --jcaller            : Joint sample variant caller; gatk, glnexus [default: gatk]
           --wgs                : If data is ehole-genome sequence (it runs whole exome - wes - by default)
           --threads            : number of computer cpus to use  [default: 11].
           --njobs              : (optional) number of jobs to submit at once [default: 10]
           --help               : print this help message.
   """
}

function svarcallusage() {
   echo -e "\nUsage: genemapngs svarcall <profile> [options] ..."
   echo """
           options:
           --------

           --alignment_dir      : (required) Path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).
           --out                : Output prefix (optional) [default: my-ngs].
           --output_dir         : (optional) [results will be saved to parent of input directory]
           --scaller            : Single sample variant caller; gatk, deepvariant, dysgu, manta [default: gatk]
           --wgs                : If data is whole-genome sequence (it runs whole exome - wes - by default)
           --threads            : number of computer cpus to use  [default: 11].
           --njobs              : (optional) number of jobs to submit at once [default: 10]
           --help               : print this help message.
   """
}

function jvarcallusage() {
   echo -e "\nUsage: genemapngs jvarcall <profile> [options] ..."
   echo """
           options:
           --------

           --gvcf_dir                  : (required if importing gVCFS to genomicsdb for the first time). 
                                         Path to directory containing gVCF files and their indexes ('.tbi').
           --update                    : Specify this flag if importing gVCF files to existing genomicsdbs
           --genomicsdb_workspace_dir  : (required if importing gVCFS to an existing genomicsdb workspace). 
                                         Path to directory containing genomicsdb workspaces (workspaces must be directories).
           --batch_size                : (optional) number of samples to read into memory by GATK sample reader per time [default: 50]
           --out                       : Output prefix (optional) [default: my-ngs].
           --output_dir                : (optional) [results will be saved to parent of input directory]
           --jcaller                   : Joint sample variant caller; gatk, glnexus [default: gatk]
           --interval                  : (optional) list containing genomic intervals to process. 
                                         one chromosome name per line and/or coordinate in bed format: <chr> <start> <stop>
           --wgs                       : Specify this flag if your data is whole-genome sequence (it runs whole exome - wes - by default)
           --threads                   : number of computer cpus to use  [default: 11].
           --njobs                     : (optional) number of jobs to submit at once [default: 10]
           --help                      : print this help message.
   """
}

function varfilterusage() {
   echo -e "\nUsage: genemapngs varfilter <profile> [options] ..."
   echo """
           options:
           --------

           --vcf_dir            : (required) Path to VCF file(s).
           --minDP              : Minimum allele depth [default: 10].
           --minGQ              : Minimun genotype quality [default: 20].
           --minAC              : Minimun allele count (to remove singletons, set to 2) [default: 1].
           --out                : Output prefix (optional) [default: my-varfilter].
           --output_dir         : (optional) [results will be saved to parent of input directory]
           --jcaller            : The tool used to perform joint variant calling; gatk, glnexus [default: gatk]
           --threads            : number of computer cpus to use  [default: 4].
           --njobs              : (optional) number of jobs to submit at once [default: 10]
           --help               : print this help message.
   """
}


############################################# CONFIGURATION FILES ####################################################
function testconfig() {
#check and remove test config file if it exists
[ -e test.config ] && rm test.config

# $indir $bpm $csv $cluster $fasta $bam $out $outdir $thrds
echo """includeConfig \"\${projectDir}/nextflow.config\"
includeConfig \"\${projectDir}/configs/profile-selector.config\"
includeConfig \"\${projectDir}/configs/test.config\"
""" >> test.config
}


function qcconfig() { #params passed as arguments
#check and remove config file if it exists
#create a unique id from date and time
id=$(date +%Y%m%d%H%M%S)
[ -e ${id}-qc.config ] && rm ${id}-qc.config

#qcconfig $input_ftype $input_dir $output_dir $threads $njobs
echo """
params {
  //genemapngs qc workflow parameters
  input_ftype = '$1'                            // required: FASTQ, BAM, CRAM     (input file type)
  input_dir = '$2'                              // (required) Path to FASTQ/BAM/CRAM files.
  output_dir = '$3'                             // optional (defaults to parent of input directory) ['input_dir/../']
  threads = ${4}				// number of computer cpus to use  [default: 4]
  njobs = ${5}                                  // (optional) number of jobs to submit at once [default: 10]
}

`setglobalparams`
""" >> ${id}-qc.config
}


function trimconfig() { #params passed as arguments
#check and remove config file if it exists
#create a unique id from date and time
id=$(date +%Y%m%d%H%M%S)
[ -e ${id}-trim.config ] && rm ${id}-trim.config

#trimconfig $ftype $input_dir $output_dir $trimmer $adapter $min_length $headcrop $crop $threads $njobs
echo """
params {
  //genemapngs trim workflow parameters
  input_ftype = '$1'                            // required: FASTQ, BAM, CRAM     (input file type)
  input_dir = '$2'                              // (required) Path to FASTQ/BAM/CRAM files.
  output_dir = '$3'                             // optional (defaults to parent of input directory) ['input_dir/../']
  trimmer = '$4'                                // options: trimmomatic, trimgalore [default: trimgalore]
  adapter = '$5'                                // required if 'trimmomatic' selected [default: NP].
  min_length = ${6}                             // minimum read leangth to keep [default: 36].
  headcrop = ${7}                               // number of bases to remove from the start of reads [default: 5].
  crop = ${8}                                   // number of bases to remove from the end of reads [default: 5].
  threads = ${9}				// number of computer cpus to use  [default: 8]
  njobs = ${10}                                 // (optional) number of jobs to submit at once [default: 10]
}

`setglobalparams`
""" >> ${id}-trim.config
}


function alignconfig() { #params passed as arguments
#check and remove config file if it exists
id=$(date +%Y%m%d%H%M%S)
[ -e ${id}-alignment.config ] && rm ${id}-alignment.config

#alignconfig $pe $exome $aligner $ftype $input_dir $output_dir $dup_marker $remove_dup $spark $threads $njobs
echo """
params {
  //genemapngs align workflow parameters
  pe = $1                                       // optional: true, false [dfault: true]   (Whether reads are paired-end or single end)
  exome = $2                                    // for manta structural variant calling, specify whether WES or WGS
  aligner = '$3'                                // optional: BWA, DRAGMAP [default: BWA]
  input_ftype = '$4'                            // required: FASTQ, BAM, CRAM     (input file type)
  input_dir = '$5'                              // (required) Path to FASTQ/BAM/CRAM files.
  output_dir = '$6'                             // optional (defaults to parent of input directory) ['input_dir/../']
  dup_marker = '$7'                             // Duplicate marker tool (optional); sambamba, samtools [default: sambamba]
  remove_dup = '$8'                             // Whether to remove duplicates (optional): true, false [default: false]
  spark = $9                                    // (optional) use GATK sprak mode for multi-threaded post-alignment processing; true, false [default: false]
  threads = ${10}				// number of computer cpus to use  [default: 11]
  njobs = ${11}                                 // (optional) number of jobs to submit at once [default: 10]
}

`setglobalparams`
""" >> ${id}-alignment.config
}


function varcallconfig() { #params passed as arguments
#check and remove config file if it exists
[ -e ${4}-varcall.config ] && rm ${4}-varcall.config

#varcallconfig $exome $input_dir $output_dir $output_prefix $scaller $jcaller $threads $njobs
echo """
params {
  //=======================================
  // genemapngs varcall workflow parameters 
  //=======================================
  exome = $1                                    // for manta structural variant calling, specify whether WES or WGS
  alignment_dir = '$2'                          // (required) Path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).
  output_dir = '$3'                             // optional (defaults to parent of input directory) ['input_dir/../']
  output_prefix = '$4'                          // required
  single_caller = '$5'                          // options: gatk, deepvariant
  joint_caller = '$6'                           // options: gatk, glnexus
  threads = ${7}                                // number of computer cpus to use  [default: 11]
  njobs = ${8}                                  // (optional) number of jobs to submit at once [default: 10]
}

`setglobalparams`
""" >> ${4}-varcall.config
}


function svarcallconfig() { #params passed as arguments
#check and remove config file if it exists
[ -e ${4}-svarcall.config ] && rm ${4}-svarcall.config

#svarcallconfig $exome $alignment_dir $output_dir $output_prefix $scaller $threads $njobs
echo """
params {
  //=======================================
  //genemapngs svarcall workflow parameters
  //=======================================
  exome = $1                                    // for manta structural variant calling, specify whether WES or WGS
  alignment_dir = '$2'                          // (required) Path to alignment (BAM/CRAM) files and their indexes (.bai/.crai)..
  output_dir = '$3'                             // optional (defaults to parent of input directory) ['input_dir/../']
  output_prefix = '$4'                          // required
  single_caller = '$5'                          // options: gatk, deepvariant
  threads = ${6}                                // number of computer cpus to use  [default: 11]
  njobs = ${7}                                  // (optional) number of jobs to submit at once [default: 10]
}

`setglobalparams`
""" >> ${4}-svarcall.config
}


function jvarcallconfig() { #params passed as arguments
#check and remove config file if it exists
[ -e ${6}-jvarcall.config ] && rm ${6}-jvarcall.config

#jvarcallconfig $exome $gvcf_dir $update $genomicsdb_workspace_dir $output_dir $output_prefix $jcaller $interval $threads $njobs ${batch_size}
echo """
params {
  //=======================================
  //genemapngs jvarcall workflow parameters
  //=======================================
  exome = $1                                    
  gvcf_dir = '$2'                               
  update = ${3}                                 
  genomicsdb_workspace_dir = '$4'               
  batch_size = ${11}                            
  output_dir = '$5'                             
  output_prefix = '$6'                          
  joint_caller = '$7'                           
  interval = '$8'                               
  threads = ${9}                                
  njobs = ${10}                                 


  /*****************************************************************************************
  -exome:
     for manta structural variant calling, specify whether WES or WGS
  -gvcf_dir:
     required if importing gVCFS to genomicsdb for the first time) 
     Path to directory containing gVCF files and their indexes ('.tbi').
  -update:
     whether to add gVCFs to existing genomicsdb workspaces. 
     If true, 'genomicsdb_workspace_dir' must be provided [defaul: false]
  -genomicsdb_workspace_dir: 
     (required if importing gVCFS to an existing genomicsdb workspace) 
     Path to directory containing genomicsdb workspaces (workspaces must be directories).
  -batch_size:
     (optional) number of samples to read into memory by GATK 
     sample reader per time [default: 50]
  -output_dir:
     optional (defaults to parent of input directory) ['gvcf_dir/../']
  -output_prefix:
     (required) name to add to output files.
  -joint_caller:
     options: gatk, deepvariant [default: gatk]
  -interval:
     (optional) list containing genomic intervals to process. 
     One chromosome name per line and/or coordinate in bed format: <chr> <start> <stop>
     If not provided, intervals will be creared from CRAM/gVCF header.
  -threads:
    (optional) number of computer cpus to use  [default: 11]
  -njobs:
      (optional) number of jobs to submit at once [default: 10]
  *******************************************************************************************/
}

`setglobalparams`
""" >> ${6}-jvarcall.config
}


function varfilterconfig() {
#check and remove config file if it exists
[ -e ${5}-varfilter.config ] && rm ${5}-varfilter.config

#varfilterconfig $vcf_dir $minDP $minGQ $minAC $out $output_dir $threads $njobs
echo """
params {
  //genemapngs varfilter workflow parameters
  vcf_dir = '${1}'                              // (required) Path to VCF file(s).
  minDP = ${2}                                  // Minimum allele depth [default: 10].
  minGQ = ${3}                                  // Minimun genotype quality [default: 20].
  minAC = ${4}                                  // Minimun allele count (to remove singletons, set to 2) [default: 1].
  output_prefix = '${5}'                        // Output prefix (optional) [default: my-varfilter].
  output_dir = '${6}'                           // (optional) [results will be saved to parent of input directory]
  joint_caller = '${7}'                         // options: gatk, glnexus
  threads = ${8}                                // number of computer cpus to use  [default: 4].
  njobs = ${9}                                  // (optional) number of jobs to submit at once [default: 10]
}

`setglobalparams`
""" >> ${5}-varfilter.config
}



if [ $# -lt 1 ]; then
   usage; exit 1;
else
   case $1 in
      test)
         profile='local,singularity,hg19'
         testconfig
      ;;
      qc)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            qcusage;
            exit 1;
         fi

         prog=`getopt -a --long "help,ftype:,input_dir:,output_dir:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #- defaults
         ftype=FASTQ                             #// optional: FASTQ, BAM, CRAM     (input file type)
         input_dir=NULL                          #// (required) Path to FASTQ/BAM/CRAM files.
         output_dir=NULL                         #// optional (defaults to parent of input directory) ['input_dir/../']
         threads=4
         njobs=10                                #// (optional) number of jobs to submit at once [default: 10]

         eval set -- "$prog"

         while true; do
            case $1 in
               --ftype) ftype="$2"; shift 2;;
               --input_dir) input_dir="$2"; shift 2;;
               --output_dir) output_dir="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; 2> /dev/null; shift 2;;
               --help) shift; qcusage; 2> /dev/null; exit 1;;
               --) shift; 2> /dev/null; break;;
               *) shift; qcusage; 2> /dev/null; exit 1;;
            esac
            continue; shift;
         done

         #- check required options
         #check_ftype $ftype
         #check_common_required_params $input_dir $output_dir $threads $njobs
         check_required_params input_dir,$input_dir
         check_output_dir $output_dir
         check_optional_params ftype,$ftype threads,$threads njobs,$njobs

         qcconfig $ftype $input_dir $output_dir $threads $njobs
         #echo `nextflow -c ${out}-qc.config run qualitycontrol.nf -profile $profile -w ${outdir}/work/`
      ;;
      trim)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            trimusage;
            exit 1;
         fi

         prog=`getopt -a --long "help,ftype:,input_dir:,output_dir:,trimmer:,adapter:,min_length:,headcrop:,crop:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #- defaults
         ftype=FASTQ                             #// optional: FASTQ, BAM, CRAM     (input file type)
         input_dir=NULL                          #// (required) Path to FASTQ/BAM/CRAM files.
         output_dir=NULL                         #// optional (defaults to parent of input directory) ['input_dir/../']
         trimmer=trimgalore                      #// options: trimmomatic, trimgalore [default: trimgalore]
         adapter=NP                              #// required if 'trimmomatic' selected [default: NP].
         min_length=36                           #// minimum read leangth to keep.
         headcrop=5                              #// number of bases to remove from the start of reads
         crop=5                                  #// number of bases to remove from the end of reads
         threads=8
         njobs=10                                #// (optional) number of jobs to submit at once [default: 10]

         eval set -- "$prog"

         while true; do
            case $1 in
               --ftype) ftype="$2"; shift 2;;
               --input_dir) input_dir="$2"; shift 2;;
               --output_dir) output_dir="$2"; shift 2;;
               --trimmer) trimmer="$2"; shift 2;;
               --adapter) adapter="$2"; shift 2;;
               --min_length) min_length="$2"; shift 2;;
               --headcrop) headcrop="$2"; shift 2;;
               --crop) crop="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; trimusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; trimusage; 1>&2; exit 1;;
            esac
            continue; shift;
         done

         #- check required options
         #check_ftype $ftype
         #check_common_required_params $input_dir $output_dir $threads $njobs
         check_required_params input_dir,$input_dir $([[ "${trimmer}" == "trimmomatic" ]] && echo "\$adapter,$adapter")
         check_output_dir $output_dir
         check_optional_params ftype,$ftype trimmer,$trimmer min_length,$min_length headcrop,$headcrop crop,$crop threads,$threads njpbs,$njobs

         trimconfig $ftype $input_dir $output_dir $trimmer $adapter $min_length $headcrop $crop $threads $njobs
         #echo `nextflow -c ${out}-qc.config run qualitycontrol.nf -profile $profile -w ${outdir}/work/`
      ;;
      align)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            alignusage; 1>&2;
            exit 1;
         fi

         prog=`getopt -a --long "help,se,wgs,spark,aligner:,ftype:,input_dir:,output_dir:,dup_marker:,remove_dup:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #defaults
         pe=true                        #// optional: true, false [dfault: ture]   (Whether reads are paired-end or single end)
         aligner=BWA                    #// optional: BWA, DRAGMAP [default: BWA]
         ftype=FASTQ                    #// required: FASTQ, BAM, CRAM     (input file type)
         input_dir=NULL                 #// required
         output_dir=NULL                #// optional [default: ${input_dir}/../]
         exome=true                     #// for manta structural variant calling, specify whether WES or WGS
         dup_marker=sambamba            #// optional: samtools, sambamba [default: sambamba]
         remove_dup=false               #// optional: true, false [default: false]
         spark=false
         threads=11
         njobs=10                       #// (optional) number of jobs to submit at once [default: 10]
         
         eval set -- "$prog"

         while true; do
            case $1 in
               --se) pe=false; shift;;
               --wgs) exome=false; shift;;
               --spark) spark=true; shift;;
               --aligner) aligner="${2}"; shift 2;;
               --ftype) ftype="$2"; shift 2;;
               --input_dir) input_dir="$2"; shift 2;;
               --output_dir) output_dir="${2}"; shift 2;;
               --dup_marker) dup_marker="$2"; shift 2;;
               --remove_dup) remove_dup="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; alignusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; alignusage; 1>&2; exit 1;;
            esac
         done

         #test arguments values
         #check_ftype $ftype
         #check_common_required_params $input_dir $output_dir $threads $njobs
         check_required_params input_dir,$input_dir
         check_output_dir $output_dir
         check_optional_params ftype,$ftype se,$pe wgs,$exome dup_marker,$dup_marker remove_dup,$remove_dup spark,$spark threads,$threads njobs,$njobs

#         if [[ $aligner == -* ]]; then
#            echo "ERROR: Invalid paramter value for option '--aligner'"
#            1&>2; exit 1;
#         fi
#         if [[ $output_prefix == -* ]]; then
#            echo "ERROR: Invalid paramter value for option '--output_prefix'"
#            1&>2; exit 1;
#         fi
#         if [[ $scaller == -* ]]; then
#            echo "ERROR: Invalid paramter value for option '--scaller'"
#            1>&2; exit 1;
#         fi
         #if [[ $gvcf == -* ]] || [[ $gvcf == NULL ]]; then
         #   echo "ERROR: Invalid paramter value for option '--gvcf_dir'"
         #   1>&2; exit 1;
         #fi

         #args=($pe $exome $spark $aligner $ftype $input_dir $output_dir $output_prefix $scaller $jcaller $threads $njobs)
         #echo ${#args[@]}
         alignconfig $pe $exome $aligner $ftype $input_dir $output_dir $dup_marker $remove_dup $spark $threads $njobs

         #echo `nextflow -c ${out}-idat2vcf.config run idat2vcf.nf -profile $profile -w ${outdir}/work/`

      ;;
      varcall)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            varcallusage; 1>&2;
            exit 1;
         fi

         prog=`getopt -a --long "help,wgs,alignment_dir:,output_dir:,out:,scaller:,jcaller:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #defaults
         #ftype=FASTQ                    #// required: FASTQ, BAM, CRAM     (input file type)
         alignment_dir=NULL             #// (optional) path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).
         output_dir=NULL                #// optional [default: ${input_dir}/../]
         output_prefix="my-ngs"         #// required [default: 'my-ngs']
         ped=NULL                       #// optional [disabled]
         scaller=gatk                   #// options: gatk, deepvariant, dysgu, manta
         exome=true                     #// for manta structural variant calling, specify whether WES or WGS
         jcaller=gatk                   #// options: gatk, glnexus
         threads=11
         njobs=10                       #// (optional) number of jobs to submit at once [default: 10]

         eval set -- "$prog"

         while true; do
            case $1 in
               --wgs) exome=false; shift;;
               #--ftype) ftype="$2"; shift 2;;
               --alignment_dir) alignment_dir="$2"; shift 2;;
               --output_dir) output_dir="${2}"; shift 2;;
               --out) output_prefix="$2"; shift 2;;
               --scaller) scaller="$2"; shift 2;;
               --jcaller) jcaller="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; varcallusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; varcallusage; 1>&2; exit 1;;
            esac
         done

         #test arguments values
         #check_common_required_params $input_dir $output_dir $threads $njobs
         check_required_params alignment_dir,$alignment_dir
         check_output_dir $output_dir
         check_optional_params output_prefix,$output_prefix scaller,$scaller jcaller,$jcaller wgs,$exome threads,$threads njobs,$njobs

#         if [[ $output_prefix == -* ]]; then
#            echo "ERROR: Invalid paramter value for option '--output_prefix'"
#            1&>2; exit 1;
#         fi
#         if [[ $scaller == -* ]]; then
#            echo "ERROR: Invalid paramter value for option '--scaller'"
#            1>&2; exit 1;
#         fi
#         if [[ $jcaller == -* ]]; then
#            echo "ERROR: Invalid paramter value for option '--jcaller'"
#            1>&2; exit 1;
#         fi

         #args=($pe $exome $spark $aligner $ftype $input_dir $output_dir $output_prefix $scaller $jcaller $threads $njobs)
         #echo ${#args[@]}
         varcallconfig $exome $alignment_dir $output_dir $output_prefix $scaller $jcaller $threads $njobs

         #echo `nextflow -c ${out}-idat2vcf.config run idat2vcf.nf -profile $profile -w ${outdir}/work/`

      ;;
      svarcall)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            svarcallusage; 1>&2;
            exit 1;
         fi

         prog=`getopt -a --long "help,wgs,alignment_dir:,output_dir:,out:,scaller:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #defaults
         #ftype=FASTQ                    #// required: FASTQ, BAM, CRAM     (input file type)
         alignment_dir=NULL                 #// (optional) path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).
         output_dir=NULL                #// optional [default: ${input_dir}/../]
         output_prefix="my-ngs"         #// required [default: 'my-ngs']
         ped=NULL                       #// optional [disabled]
         scaller=gatk                   #// options: gatk, deepvariant, dysgu, manta
         exome=true                     #// for manta structural variant calling, specify whether WES or WGS
         threads=11
         njobs=10                       #// (optional) number of jobs to submit at once [default: 10]

         eval set -- "$prog"

         while true; do
            case $1 in
               --wgs) exome=false; shift;;
               #--ftype) ftype="$2"; shift 2;;
               --alignment_dir) alignment_dir="$2"; shift 2;;
               --output_dir) output_dir="${2}"; shift 2;;
               --out) output_prefix="$2"; shift 2;;
               --scaller) scaller="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; svarcallusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; svarcallusage; 1>&2; exit 1;;
            esac
         done

         #test arguments values
         #check_common_required_params $input_dir $output_dir $threads $njobs
         check_required_params alignment_dir,$alignment_dir
         check_output_dir $output_dir
         check_optional_params output_prefix,$output_prefix scaller,$scaller wgs,$exome threads,$threads njobs,$njobs

         svarcallconfig $exome $alignment_dir $output_dir $output_prefix $scaller $threads $njobs

         #echo `nextflow -c ${out}-idat2vcf.config run idat2vcf.nf -profile $profile -w ${outdir}/work/`

      ;;
      jvarcall)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            jvarcallusage; 1>&2;
            exit 1;
         fi

         prog=`getopt -a --long "help,wgs,gvcf_dir:,update,genomicsdb_workspace_dir:,batch_size:,output_dir:,out:,jcaller:,interval:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #defaults
         #ftype=FASTQ                    #// required: FASTQ, BAM, CRAM     (input file type)
         output_dir=NULL                #// optional [default: ${input_dir}/../]
         output_prefix="my-ngs"         #// required [default: 'my-ngs']
         ped=NULL                       #// optional [disabled]
         exome=true                     #// for manta structural variant calling, specify whether WES or WGS
         jcaller=gatk                   #// options: gatk, glnexus
         interval=NULL                  #// optional: list containing interval to process
         gvcf_dir=NULL                  #// (required if creating new genomicsdb workspaces) path to pre-existing gVCF files
         update=false                   #// whether to add gVCF files to existing genomicsdb workspaces
         genomicsdb_workspace_dir=NULL  #// (required if updating exisiting enomicsdb workspaces)
         batch_size=50                  #// (optional) number of samples to read into memory by GATK sample reader per time [default: 50]
         threads=11
         njobs=10                       #// (optional) number of jobs to submit at once [default: 10]

         eval set -- "$prog"

         while true; do
            case $1 in
               --wgs) exome=false; shift;;
               --update) update=true; shift;;
               #--ftype) ftype="$2"; shift 2;;
               --output_dir) output_dir="${2}"; shift 2;;
               --out) output_prefix="$2"; shift 2;;
               --jcaller) jcaller="$2"; shift 2;;
               --interval) interval="$2"; shift 2;;
               --gvcf_dir) gvcf_dir="$2"; shift 2;;
               --genomicsdb_workspace_dir) genomicsdb_workspace_dir="$2"; shift 2;;
               --batch_size) batch_size="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; jvarcallusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; jvarcallusage; 1>&2; exit 1;;
            esac
         done

         #test arguments values
         #check_common_required_params $input_dir $output_dir $threads $njobs
         if [[ "${update}" == "true" ]]; then
            check_required_params gvcf_dir,$gvcf_dir genomicsdb_workspace_dir,$genomicsdb_workspace_dir
         else
            check_required_params gvcf_dir,$gvcf_dir
         fi
         check_output_dir $output_dir
         check_optional_params \
             output_prefix,$output_prefix \
             jcaller,$jcaller \
             interval,$interval \
             wgs,$exome \
             genomicsdb_workspace_dir,$genomicsdb_workspace_dir \
             threads,$threads \
             njobs,$njobs \
             batch_size,$batch_size

         jvarcallconfig \
             $exome \
             $gvcf_dir \
             $update \
             $genomicsdb_workspace_dir \
             $output_dir \
             $output_prefix \
             $jcaller \
             $interval \
             $threads \
             $njobs \
             $batch_size

         #echo `nextflow -c ${out}-idat2vcf.config run idat2vcf.nf -profile $profile -w ${outdir}/work/`

      ;;
      varfilter)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            varfilterusage;
            exit 1;
         fi

         prog=`getopt -a --long "help,vcf_dir:,minDP:,minGQ:,minAC:,out:,output_dir:,jcaller:,threads:,njobs:" -n "${0##*/}" -- "$@"`;

         #- defaults
         vcf_dir=NULL                            #// (required) Path to FASTQ/BAM/CRAM files.
         minDP=10                                #// Minimum allele depth [default: 10].
         minGQ=20                                #// Minimun genotype quality [default: 20].
         minAC=1                                 #// Minimun allele count (to remove singletons, set to 2) [default: 1].
         out="my-varfilter"                      #// Output prefix (optional) [default: my-ngs].
         output_dir=NULL                         #// (optional) [results will be saved to parent of input directory]
         jcaller=gatk                            #// options: gatk, glnexus
         threads=4                               #// number of computer cpus to use  [default: 11].
         njobs=10                                #// (optional) number of jobs to submit at once [default: 11]

         eval set -- "$prog"

         while true; do
            case $1 in
               --vcf_dir) vcf_dir="$2"; shift 2;;
               --minDP) minDP="$2"; shift 2;;
               --minGQ) minGQ="$2"; shift 2;;
               --minAC) minAC="$2"; shift 2;;
               --out) out="$2"; shift 2;;
               --output_dir) output_dir="$2"; shift 2;;
               --jcaller) jcaller="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; varfilterusage; exit 1;;
               --) shift; break;;
               *) shift; varfilterusage; exit 1;;
            esac
            continue; shift;
         done

         #- check required options
         if [[ $vcf_dir == -* ]]  || [[ $vcf_dir == NULL ]]; then
            echo "ERROR: Invalid paramter value for option '--vcf_dir'"
            exit 1;
         fi
         if [[ $output_dir == -* ]]  || [[ $output_dir == NULL ]]; then
            output_dir="${vcf_dir}/../"
            if [[ $output_dir == NULL* ]]; then
               echo "ERROR: Invalid paramter value for option '--output_dir'"
               exit 1;
            fi
         fi
         if [[ $threads == -* ]]; then
            echo "ERROR: Invalid paramter value for option '--threads'"
            exit 1;
         fi
         if [[ $njobs == -* ]]; then
            echo "ERROR: Invalid paramter value for option '--njobs'"
            exit 2;
         fi

         varfilterconfig $vcf_dir $minDP $minGQ $minAC $out $output_dir $jcaller $threads $njobs
         #echo `nextflow -c ${out}-qc.config run qualitycontrol.nf -profile $profile -w ${outdir}/work/`
      ;;

      assoc) echo "assoc"; shift ;;
      help) shift; varfilterusage; exit 1;;
      *) shift; usage; exit 1;;
   esac

   #echo -e "\nRunning ${comd}...\n"


fi


