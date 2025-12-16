import os
import sys
import importlib
import argparse
import pathlib
import shutil
import json
import psutil
import random
import string
import time

def test_fxn1(name=None):
    print(f"NAME: {name}")

def test_fxn2(dob=None):
    print(f"FOUNDED: {dob}")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CHECK REQUIRED PACKAGES
def is_bash_tool_installed(tool_name):
    """
    Checks if a given Bash tool is installed and accessible in the system's PATH.
    """
    return shutil.which(tool_name) is not None

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CHECK PYTHON PACKAGE INSTALLATION
# A dictionary to hold the imported modules
imported_modules = {}

def check_python_package(module_list=None):
    for module_name in module_list:
        try:
            # Import the module using importlib.import_module()
            module = importlib.import_module(module_name)
            # Store the imported module in the dictionary
            imported_modules[module_name] = module
            #print(f"Successfully imported {module_name}")
        except ImportError as e:
            print(f"Error importing {module_name}: {e}")    


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE UNIQUE UD FOR EACH RUN
def get_run_id(length):
    # Define the possible characters: only letters and digits in this example
    characters = string.ascii_letters + string.digits
    # Use random.choices() which is efficient for generating multiple selections with replacement
    random_string = ''.join(random.choices(characters, k=length))
    return random_string


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MONITOR RESOURCES
def monitor_resources():
    while True:
        # CPU Usage
        cpu_percent = psutil.cpu_percent(interval=1)  # Interval for sampling CPU usage
        print(f"CPU Usage: {cpu_percent}%")

        # Memory Usage
        memory = psutil.virtual_memory()
        print(f"Memory Usage: {memory.percent}% (Used: {memory.used / (1024**3):.2f} GB / Total: {memory.total / (1024**3):.2f} GB)")

        # Disk Usage (for the root partition)
        disk = psutil.disk_usage('/')
        print(f"Disk Usage: {disk.percent}% (Used: {disk.used / (1024**3):.2f} GB / Total: {disk.total / (1024**3):.2f} GB)")

        # Network Usage (example for total bytes sent/received)
        net_io = psutil.net_io_counters()
        print(f"Network (Total): Sent: {net_io.bytes_sent / (1024**2):.2f} MB, Received: {net_io.bytes_recv / (1024**2):.2f} MB")

        time.sleep(5)  # Wait for 5 seconds before the next check

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USAGE
def get_usage(
    script_name=None,
    version=None     
):
#    wrapper_script = sys.argv[0]
#    script_name = os.path.basename(wrapper_script)
#    script_path = os.path.dirname(wrapper_script)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # colors
    ANSIRESET = '\033[0m'
    ANSIRED = '\033[31m'
    ANSIGRN = '\033[32m'
    ANSIYLW = '\033[33m'
    ANSIBLU = '\033[34m'
    ANSILBL = '\033[36m'
    ANSIPPL = '\033[35m'
    ANSIGRY = '\033[2m'

    grey_pipe = f"{ANSIGRY}|{ANSIRESET}"

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # define banner and usage messages
    descmsg = f"""
        {ANSIGRY}#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#{ANSIRESET}
        {grey_pipe}  {ANSIYLW}WRAPPER FOR THE {ANSILBL}ESOH{ANSIRED}INFORMATICS {ANSIPPL}NGS {ANSIGRN}NEXTFLOW {ANSIYLW}WORKFLOW{ANSIRESET}  {grey_pipe}
        {grey_pipe}                    Kevin Esoh, 2025                     {grey_pipe}
        {grey_pipe}                   kesohku1@jhmi.edu                     {grey_pipe}
        {ANSIGRY}#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#{ANSIRESET}
    """

#    version="0.1 (beta)"

    usage = f"""
    {descmsg}

    LICENSCE: GPL-3.0 license
    VERSION: {version}

    Usage: {script_name} <command> [-v/--version] [-h/--help] <options>

    Commands:
            test: Run test to see if workfow installed correctly.
                qc: Check FASTQ or Alignment (BAM/CRAM) quality.
            trim: Trim adapters and poor quality bases from reads (input is FASTQ or BAM/CRAM).
            align: Align/map reads to reference and post-alignment processing (input is FASTQ or BAM/CRAM).
        mergealign: Megre Alignment (BAM/CRAM) files.
            call: Perform variant calling (both single and joint sample) in one run.
            scall: Perform only sinlge sample variant calling to generate gVCF files.
            jcall: Perform only joint (multi-sample) variant calling with pre-existing gVCF files.
            filter: Filter variant calls in VCF/BCF files.
        annotate: Annotate variants with ANNOVAR
            raw2vcf: Run entire workflow from 'align' to 'call' (i.e., 'scall' + 'jcall')
    raw2annotate: Run entire workflow from 'align' to 'annotate'

        [NOTE] For 'trim' and 'align', BAM/CRAM input is first converted to FASTQ

    """
    return usage

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DEFINE ARGUMENTS
def get_arguments(descmsg=None, prog=None, version=None):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # GROUP ARGUMENTS SHARED BY ALL COMMANDS
    parent_parser = argparse.ArgumentParser(
        description=descmsg,
        formatter_class=argparse.RawTextHelpFormatter,
        #formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False
    )

    required_options = parent_parser.add_argument_group('Required')
    optional_arguments = parent_parser.add_argument_group('Optional')

    required_options.add_argument(
        "--output_dir",
        help="Where to store results. This is required because output files from NGS analyses can be very large.",
        required=True,
        type=pathlib.Path,
        metavar="<path>"
    )

    optional_arguments.add_argument(
        "-h", "--help",
        action="help",
        help="Show help message and exit."
    )

    optional_arguments.add_argument(
        "-v", "--version",
        action="version",
        version=f"{version}"
    )

    optional_arguments.add_argument(
        '--wgs',
        help="Specify this flag if your data is whole-genome sequencing (WGS) (whole-exome - WES - is assumed by default)",
        action="store_true"
    )

    optional_arguments.add_argument(
        "--profile",
        #choices=[
        #    "local", "singularity", "docker", "apptainer", 
        #    "singularity,slurm"
        #],
        help="Select a profile to execute the commands with [defaul: singularity]",
        required=False,
        default="singularity",
        type=str,
        metavar="<text>"
    )

    optional_arguments.add_argument(
        "--threads",
        help="Number of computer cpus to use [default: 1]",
        required=False,
        default=1,
        type=int,
        metavar="<integer>"        
    )

    optional_arguments.add_argument(
        "--njobs",
        help="Number of jobs to run simultaneously [default: 1]",
        required=False,
        default=1,
        type=int,
        metavar="<integer>"
    )

    optional_arguments.add_argument(
        "--resume",
        help="""
        Use this flag to resume the nextflow execution. Either the flag without any value (--resume) to resume from last
        exit point or add the flag with the specific RUN NAME as can be obtained using `nexflow log` to resume from
        a specific point (--resume kickass_church).
        """,
        nargs='?',
        const=False,
        default=False,
        required=False,
        metavar="uuid"
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ALL ARGUMENTS SHARED BY QC, TRIM AND ALIGN COMMANDS GO HERE
    # THESE ALSO INHERIT THE COMMON ARGUMENTS ABOVE
    qc_trim_align_common_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[parent_parser],
        formatter_class=argparse.RawTextHelpFormatter,
        #formatter_class=argparse.RawDescriptionHelpFormatter
  )

    qc_trim_align_common_required = qc_trim_align_common_parser.add_argument_group("Required")
    qc_trim_align_common_optional = qc_trim_align_common_parser.add_argument_group("Optional")

    qc_trim_align_common_required.add_argument(
        "--input_dir",
        help="Path to FASTQ/BAM/CRAM input files.",
        required=True,
        type=pathlib.Path,
        metavar="<path>"
    )

    qc_trim_align_common_optional.add_argument(
        "--ftype",
        help="Input file type; 'FASTQ', 'BAM', 'CRAM', 'VCF' [default: FASTQ]",
        required=False,
        default="FASTQ",
        type=str,
        metavar="<text>"
    )


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # REFERENCE GENOME PARSER
    ref_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[parent_parser],
        formatter_class=argparse.RawTextHelpFormatter
    )

    ref_parser_required = ref_parser.add_argument_group("Required")

    ref_parser_required.add_argument(
        "--build",
        help="Reference build to use [default: %(default)s]",
        choices=["hg19", "hg38", "t2t"],
        default="hg19",
        type=str,
        required=True
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FILTER and ANNOTATE SHARED ARGUMENTS
    filter_annotate_common_parser = argparse.ArgumentParser(
        add_help=False,
        #parents=[parent_parser],
        parents=[ref_parser],
        formatter_class=argparse.RawTextHelpFormatter,
        #formatter_class=argparse.RawDescriptionHelpFormatter
    )

    filter_annotate_common_required = filter_annotate_common_parser.add_argument_group("Required")
    filter_annotate_common_optional = filter_annotate_common_parser.add_argument_group("Optional")

    filter_annotate_common_required.add_argument(
        "--vcf_dir",
        help="""
        Path containing VCF file(s).
        """,
        required=True,
        type=pathlib.Path,
        metavar="<path>"
    )

    filter_annotate_common_optional.add_argument(
        "--left_norm",
        help="""
        Add this flag to left-normalized variants such as is recommended by ANNOVAR.
        """,
        required=False,
        action="store_true"
    )

    filter_annotate_common_optional.add_argument(
        "--out",
        help="""
        Output prefix [default: my-ngs-vcf].
        """,
        required=False,
        default='my-ngs-vcf',
        type=str,
        metavar="<text>"
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # TRIM-SPECIFIC OPTIONS
    trim_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[qc_trim_align_common_parser],
        #formatter_class=argparse.RawDescriptionHelpFormatter
    )

    trim_required = trim_parser.add_argument_group("Required")
    trim_optional = trim_parser.add_argument_group("Optional")

    trim_optional.add_argument(
        "--trimmer",
        help="""
        Tool to use for trimming.
        options: trimmomatic, trimgalore [default: trimgalore]
        """,
        required=False,
        metavar="<text>"
    )

    trim_optional.add_argument(
        "--adapter",
        help="""
        Adapter sequence to use to clip adapters from reads [default: NP].

        Note: 
            - This is only required if 'trimmomatic' trimmer is selected.
            - 'trimgalore' will auto-detect adapters, hence suitable to process data from different sequencing companies in one batch

            options:
                    NP   --> NexteraPE-PE.fa
                    T3U  --> TruSeq3-PE-2.fa [Illumina universal]
                    T2P  --> TruSeq2-PE.fa
                    T2S  --> TruSeq2-SE.fa
                    T3P  --> TruSeq3-PE.fa
                    T3S  --> TruSeq3-SE.fa     
        """,
        required=False,
        type=str,
        default="NP",
        metavar="<text>"
    )

    trim_optional.add_argument(
        "--min_length",
        help="Minimum read leangth to keep [default: 36]",
        required=False,
        default=36,
        type=int,
        metavar="<integer>"
    )

    trim_optional.add_argument(
        "--headcrop",
        help="Number of bases to remove from the start of reads [default: 5]",
        required=False,
        default=5,
        type=int,
        metavar="<integer>"
    )

    trim_optional.add_argument(
        "--crop",
        help="Number of bases to remove from the end of reads [default: 5]",
        required=False,
        default=5,
        type=int,
        metavar="<integer>"
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ALIGN-SPECIFIC OPTIONS
    align_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[qc_trim_align_common_parser],
    )

    align_required = align_parser.add_argument_group("Required")
    align_optional = align_parser.add_argument_group("Optional")

    align_required.add_argument(
        "--build",
        help="Reference build to use [default: %(default)s]",
        choices=["hg19", "hg38", "t2t"],
        default="hg19",
        type=str,
        required=True
    )

    align_required.add_argument(
        "--aligner",
        help="Alignment tool; 'BWA', 'DRAGMAP' [default: BWA]",
        required=False,
        default="BWA",
        type=str,
        metavar="<text>"
    )

    align_optional.add_argument(
        "--se",
        help="Add this flag if your data is single-end (SE) reads. By default, paired-end (PE) reads are assumed.",
        required=False,
        action="store_true"
    )

    align_optional.add_argument(
        "--dup_marker",
        help="""
        Duplicate marker tool; 'sambamba', 'samtools' [default: sambamba].
        NOTE: New versions of sambamba do not support CRAM output. So if would like to benefit from the gain in 
        storage space that CRAM output offers, select 'samtools'.
        """,
        required=False,
        default="sambamba",
        type=str,
        metavar="<text>"
    )

    align_optional.add_argument(
        "--remove_dup",
        help="""
        Add this flag to remove duplicate reads.
        """,
        required=False,
        action="store_true"
    )

    align_optional.add_argument(
        "--spark",
        help="""
        Add this flag to use GATK multi-threaded SPARK mode for post-alignment processing.
        """,
        required=False,
        action="store_true"
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MERGE-ALIGN-SPECIFIC OPTIONS
    mergealign_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[parent_parser],
    )

    mergealign_optional = mergealign_parser.add_argument_group("Optional")

    mergealign_optional.add_argument(
        "--sort_order",
        help="Sort order of input alignment files; 'name', 'coordinate' [default: coordinate].",
        required=False,
        default='coordinate',
        type=str,
        metavar="<text>"
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # CALL-SPECIFIC OPTIONS
    call_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[parent_parser],
    )

    call_required = call_parser.add_argument_group("Required")
    call_optional = call_parser.add_argument_group("Optional")

    call_required.add_argument(
        "--alignment_dir",
        help="Path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).",
        required=True,
        type=pathlib.Path,
        metavar="<path>"
    )

    call_optional.add_argument(
        "--scaller",
        help="""
        Single sample variant caller; 'gatk', 'deepvariant' [default: gatk].
        'gatk' will use GATK HaplotypeCaller for single sample calling.
        For single sample calling of structural variants, use the 'scall' command.
        """,
        required=False,
        default='gatk',
        type=str,
        metavar="<text>"
    )

    call_optional.add_argument(
        "--jcaller",
        help="""
        Joint sample variant caller; 'gatk', 'glnexus' [default: gatk].
        """,
        required=False,
        default='gatk',
        type=str,
        metavar="<text>"
    )

    call_optional.add_argument(
        "--batch_size",
        help="""
        Number of samples to read into memory by GATK sample reader each time [default: 50].
        """,
        required=False,
        default=50,
        type=int,
        metavar="<integer>"
    )

    call_optional.add_argument(
        "--interval",
        help="""
        List containing genomic intervals, one chromosome name per line and/or coordinate in bed format: <chr> <start> <stop>.
        NB: Ensure that your chromosome names are the same as in the reference (e.g. chr1 or 1).
		If not provided, non-overlapping intervals (5M bp) will be generated from gVCF header..
        """,
        required=False,
        type=pathlib.Path,
        metavar="<file>"
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # SCALL-SPECIFIC OPTIONS
    scall_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[parent_parser],
    )

    scall_required = scall_parser.add_argument_group("Required")
    scall_optional = scall_parser.add_argument_group("Optional")

    scall_required.add_argument(
        "--alignment_dir",
        help="Path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).",
        required=True,
        type=pathlib.Path,
        metavar="<path>"
    )

    scall_required.add_argument(
        "--vcf_dir",
        help="Path to VCF files and their indexes (.tbi) for 'dysgu merge' only!. NOTE: You ust specify '--ftype VCF'",
        required=False, # will be conditionally required if --ftype is VCF
        type=pathlib.Path,
        metavar="<path>"
    )

    scall_optional.add_argument(
        "--scaller",
        help="""
        Joint sample variant caller; gatk-hap, gatk-som, gatk-mt, deepvariant, dysgu, manta [default: gatk-hap]
        NOTE: gatk-hap -> GATK haplotypeCaller, gatk-som -> GATK Mutect2 (Somatic caller), gatk-mt -> GATK Mitochondria caller.
        """,
        required=False,
        default='gatk',
        type=str,
        metavar="<text>"
    )

    scall_optional.add_argument(
        "--interval",
        help="""
        List containing genomic intervals, one chromosome name per line and/or coordinate in bed format: <chr> <start> <stop>.
        This is only relevant for 'Dysgu'. HaplotypeCaller or DeepVariant calling per interval is not very efficient thousands of 
        intervals would need to be generated per sample.
        NB: Ensure that your chromosome names are the same as in the reference (e.g. chr1 or 1).
		If not provided, intervals of 5M bp with 1kb overlaps will be generated from VCF header.
        """,
        required=False,
        type=pathlib.Path,
        metavar="<file>"
    )

    scall_optional.add_argument(
        "--out",
        help="""
        Output prefix [default: my-ngs].
        """,
        required=False,
        default='my-ngs',
        type=str,
        metavar="<text>"
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # JCALL-SPECIFIC OPTIONS
    jcall_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[parent_parser],
    )

    jcall_required = jcall_parser.add_argument_group("Required")
    jcall_optional = jcall_parser.add_argument_group("Optional")

    jcall_required.add_argument(
        "--genomicsdb_workspace_dir",
        help="""
        Path containiner GenomicsDB workspaces. Required if calling variants from, or imprting gVCFS to, 
        existing genomicsdb workspaces.
        """,
        required=False, # conditionally required if '--imprt' is NOT set.
        type=pathlib.Path,
        metavar="<path>"
    )

    jcall_optional.add_argument(
        "--imprt",
        help="""
        Add this flag if importing gVCF files to new genomicsdb workspaces. Cannot be used with '--update'.
        """,
        required=False,
        action="store_true"
    )

    jcall_optional.add_argument(
        "--update",
        help="""
        Add this flag if importing gVCF files to existing genomicsdb workspaces. Cannot be used with '--imprt'.
        """,
        required=False,
        action="store_true"
    )

    jcall_optional.add_argument(
        "--gvcf_dir",
        help="""
        Path containing gVCF files and their indexes ('.tbi').
        """,
        required=False, # conditionally required if '--update' is set
        type=pathlib.Path,
        metavar="<path>"
    )

    jcall_optional.add_argument(
        "--jcaller",
        help="""
        Joint sample variant caller; 'gatk', 'glnexus' [default: gatk].
        """,
        required=False,
        default='gatk',
        type=str,
        metavar="<text>"
    )

    jcall_optional.add_argument(
        "--interval",
        help="""
        List containing genomic intervals, one chromosome name per line and/or coordinate in bed format: <chr> <start> <stop>.
        This is only relevant for 'Dysgu'. HaplotypeCaller or DeepVariant calling per interval is not very efficient thousands of 
        intervals would need to be generated per sample.
        NB: Ensure that your chromosome names are the same as in the reference (e.g. chr1 or 1).
		If not provided, non-overlapping intervals of 5M bp will be generated from gVCF header.
        If '--update' is specified and interval list is provided, it must be the same as that used to import.
        """,
        required=False,
        type=pathlib.Path,
        metavar="<file>"
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FILTER-SPECIFIC OPTIONS
    filter_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[filter_annotate_common_parser],
    )

    filter_required = filter_parser.add_argument_group("Required")
    filter_optional = filter_parser.add_argument_group("Optional")

    filter_optional.add_argument(
        "--minDP",
        help="""
        Minimum allele depth [default: 10].
        """,
        required=False,
        default=10,
        type=int,
        metavar="<integer>"
    )

    filter_optional.add_argument(
        "--minGQ",
        help="""
        Minimum genotype quality [default: 20].
        """,
        required=False,
        default=20,
        type=int,
        metavar="<integer>"
    )

    filter_optional.add_argument(
        "--jcaller",
        help="""
        The tool that was used to generate joint call VCF file; 'gatk', 'glnexus' [default: gatk].
        """,
        required=False,
        default='gatk',
        type=str,
        metavar="<text>"
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ANNOTATE-SPECIFIC OPTIONS
    annotate_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[filter_annotate_common_parser],
    )

    annotate_required = annotate_parser.add_argument_group("Required")
    annotate_optional = annotate_parser.add_argument_group("Optional")

    annotate_optional.add_argument(
        "--interval",
        help="""
        List containing genomic intervals. E.g. one chromosome name per line and/or coordinate in bed format: <chr> <start> <stop>.
        NB: Ensure that your chromosome names are the same as in the reference (e.g. chr1 or 1).
        If not provided, full VCF file provided will be processed in one run. This will take
        longer for large files. It is recommended to provide interval list.
        """,
        required=False,
        type=pathlib.Path,
        metavar="<file>"
    )

    annotate_optional.add_argument(
        "--minimal",
        help="""
        Add this flag to perform only minimal annotation. Only a minimal set of databases are used 
        e.g., refGene, knownGene, cytoBand, avsnp156. This is useful for annotating GWAS results.
        """,
        required=False,
        action="store_true"
    )

    annotate_optional.add_argument(
        "--alt_contig",
        help="""
        Add this flag to include alternate contigs in the annotation. Alternate contigs are removed by default. Only chomosomes 1-22,M,X,Y are processed.
        """,
        required=False,
        action="store_true"
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MAIN PAERSER
    parser = argparse.ArgumentParser(
        description=descmsg,
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )

    parser_optional = parser.add_argument_group("Help")

    parser_optional.add_argument(
        "-h", "--help",
        action="help",
        help=argparse.SUPPRESS
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # SUBPARSER
    # add subparsers for commands
    subparsers = parser.add_subparsers(
        title="Commands", 
        dest='command',
        help=argparse.SUPPRESS,
        metavar="",
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # TEST SUBPARSER 
    subparsers.add_parser(
        "test",
        help=argparse.SUPPRESS
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # QC SUBPARSER    
    subparsers.add_parser(
        "qc",
        prog=prog,
        usage="%(prog)s qc [-h/--help] <options>",
        parents=[qc_trim_align_common_parser],
        description="READS QUALITY ASSESSMENT",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # TRIM SUBPARSER 
    subparsers.add_parser(
        "trim",
        prog=prog,
        usage="%(prog)s trim [-h/--help] <options>",
        parents=[trim_parser],
        description="READS TRIMMING",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ALIGN SUBPARSER     
    subparsers.add_parser(
        "align",
        prog=prog,
        usage="%(prog)s align [-h/--help] <options>",
        parents=[align_parser],
        description="READS ALIGNMENT",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MERGE ALIGN SUBPARSER
    subparsers.add_parser(
        "mergealign",
        prog=prog,
        usage="%(prog)s mergealign [-h/--help] <options>",
        parents=[mergealign_parser],
        description="ALIGMENT FILES MERGE",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # call SUBPARSER
    subparsers.add_parser(
        "call",
        prog=prog,
        usage="%(prog)s call [-h/--help] <options>",
        parents=[call_parser],
        description="ONE-RUN SINGLE AND JOINT VARIANT CALLING",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Scall SUBPARSER
    subparsers.add_parser(
        "scall",
        prog=prog,
        usage="%(prog)s scall [-h/--help] <options>",
        parents=[scall_parser],
        description="SINGLE SAMPLE VARIANT CALLING",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Jcall SUBPARSER
    subparsers.add_parser(
        "jcall",
        prog=prog,
        usage="%(prog)s jcall [-h/--help] <options>",
        parents=[jcall_parser],
        description="JOINT-SAMPLE VARIANT CALLING",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # filter SUBPARSER
    subparsers.add_parser(
        "filter",
        prog=prog,
        usage="%(prog)s filter [-h/--help] <options>",
        parents=[filter_parser],
        description="VARIANT FILTERATION",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ANNOTATE SUBPARSER
    subparsers.add_parser(
        "annotate",
        prog=prog,
        usage="%(prog)s annotate [-h/--help] <options>",
        parents=[annotate_parser],
        description="VARIANT ANNOTATION",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # RAW TO VCF SUBPARSER
    subparsers.add_parser(
        "raw2vcf",
        prog=prog,
        usage="%(prog)s raw2vcf [-h/--help] <options>",
        parents=[align_parser],
        description="""
        READS ALIGNMENT TO VARIANT CALLING:
        This workflow will run from alignment to variant calling
        """,
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # RAW TO ANNOTATED SUBPARSER
    subparsers.add_parser(
        "raw2annotate",
        prog=prog,
        usage="%(prog)s raw2annotate [-h/--help] <options>",
        parents=[align_parser],
        description="""
        READS ALIGNMENT TO ANNOTATED VCF:    
        This workflow will run from alignment through variant calling to generate an annotated VCF
        """,
        add_help=False
    )

    args = parser.parse_args()

    return args


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET ONE-TIME SYSTEM SETTINGS
"""
We will either use the template 'system.json' file to provide global systems configurations by
copying the file to a 'nextflow.json' file and filling in the parameter values, or by simply running 
the 'ei-ngs.py' script on the command line and responding to the prompts to provide values to the required
system parameter settings.
"""
def system_settings(
        account=None, 
        partition=None,
        project_name=None,
        containsers_dir=None,
        workspace=None,
        email=None,
        project_dir=None,
    ):
    # check if 'nextflow.config' exists and get information from it
    main_config_file = pathlib.Path('nextflow.json')
    if os.path.exists(main_config_file):
        main_config_json = open(main_config_file, 'r')
        main_config = json.load(main_config_json)
        account = main_config['account']
        partition = main_config['partition']
        project_name = main_config['project_name']
        containsers_dir = main_config['containsers_dir']
        workspace = main_config['workspace']
        email = main_config['email']
        if account == None or account == "NULL" or account == "" or \
            partition == None or partition == "NULL" or partition == "" or \
            project_name == None or project_name == "NULL" or project_name == "" or \
            containsers_dir == None or containsers_dir == "NULL" or containsers_dir == "" or \
            workspace == None or workspace == "NULL" or workspace == "":
            logger.info("Please provide one-time system settings")
            account = str(input("Enter your cluster group account name: "))
            partition = str(input("Enter the partition/queue you want jobs to be submitted to: "))
            project_name = str(input("Enter the name of this project.\n[Note, this name will be used to construct output file names. Therefore, set a different name for different projects]: "))
            containsers_dir = str(input("Enter path where containers will be stored:\n[Note, this should have sufficient storage capacity since all containers needed could be many gigabytes in size]: "))
            workspace = str(input("Enter path where all analyses will be staged\n[Note, this must sufficient storage capacity since NGS intermediate files and outputs could be many terabytes in size]: "))
            email = str(input("Provide your email address to receive notifications of job status.\n[Leave empty if you do not wish to receive email notifications]: "))
            if email == "":
                email = "NULL"
            main_config_dic = {
                "account": account,
                "partition": partition,
                "project_name": project_name,
                "containsers_dir": containsers_dir,
                "workspace": workspace,
                "email": email
            }
            with open(main_config_file, 'w') as main_config_json:
                json.dump(main_config_dic, main_config_json, indent=4)
    else:
        logger.info("Please provide base system settings")
        account = input("Enter your cluster group account name: ")
        partition = input("Enter the partition/queue you want jobs to be submitted to: ")
        project_name = input("Enter the name of this project.\n[Note, this name will be used to construct output file names. Therefore, set a different name for different projects]: ")
        containsers_dir = str(input("Enter path where containers will be stored:\n[Note, this should have sufficient storage capacity since all containers needed could be many gigabytes in size]: "))
        workspace = str(input("Enter path where all analyses will be staged\n[Note, this must sufficient storage capacity since NGS intermediate files and outputs could be many terabytes in size]: "))
        email = str(input("Provide your email address to receive notifications of job status. \n[Leave empty if you do not wish to receive email notifications]: "))
        if email == "":
            email = "NULL"
        main_config_dic = {
            "account": account,
            "partition": partition,
            "project_name": project_name,
            "containsers_dir": containsers_dir,
            "workspace": workspace,
            "email": email
        }
        with open(main_config_file, 'w') as main_config_json:
            json.dump(main_config_dic, main_config_json, indent=4)

    # GENERATE NEXTFLOW CONFIG FROM MAIN CONFIG DICTIONARY
    nf_config_name = 'nextflow.config'
    nf_config = open(nf_config_name, 'w')
    nf_config.write('params { // project-specific one-time system configuration //' + "\n")
    nf_config.write(f"  account = '{account}'" + "\n")
    nf_config.write(f"  queue = '{partition}'" + "\n")
    nf_config.write(f"  project_name = '{project_name}'" + "\n")
    nf_config.write(f"  containers_dir = '{containsers_dir}'" + "\n")
    nf_config.write(f"  workspace = '{workspace}'" + "\n")
    nf_config.write(f"  email = '{email}'" + "\n")
    nf_config.write("}" + "\n")

def get_project_config(
        dtype=None, 
        cmd=None
    ):
    # GENERATE MAIN PROJECT CONFIG
    main_config = json.load(open('nextflow.json', 'r'))
    workspace = main_config['workspace'] + f'/{main_config["project_name"]}'
    project_name = main_config['project_name']
    project_config_name = project_name + '.config'
    project_config = open(project_config_name, 'w')
    project_config.write("includeConfig \"${launchDir}/nextflow.config\"\n")
    project_config.write("params {" + "\n")
    project_config.write(f"  wkflow = '{cmd}'" + "\n")
    project_config.write("}" + "\n")
    project_config.write("includeConfig \"${projectDir}/configs/profile-selector.config\"\n")
    if dtype.upper() == "WGS":
        rselector = "includeConfig \"${projectDir}/configs/resourceselector/resource-selector-wgs.config\""
        project_config.write(f"{rselector}\n")
    else:
        rselector = "includeConfig \"${projectDir}/configs/resourceselector/resource-selector-wes.config\""
        project_config.write(f"{rselector}\n")

    return workspace, project_name, project_config_name

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAKE NEXTFLOW CONFIGURATION FILES
# MAKE TEST CONFIG FILE
def get_test_config():
    #check if exist and remove
    test_config_file = pathlib.Path('test.config')
    if os.path.exists(test_config_file):
        os.remove(test_config_file)
    
    #create new
    with open(test_config_file, 'a') as f:
        f.write("includeConfig \"${launchDir}/nextflow.config\"\n")
        f.write("params { wkflow = 'test' }\n")
        f.write("includeConfig \"${projectDir}/configs/profile-selector.config\"\n")
        f.write("includeConfig \"${projectDir}/configs/test.config\"\n")

