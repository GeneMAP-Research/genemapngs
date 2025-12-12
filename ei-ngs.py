#!/usr/bin/env python

from runpy import run_path
import subprocess
import pathlib
import json
from email.policy import default
import importlib
import logging
import os
import warnings
import shutil
import sys
import argparse
import psutil
import time
import textwrap
import nextflow
import random
import string
#import pandas as pd
#import numpy as np
#from pandas.api.types import is_integer_dtype
#from pandas.api.types import is_float_dtype
#from pandas.api.types import is_string_dtype
#import statsmodels.stats.multitest as smm

warnings.filterwarnings('ignore')

"""
REQUIREMTNS:
    - psutil
"""

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
# Configure basic logging to console
#logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(message)s')

# Get a logger instance
logger = logging.getLogger(__name__)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE UNIQUE UD FOR EACH RUN
def generate_random_string(length):
    # Define the possible characters: only letters and digits in this example
    characters = string.ascii_letters + string.digits
    # Use random.choices() which is efficient for generating multiple selections with replacement
    random_string = ''.join(random.choices(characters, k=length))
    return random_string

run_id = generate_random_string(length=10)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get script name
script_path = sys.argv[0]
script_name = os.path.basename(script_path)

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
# CHECK REQUIRED PACKAGES
def is_bash_tool_installed(tool_name):
    """
    Checks if a given Bash tool is installed and accessible in the system's PATH.
    """
    return shutil.which(tool_name) is not None

# ADD BASH TOOLS THAT NEED TO BE CHECKED HERE
tools = [
    "nextflow",
    "singularity",
    "ls"
]

for tool in tools:
    if is_bash_tool_installed(tool):
        next
        #print("parallel is installed/loaded.")
    else:
        sys.exit(f"\nError: {tool} is not installed/loaded. Exiting...\n")

# CHECK PYTHON PACKAGE INSTALLATION
# List of module names you want to import
# Assumes these modules are in the same directory or on the Python path
module_names = ['psutil']

# A dictionary to hold the imported modules
imported_modules = {}

for module_name in module_names:
    try:
        # Import the module using importlib.import_module()
        module = importlib.import_module(module_name)
        # Store the imported module in the dictionary
        imported_modules[module_name] = module
        #print(f"Successfully imported {module_name}")
    except ImportError as e:
        print(f"Error importing {module_name}: {e}")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define banner and usage messages
descmsg = f"""
    {ANSIGRY}#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#{ANSIRESET}
    {grey_pipe}  {ANSIYLW}WRAPPER FOR THE {ANSILBL}ESOH{ANSIRED}INFORMATICS {ANSIPPL}NGS {ANSIGRN}NEXTFLOW {ANSIYLW}WORKFLOW{ANSIRESET}  {grey_pipe}
    {grey_pipe}                    Kevin Esoh, 2025                     {grey_pipe}
    {grey_pipe}                   kesohku1@jhmi.edu                     {grey_pipe}
    {ANSIGRY}#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#{ANSIRESET}
"""

version="0.1 (beta)"

usage = f"""
{descmsg}

LICENSCE: GNU GPLv3+, due to use of the GNU Scientific Library
VERSION: {version}

Usage: {script_name} <command> [-v/--version] [-h/--help] <options>

Commands:
           test: Run test to see if workfow installed correctly.
             qc: Check FASTQ or Alignment (BAM/CRAM) quality.
           trim: Trim adapters and poor quality bases from reads.
          align: Align/map reads to reference and post-alignment processing.
     mergealign: Megre Alignment (BAM/CRAM) files.
        varcall: Perform variant calling (both single and joint sample) in one run.
       svarcall: Perform only sinlge sample variant calling to generate gVCF files.
       jvarcall: Perform only joint (multi-sample) variant calling with pre-existing gVCF files.
      varfilter: Filter variant calls in VCF/BCF files.
       annotate: Annotate variants with ANNOVAR
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DEFINE ARGUMENTS
def get_arguments(descmsg=None, prog=None):
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
        "-resume","--resume",
        help="""
        Whether to resume the nextflow execution. Either the flag without any value (--resume) to resume from last
        exit point or add the flag with the specific RUN NAME as can be obtained using `nexflow log` to resume from
        a specific point (--resume kickass_church).
        """,
        required=False,
        nargs='?',
        const="True",
        default="False",
        type=str
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
    # VARFILTER and ANNOTATE SHARED ARGUMENTS
    varfilter_annotate_common_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[parent_parser],
        formatter_class=argparse.RawTextHelpFormatter,
        #formatter_class=argparse.RawDescriptionHelpFormatter
    )

    varfilter_annotate_common_required = varfilter_annotate_common_parser.add_argument_group("Required")
    varfilter_annotate_common_optional = varfilter_annotate_common_parser.add_argument_group("Optional")

    varfilter_annotate_common_required.add_argument(
        "--vcf_dir",
        help="""
        Path containing VCF file(s).
        """,
        required=True,
        type=pathlib.Path,
        metavar="<path>"
    )

    varfilter_annotate_common_optional.add_argument(
        "--left_norm",
        help="""
        Add this flag to left-normalized variants such as is recommended by ANNOVAR.
        """,
        required=False,
        action="store_true"
    )

    varfilter_annotate_common_optional.add_argument(
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
    # VARCALL-SPECIFIC OPTIONS
    varcall_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[parent_parser],
    )

    varcall_required = varcall_parser.add_argument_group("Required")
    varcall_optional = varcall_parser.add_argument_group("Optional")

    varcall_required.add_argument(
        "--alignment_dir",
        help="Path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).",
        required=True,
        type=pathlib.Path,
        metavar="<path>"
    )

    varcall_optional.add_argument(
        "--scaller",
        help="""
        Single sample variant caller; 'gatk', 'deepvariant' [default: gatk].
        'gatk' will use GATK HaplotypeCaller for single sample calling.
        For single sample calling of structural variants, use the 'svarcall' command.
        """,
        required=False,
        default='gatk',
        type=str,
        metavar="<text>"
    )

    varcall_optional.add_argument(
        "--jcaller",
        help="""
        Joint sample variant caller; 'gatk', 'glnexus' [default: gatk].
        """,
        required=False,
        default='gatk',
        type=str,
        metavar="<text>"
    )

    varcall_optional.add_argument(
        "--batch_size",
        help="""
        Number of samples to read into memory by GATK sample reader each time [default: 50].
        """,
        required=False,
        default=50,
        type=int,
        metavar="<integer>"
    )

    varcall_optional.add_argument(
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
    # SVARCALL-SPECIFIC OPTIONS
    svarcall_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[parent_parser],
    )

    svarcall_required = svarcall_parser.add_argument_group("Required")
    svarcall_optional = svarcall_parser.add_argument_group("Optional")

    svarcall_required.add_argument(
        "--alignment_dir",
        help="Path to alignment (BAM/CRAM) files and their indexes (.bai/.crai).",
        required=True,
        type=pathlib.Path,
        metavar="<path>"
    )

    svarcall_required.add_argument(
        "--vcf_dir",
        help="Path to VCF files and their indexes (.tbi) for 'dysgu merge' only!. NOTE: You ust specify '--ftype VCF'",
        required=False, # will be conditionally required if --ftype is VCF
        type=pathlib.Path,
        metavar="<path>"
    )

    svarcall_optional.add_argument(
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

    svarcall_optional.add_argument(
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

    svarcall_optional.add_argument(
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
    # JVARCALL-SPECIFIC OPTIONS
    jvarcall_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[parent_parser],
    )

    jvarcall_required = jvarcall_parser.add_argument_group("Required")
    jvarcall_optional = jvarcall_parser.add_argument_group("Optional")

    jvarcall_required.add_argument(
        "--genomicsdb_workspace_dir",
        help="""
        Path containiner GenomicsDB workspaces. Required if calling variants from, or imprting gVCFS to, 
        existing genomicsdb workspaces.
        """,
        required=False, # conditionally required if '--imprt' is NOT set.
        type=pathlib.Path,
        metavar="<path>"
    )

    jvarcall_optional.add_argument(
        "--imprt",
        help="""
        Add this flag if importing gVCF files to new genomicsdb workspaces. Cannot be used with '--update'.
        """,
        required=False,
        action="store_true"
    )

    jvarcall_optional.add_argument(
        "--update",
        help="""
        Add this flag if importing gVCF files to existing genomicsdb workspaces. Cannot be used with '--imprt'.
        """,
        required=False,
        action="store_true"
    )

    jvarcall_optional.add_argument(
        "--gvcf_dir",
        help="""
        Path containing gVCF files and their indexes ('.tbi').
        """,
        required=False, # conditionally required if '--update' is set
        type=pathlib.Path,
        metavar="<path>"
    )

    jvarcall_optional.add_argument(
        "--jcaller",
        help="""
        Joint sample variant caller; 'gatk', 'glnexus' [default: gatk].
        """,
        required=False,
        default='gatk',
        type=str,
        metavar="<text>"
    )

    jvarcall_optional.add_argument(
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
    # VARFILTER-SPECIFIC OPTIONS
    varfilter_parser = argparse.ArgumentParser(
        add_help=False,
        parents=[varfilter_annotate_common_parser],
    )

    varfilter_required = varfilter_parser.add_argument_group("Required")
    varfilter_optional = varfilter_parser.add_argument_group("Optional")

    varfilter_optional.add_argument(
        "--minDP",
        help="""
        Minimum allele depth [default: 10].
        """,
        required=False,
        default=10,
        type=int,
        metavar="<integer>"
    )

    varfilter_optional.add_argument(
        "--minGQ",
        help="""
        Minimum genotype quality [default: 20].
        """,
        required=False,
        default=20,
        type=int,
        metavar="<integer>"
    )

    varfilter_optional.add_argument(
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
        parents=[varfilter_annotate_common_parser],
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
    # VARCALL SUBPARSER
    subparsers.add_parser(
        "varcall",
        prog=prog,
        usage="%(prog)s varcall [-h/--help] <options>",
        parents=[varcall_parser],
        description="ONE-RUN SINGLE AND JOINT VARIANT CALLING",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # SVARCALL SUBPARSER
    subparsers.add_parser(
        "svarcall",
        prog=prog,
        usage="%(prog)s svarcall [-h/--help] <options>",
        parents=[svarcall_parser],
        description="SINGLE SAMPLE VARIANT CALLING",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # JVARCALL SUBPARSER
    subparsers.add_parser(
        "jvarcall",
        prog=prog,
        usage="%(prog)s jvarcall [-h/--help] <options>",
        parents=[jvarcall_parser],
        description="JOINT-SAMPLE VARIANT CALLING",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # VARFILTER SUBPARSER
    subparsers.add_parser(
        "varfilter",
        prog=prog,
        usage="%(prog)s varfilter [-h/--help] <options>",
        parents=[varfilter_parser],
        description="VARIANT FILTERATION",
        add_help=False
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ANNOTATE SUBPARSER
    subparsers.add_parser(
        "annotate",
        prog=prog,
        usage="%(prog)s qc [-h/--help] <options>",
        parents=[annotate_parser],
        description="VARIANT ANNOTATION",
        add_help=False
    )

    args = parser.parse_args()

    return args


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET ONE-TIME SYSTEM SETTINGS
def system_settings(
        account=None, 
        partition=None,
        project_name=None,
        containsers_dir=None,
        workspace=None,
        email=None,
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

    # GENERATE NEXTFLOW CONFIG
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

def get_project_config(dtype=None):
    # GENERATE MAIN PROJECT CONFIG
    main_config = json.load(open('nextflow.json', 'r'))
    project_config_name = main_config['project_name'] + '.config'
    project_config = open(project_config_name, 'w')
    project_config.write("includeConfig \"${projectDir}/nextflow.config\"\n")
    project_config.write("includeConfig \"${projectDir}/configs/profile-selector.config\"\n")            
    if dtype.upper() == "WGS":
        rselector = 'includeConfig "${projectDir}/configs/resourceselector/resource-selector-wgs.config"'
        project_config.write(f"{rselector}\n")
    else:
        rselector = 'includeConfig "${projectDir}/configs/resourceselector/resource-selector-wes.config"'
        project_config.write(f"{rselector}\n")    



# MAKE NEXTFLOW CONFIGURATION FILES
# MAKE TEST CONFIG FILE
def get_test_config():
    #check if exist and remove
    test_config_file = pathlib.Path('test.config')
    if test_config_file.is_file:
        test_config_file.unlink()
    
    #create new
    with open(test_config_file, 'a') as f:
        f.write("includeConfig \"${projectDir}/nextflow.config\"\n")
        f.write("includeConfig \"${projectDir}/configs/profile-selector.config\"\n")
        f.write("includeConfig \"${projectDir}/configs/test.config\"\n")


if __name__ == "__main__":
    #monitor_resources()
    args = get_arguments(
        descmsg=usage, 
        prog=script_name
    )

    system_settings() 

    if not args.command:
        print(usage)
    elif args.command == 'test':
        print("Testing if ei-ngs nextflow workflow installed successfully...")
        get_test_config()
        nf_test = subprocess.run(
            [
                str("nextflow"),
                str("-c"),
                str("test.config"),
                str("run"),
                str("test.nf"),
                str("-params-file"),
                str("nextflow.json"),
                str("-profile"), 
                str("singularity")

            ], 
            capture_output=True,
            text = True
        )
        if nf_test.stdout:
            print("\n!!!WORKFLOW TEST SUCCESSFUL!!!\n")
            print(nf_test.stdout)
        else:
            logger.error("Workflow test terminated with an error.")
            print(nf_test.stderr)
    else:
        if args.wgs:
            dtype = "WGS"
        else:
            dtype = "WES"

        get_project_config(dtype=dtype)

        # get workspace from 'nextflow.json' to use as workdir
        main_config = json.load(open('nextflow.json', 'r'))
        workspace = main_config['workspace'] + f'/{main_config["project_name"]}'
        # get project name to build project config from
        project_name = main_config['project_name']
        project_config = project_name + '.config'

        if args.command == 'qc':
            print("READS QUALITY ASSESSMENT")
            qc_params = {
                "input_ftype": f"{args.ftype}",
                "input_dir": f"{args.input_dir}",
                "output_dir": f"{args.output_dir}",
                "threads": f"{args.threads}",
                "njobs": f"{args.njobs}"
            }
            """
            qc_job = subprocess.run(
                [
                    str("nextflow"),
                    str("-c"),
                    str(f"{project_config}"),
                    str("run"),
                    str("getQualityReports.nf"),
                    str("-profile"),
                    str(f"{args.profile}"),
                    str("-w"),
                    str(f"{workspace}"),
                    str("-with-report"),
                    str("--input_ftype"),
                    str(f"{args.ftype}"),
                    str("--input_dir"),
                    str(f"{args.input_dir}"),
                    str("--output_dir"),
                    str(f"{args.output_dir}"),
                    str("--threads"),
                    str(f"{args.threads}"),
                    str("--njobs"),
                    str(f"{args.njobs}")
                ],
                capture_output = True,
                text = True,
                check = True
                )
            if qc_job.stdout:
                print("\n!!!QC WORKFLOW SUCCESSFUL!!!\n")
                print(qc_job.stdout)
            else:
                logger.error("QC workflow terminated with an error.")
                print(qc_job.stderr)
            """
            os.makedirs(
                workspace,
                exist_ok=True
            )
            execution = nextflow.run(
                "getQualityReports.nf",
                params=qc_params,
                run_path=".",
                output_path=f"{workspace}",
                profiles=[
                    f"{args.profile}"
                ],
                configs=[
                    f"{project_config}"
                ],
                resume=f"{args.resume}",
                report=f"{project_name}-{run_id}-report.html",
                timeline=f"{project_name}-{run_id}-timeline.html",
                dag=f"{project_name}-{run_id}-dag.html",
                trace=f"{project_name}-{run_id}-trace.txt"
            )
            print(f"SESSION STATUS: {execution.status}")
            print(f"SESSION UNIQUE ID (UUID): {execution.uuid}")
            print(f"SESSION START: {execution.start}")
            print(f"SESSION FINISHED: {execution.finished}")
            print(f"SESSION COMMAND: {execution.command}")
            print(f"SESSION EXECUTION PATH: {execution.path}")
            print(execution.stderr)
            print(execution.stdout)
