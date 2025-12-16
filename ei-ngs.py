#!/usr/bin/env python

from runpy import run_path
from email.policy import default
import logging
import os
import warnings
import sys
#import subprocess
#import pathlib
#import json
#import importlib
#import shutil
#import argparse
#import psutil
#import time
#import textwrap
#import nextflow
#import random
#import string
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
    - nextflowpy
"""

# Get the absolute path of the directory containing the module based on installation path
wrapper_script = sys.argv[0]
script_name = os.path.basename(wrapper_script)
script_path = os.path.dirname(wrapper_script)
module_dir = os.path.abspath(f'{script_path}/modules')

# Add the directory to sys.path
sys.path.append(module_dir)

# Version
version="0.1 (beta)"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Configure basic logging to console
#logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(message)s')

# Get a logger instance
logger = logging.getLogger(__name__)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get script name
wrapper_script = sys.argv[0]
script_name = os.path.basename(wrapper_script)
script_path = os.path.dirname(wrapper_script)


# Bash tools to check
tools = [
    "nextflow",
    "singularity",
    "ls"
]

# List of module names to import
module_names = [
    'psutil',
    'nextflow'
]


if __name__ == "__main__":
    # import the custom module
    from wfstaging import ( 
        get_run_id, 
        is_bash_tool_installed, 
        check_python_package, 
        get_usage, 
        get_arguments, 
        system_settings, 
        get_project_config, 
    )

    from wfexecution import ( 
        test_workflow, 
        qc_workflow,
        trim_workflow
    )

    # Get RUN ID
    run_id = get_run_id(length=10)

    # CHECK BASH TOOL INSTALLATION
    for tool in tools:
        if not is_bash_tool_installed(tool):
            sys.exit(f"\nError: {tool} is not installed/loaded. Exiting...\n")

    # CHECK PYTHON PACKAGE INSTALLATION
    check_python_package(module_list=module_names)
    #monitor_resources()

    usage = get_usage(
        script_name=script_name,
        version=version        
    )

    args = get_arguments(
        descmsg=usage, 
        prog=script_name,
        version=version
    )

    system_settings() 

    if not args.command:
        print(usage)
    elif args.command == 'test':
        print("Testing if ei-ngs nextflow workflow installed successfully...")
        test_workflow(
            project_dir=script_path
        )
    else:
        if args.wgs:
            dtype = "WGS"
        else:
            dtype = "WES"

        # GET VALUES FROM 'get_project_config' FUNCTION
        workspace, project_name, project_config = get_project_config(dtype=dtype, cmd=args.command)
    
        os.makedirs(
            workspace,
            exist_ok=True
        )

        print(f"PROJECT NAME: {project_name}")
        print(f"JOB ID: {run_id}")

        if args.command == 'qc':
            print("READS QUALITY ASSESSMENT")
            qc_workflow(
                args=args,
                run_id=run_id,                
                workspace=workspace, 
                project_name=project_name, 
                project_config=project_config,
                project_dir=script_path
            )
            
        if args.command == 'trim':
            print("READS TRIMMING")
            trim_workflow(
                args=args,
                run_id=run_id,                
                workspace=workspace, 
                project_name=project_name, 
                project_config=project_config,
                project_dir=script_path
            )


