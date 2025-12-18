import subprocess
import wfstaging
import nextflow
from pandas.api.types import is_string_dtype

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DEFINE WORKFLOW FUNCTIONS
#
# TEST WORKFLOW
def test_workflow(
        project_dir=None,
        logger=None
    ):
    wfstaging.get_test_config()

    # Use subprocess to run test workflow 
    nf_test = subprocess.run(
        [
            str("nextflow"),
            str("-c"),
            str("test.config"),
            str("run"),
            #str(f"{project_dir}/workflows/test.nf"),
            str(f"{project_dir}/main.nf"),
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
        #print(nf_test.stdout)
    else:
        logger.error("Workflow test terminated with an error.")
        print(nf_test.stderr)



# QC WORKFLOW
def qc_workflow(
        args=None,
        run_id=None,
        workspace=None, 
        project_name=None, 
        project_config=None,
        project_dir=None
    ):
    qc_params = {
        "input_ftype": f"{args.ftype}",
        "input_dir": f"{args.input_dir}",
        "output_dir": f"{args.output_dir}",
        "threads": f"{args.threads}",
        "njobs": f"{args.njobs}"
    }

    nf_qc = nextflow.run(
        f"{project_dir}/main.nf",
        params=qc_params,
        run_path=".",
        output_path=f"{workspace}",
        log_path=".",
        profiles=[
            f"{args.profile}"
        ],
        configs=[
            f"{project_config}"
        ],
        resume=f"{args.resume}",
        report=f"{project_name}-qc-{run_id}-report.html",
        timeline=f"{project_name}-qc-{run_id}-timeline.html",
        dag=f"{project_name}-qc-{run_id}-dag.html",
        trace=f"{project_name}-qc-{run_id}-trace.txt"
    )
    print(f"NEXTFLOW SESSION STATUS: {nf_qc.status}")
    print(nf_qc.stderr)
    print(nf_qc.stdout)


# TRIM WORKFLOW
def trim_workflow(
        args=None,
        run_id=None,
        workspace=None, 
        project_name=None, 
        project_config=None,
        project_dir=None
    ):
    trim_params = {
        "input_ftype": f"{args.ftype}",
        "input_dir": f"{args.input_dir}",
        "output_dir": f"{args.output_dir}",
        "threads": f"{args.threads}",
        "njobs": f"{args.njobs}",
        "trimmer": f"{args.trimmer}",
        "adapter": f"{args.adapter}",
        "min_length": f"{args.min_length}",
        "headcrop": f"{args.headcrop}",
        "crop": f"{args.crop}",
    }

    # delete adapter from trim parameters if trimgalore is selected
    if args.trimmer == 'trimgalore':
        del trim_params['adapter']
        # hard code max threads for trimgalore to 8,
        # maximum recommended by Babraham Institute
        if args.threads > 8:
            trim_params['threads'] = 8

    if args.headcrop == 0:
        trim_params['headcrop'] = "NULL"
    if args.crop == 0:
        trim_params['crop'] = "NULL"

    if args.resume == 'false':
        nf_trim = nextflow.run(
            f"{project_dir}/main.nf",
            params=trim_params,
            run_path=".",
            output_path=f"{workspace}",
            log_path=".",
            profiles=[
                f"{args.profile}"
            ],
            configs=[
                f"{project_config}"
            ],
            report=f"{project_name}-trim-{run_id}-report.html",
            timeline=f"{project_name}-trim-{run_id}-timeline.html",
            dag=f"{project_name}-trim-{run_id}-dag.html",
            trace=f"{project_name}-trim-{run_id}-trace.txt"
        )
    else:
        if args.resume == 'true':
            resume_val = True
        else:
            resume_val = f"{args.resume}"

        nf_trim = nextflow.run(
            f"{project_dir}/main.nf",
            params=trim_params,
            run_path=".",
            output_path=f"{workspace}",
            log_path=".",
            profiles=[
                f"{args.profile}"
            ],
            configs=[
                f"{project_config}"
            ],
            resume=resume_val,
            report=f"{project_name}-trim-{run_id}-report.html",
            timeline=f"{project_name}-trim-{run_id}-timeline.html",
            dag=f"{project_name}-trim-{run_id}-dag.html",
            trace=f"{project_name}-trim-{run_id}-trace.txt"
        )

    print(f"NEXTFLOW SESSION STATUS: {nf_trim.status}")
    print(nf_trim.stderr)
    print(nf_trim.stdout)


# ALIGN WORKFLOW
def align_workflow(
        args=None,
        run_id=None,
        workspace=None, 
        project_name=None, 
        project_config=None,
        project_dir=None
    ):

    if args.se:
        pe = 'false'
    else:
        pe = 'true'

    align_params = {
        "input_ftype": f"{args.ftype}",
        "input_dir": f"{args.input_dir}",
        "output_dir": f"{args.output_dir}",
        "threads": f"{args.threads}",
        "njobs": f"{args.njobs}",
        "aligner": f"{args.aligner}",
        "pe": f"{pe}",
        "dup_marker": f"{args.dup_marker}",
        "remove_dup": f"{args.remove_dup}",
        "spark": "true",
    }

    if not args.spark:
        del align_params['spark']

    nf_align = nextflow.run(
        f"{project_dir}/main.nf",
        params=align_params,
        run_path=".",
        output_path=f"{workspace}",
        profiles=[
            f"{args.profile}"
        ],
        configs=[
            f"{project_config}"
        ],
        resume=f"{args.resume}",
        report=f"{project_name}-align-{run_id}-report.html",
        timeline=f"{project_name}-align-{run_id}-timeline.html",
        dag=f"{project_name}-align-{run_id}-dag.html",
        trace=f"{project_name}-align-{run_id}-trace.txt"
    )
    print(f"NEXTFLOW SESSION STATUS: {nf_align.status}")
    print(nf_align.stderr)
    print(nf_align.stdout)


# MERGEALIGN WORKFLOW
def mergealign_workflow(
        args=None,
        run_id=None,
        workspace=None, 
        project_name=None, 
        project_config=None,
        project_dir=None
    ):

    mergealign_params = {
        "input_dir": f"{args.input_dir}",
        "output_dir": f"{args.output_dir}",
        "threads": f"{args.threads}",
        "njobs": f"{args.njobs}",
        "sort_order": f"{args.sort_order}"
    }

    nf_mergealign = nextflow.run(
        f"{project_dir}/main.nf",
        params=mergealign_params,
        run_path=".",
        output_path=f"{workspace}",
        profiles=[
            f"{args.profile}"
        ],
        configs=[
            f"{project_config}"
        ],
        resume=f"{args.resume}",
        report=f"{project_name}-mergealign-{run_id}-report.html",
        timeline=f"{project_name}-mergealign-{run_id}-timeline.html",
        dag=f"{project_name}-mergealign-{run_id}-dag.html",
        trace=f"{project_name}-mergealign-{run_id}-trace.txt"
    )
    print(f"NEXTFLOW SESSION STATUS: {nf_mergealign.status}")
    print(nf_mergealign.stderr)
    print(nf_mergealign.stdout)