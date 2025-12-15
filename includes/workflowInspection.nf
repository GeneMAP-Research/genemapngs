def jobCompletionMessage() {
   return msg = """
          Pipeline execution summary
          ---------------------------
          Completed at: ${workflow.complete}
          Duration    : ${workflow.duration}
          exit status : ${workflow.success ? 'OK' : 'failed' }
          Success     : ${workflow.success}
          workDir     : ${workflow.workDir}
          outputDir   : ${params.output_dir}
          """
          .stripIndent()
}
