task gs_run_lcnv_pipeline {

    Array[String]? samples
    File profilesArchive
    File gsReferenceBundle

    Int memory
    Int diskSize
    Int numThreads
    Int numPreempt
 
    command {
        $SV_DIR/scripts/firecloud/gs_run_lcnv_pipeline.sh ${write_lines(samples)} ${profilesArchive} ${gsReferenceBundle}
    }

    output {
        File lcnvCalls = "lcnv_calls.tar.gz"
    }

    runtime {
        docker: "skashin/genome-strip:dev"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        cpu: "${numThreads}"
        preemptible: "${numPreempt}"
    }


    meta {
        author: "Seva Kashin"
    }
}

workflow gs_run_lcnv_pipeline_workflow {
    call gs_run_lcnv_pipeline
}
