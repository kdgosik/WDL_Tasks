task gs_merge_metadata {
    Array[File] mdArchives
    File gsReferenceBundle

    Int memory
    Int diskSize
    Int numThreads
    Int numPreempt

    String mdArchivesList = "md_archives.list"

    command {
        echo ${sep=' ' mdArchives} | sed 's/ /\n/g' > ${mdArchivesList}
        $SV_DIR/scripts/firecloud/gs_merge_metadata.sh ${mdArchivesList} ${gsReferenceBundle}
    }

    output {
        File zippedMetadata = "metadata.zip"
    }

    runtime {
        docker: "skashin/genome-strip:latest"

        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        cpu: "${numThreads}"
        preemptible: "${numPreempt}"
    }

    meta {
        author: "Seva Kashin"
    }
}

workflow gs_merge_metadata_workflow {
    call gs_merge_metadata
}