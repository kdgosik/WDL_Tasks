task preprocess {
    File bamFile
    File bamIndex
    File referenceBundle
    String sampleId

    Int memory
    Int numPreempt

    Int diskSize = round(size(bamFile, "G")) + 40

    command {
        $SV_DIR/scripts/firecloud/preprocess.sh ${sampleId} ${bamFile} ${referenceBundle}
    }
    
    output {
        File mdPath = "metadata.zip"
    }

    runtime {
        docker: "skashin/genome-strip:latest"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

workflow gs_preprocessing_wf {
    call preprocess
}
