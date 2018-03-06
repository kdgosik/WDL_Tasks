task VIZTASK{
    File pvcf
    File ref
    File tumorBAM
    File controlBAM
    File tumorBAMIndex
    File controlBAMIndex
    String samplename
    
    runtime{
        docker : "erictdawson/svviz"
        memory : "28 GB"
        cpu : "4"
        disks : "local-disk 1000 HDD"
    }


    command <<<
        mkdir results_${samplename} && svviz --export results_${samplename}/ --format pdf --type batch --summary ${samplename}.summary.tsv -b ${tumorBAM} ${ref} ${pvcf} && tar cvzf results_${samplename}.tgz results_${samplename}/
    >>>

    output{
        File image_tarball = "results_${samplename}.tgz"
    }
}

workflow svvizflow{
    File pvcf
    File ref
    File tumorBAM
    File tumorBAMIndex
    File controlBAM
    File controlBAMIndex
    String samplename

    call VIZTASK{
        input:
            pvcf=pvcf,
            ref=ref,
            tumorBAM=tumorBAM,
            controlBAM=controlBAM,
            tumorBAMIndex=tumorBAMIndex,
            controlBAMIndex=controlBAMIndex,
            samplename=samplename
    }
}