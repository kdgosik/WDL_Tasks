task gs_build_profiles {
    Int binSize
    Array[File] mdArchives
    File gsReferenceBundle

    Int memory
    Int diskSize
    Int numThreads
    Int numPreempt
 
    command {
        echo ${sep=' ' mdArchives} | sed 's/ /\n/g' > md_path.list
        $SV_DIR/scripts/firecloud/gs_build_profiles.sh md_path.list ${binSize} ${gsReferenceBundle}
    }

    output {
        File profiles = "profiles_${binSize}.tar.gz"
    }

    runtime {
        docker: "skashin/genome-strip:latest"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        cpu: "${numThreads}"
        preemptible: "${numPreempt}"
    }
}

workflow gs_build_profiles_wf {
    call gs_build_profiles
}
