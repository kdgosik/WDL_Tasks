task setupDiscovery {

    Int windowSize
    Int windowOverlap
    String intervals
    File referenceBundle

    Int diskSize
    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/cnv/compute_partitions.sh ${windowSize} ${windowOverlap} "${intervals}" partitions.dat ${referenceBundle} || exit 1
        cut -f 1 partitions.dat > partitions.list
    }

    output {
        File partitionFile = "partitions.dat"
        Array[String] partitionList = read_lines("partitions.list")
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task runParallelDiscovery {
    File partitionFile
    String partitionName
    File mdPath
    File referenceBundle
    File credentialsKeyFile

    Int memory
    Int diskSize
    Int numThreads
    Int numPreempt

    Int adjustedDiskSize = diskSize + round(size(mdPath, "G")) + 10

    command {
        cat ${partitionFile} | awk -v partitionName=${partitionName} '$1 == partitionName' | cut -f 2 > interval.list

        $SV_DIR/scripts/firecloud/cnv/run_parallel_discovery.sh ${partitionName} ${mdPath} interval.list ${referenceBundle} ${credentialsKeyFile}
    }

    output {
        File partitionOutput = "cnv_output.tar.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        memory: "${memory}GB"
        disks: "local-disk ${adjustedDiskSize} HDD"
        cpu: "${numThreads}"
        preemptible: "${numPreempt}"
    }
}

task mergeDiscovery {
    String mdPath
    Array[File] partitionOutputList
    File referenceBundle
    File credentialsKeyFile

    Int memory
    Int numThreads
    Int diskSize
    Int numPreempt

    command {
        echo ${sep=' ' partitionOutputList} | sed 's/ /\n/g' > partition_archives.list

        $SV_DIR/scripts/firecloud/cnv/merge_parallel_discovery.sh ${mdPath} partition_archives.list ${referenceBundle} ${credentialsKeyFile}
    }

    output {
        File cnvOutput = "cnv_output.tar.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
        memory: "${memory}GB"
        cpu: "${numThreads}"
    }
}

workflow gs_run_cnv_pipeline_wf {

    String mdPath
    Int windowSize
    Int windowOverlap
    String intervals
    File referenceBundle
    File credentialsKeyFile

    Int memory
    Int diskSize
    Int numThreads
    Int numPreempt

    String filePrefix = "gs_cnvs"

    call setupDiscovery {
        input:
            windowSize = windowSize,
            windowOverlap = windowOverlap,
            intervals = intervals,
            referenceBundle = referenceBundle,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    scatter(partitionName in setupDiscovery.partitionList) {
        call runParallelDiscovery {
            input: 
                partitionFile = setupDiscovery.partitionFile,
                partitionName = partitionName,
                mdPath = mdPath,
                referenceBundle = referenceBundle,
                credentialsKeyFile = credentialsKeyFile,
                memory = memory,
                diskSize = diskSize,
                numThreads = numThreads,
                numPreempt = numPreempt
        }
    }

    call mergeDiscovery {
        input:
            mdPath = mdPath,
            partitionOutputList = runParallelDiscovery.partitionOutput,
            referenceBundle = referenceBundle,
            credentialsKeyFile = credentialsKeyFile,
            numPreempt = numPreempt
    }

    output {
        File cnvOutput = mergeDiscovery.cnvOutput
    }
}
