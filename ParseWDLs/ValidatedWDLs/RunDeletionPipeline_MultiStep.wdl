task setupDiscovery {

    String searchIntervals
    File referenceBundle

    Int diskSize
    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/del/compute_partitions.sh "${searchIntervals}" partitions.dat ${referenceBundle} || exit 1
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
    Boolean storeReadPairFile
    String mdPath
    Array[String] bamFileList
    File referenceBundle
    File repeatTrackFile
    File credentialsKeyFile

    Int memory
    Int diskSize
    Int numThreads
    Int numPreempt

    command {
        cat ${partitionFile} | awk -v partitionName=${partitionName} '$1 == partitionName' | cut -f 2- | sed 's/\t/ /g' > partition.args
        echo ${sep=' ' bamFileList} | sed 's/ /\n/g' > bam_files.list

        $SV_DIR/scripts/firecloud/del/run_parallel_discovery.sh ${partitionName} "$(cat partition.args)" ${storeReadPairFile} ${mdPath} bam_files.list ${referenceBundle} ${repeatTrackFile} ${credentialsKeyFile}
    }

    output {
        File partitionOutput = "del_output.tar.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        cpu: "${numThreads}"
        preemptible: "${numPreempt}"
    }
}

task mergeDiscovery {
    String filePrefix
    Array[File] partitionOutputList
    Boolean storeReadPairFile
    File referenceBundle
    File repeatTrackFile

    Int diskSize
    Int numPreempt

    command {
        echo ${sep=' ' partitionOutputList} | sed 's/ /\n/g' > partition_archives.list

        $SV_DIR/scripts/firecloud/del/merge_parallel_discovery.sh ${filePrefix} partition_archives.list ${storeReadPairFile} ${referenceBundle} ${repeatTrackFile}
    }

    output {
        File delOutput = "del_output.tar.gz"
        File genotypingSitesVcf = "del_output/genotyping/gs_dels.sites.vcf.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task setupGenotyping {

    File vcfFile
    Int parallelRecords

    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/compute_vcf_partitions.sh ${vcfFile} ${parallelRecords} partitions.dat
        cut -f 1 partitions.dat > partitions.list
    }

    output {
        File partitionFile = "partitions.dat"
        Array[String] partitionList = read_lines("partitions.list")
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        preemptible: "${numPreempt}"
    }
}

task runParallelGenotyper {
    File vcfFile
    File partitionFile
    String partitionName
    String mdPath
    File delOutput
    #Array[String] bamFileList
    File referenceBundle
    File credentialsKeyFile

    Int memory
    Int diskSize
    Int numPreempt

    String bamFileList = "del_output/gs_dels.pairs.bam"

    command {
        cat ${partitionFile} | awk -v partitionName=${partitionName} '$1 == partitionName' | cut -f 2 > partition.arg
        #echo ${sep=' ' bamFileList} | sed 's/ /\n/g' > bam_files.list

        tar -xvzf ${delOutput} || exit 1

        $SV_DIR/scripts/firecloud/genotyping/run_parallel_genotyper.sh ${vcfFile} ${partitionName} "$(cat partition.arg)" ${mdPath} ${bamFileList} ${referenceBundle} false ${credentialsKeyFile}
    }

    output {
        File partitionOutput = "gtrun.tar.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task mergeGenotyping {
    String filePrefix
    Array[File] partitionOutputList
    File referenceBundle

    Int diskSize
    Int numPreempt

    command {
        echo ${sep=' ' partitionOutputList} | sed 's/ /\n/g' > partition_archives.list

        $SV_DIR/scripts/firecloud/genotyping/merge_parallel_genotyping.sh ${filePrefix} partition_archives.list ${referenceBundle}
    }

    output {
        File gtRun = "gtrun.tar.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task mergeResults {
    File gtRun
    File delOutput

    Int diskSize
    Int numPreempt

    command {
        tar -xvzf ${delOutput} || exit 1
        tar -xvzf ${gtRun} -C del_output/genotyping/ --strip-components 1 || exit 1
        tar -cvzf del_output.tar.gz del_output
    }

    output {
        File mergedOutput = "del_output.tar.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

workflow gs_run_del_pipeline_multi_wf {

    String mdPath
    Array[String] bamFileList
    Boolean storeReadPairFile
    String searchIntervals
    File referenceBundle
    File repeatTrackFile
    File credentialsKeyFile

    Int memory
    Int diskSize
    Int numThreads
    Int numPreempt

    String filePrefix = "gs_dels"

    call setupDiscovery {
        input:
            searchIntervals = searchIntervals,
            referenceBundle = referenceBundle,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    scatter(partitionName in setupDiscovery.partitionList) {
        call runParallelDiscovery {
            input: 
                partitionFile = setupDiscovery.partitionFile,
                partitionName = partitionName,
                storeReadPairFile = storeReadPairFile,
                mdPath = mdPath,
                bamFileList = bamFileList,
                referenceBundle = referenceBundle,
                repeatTrackFile = repeatTrackFile,
                credentialsKeyFile = credentialsKeyFile,
                memory = memory,
                diskSize = diskSize + 15,
                numThreads = numThreads,
                numPreempt = numPreempt
        }
    }

    call mergeDiscovery {
        input:
            filePrefix = filePrefix,
            partitionOutputList = runParallelDiscovery.partitionOutput,
            storeReadPairFile = storeReadPairFile,
            referenceBundle = referenceBundle,
            repeatTrackFile = repeatTrackFile,
            diskSize = diskSize + 50,
            numPreempt = numPreempt
    }

    call setupGenotyping {
        input:
            vcfFile = mergeDiscovery.genotypingSitesVcf,
            parallelRecords = 100,
            numPreempt = numPreempt
    }

    scatter(partitionName in setupGenotyping.partitionList) {
        call runParallelGenotyper {
            input: 
                vcfFile = mergeDiscovery.genotypingSitesVcf,
                partitionFile = setupGenotyping.partitionFile,
                partitionName = partitionName,
                mdPath = mdPath,
                delOutput = mergeDiscovery.delOutput,
                #bamFileList = bamFileList,
                referenceBundle = referenceBundle,
                credentialsKeyFile = credentialsKeyFile,
                memory = memory,
                diskSize = diskSize + 30,
                numPreempt = numPreempt
        }
    }

    call mergeGenotyping {
        input:
            filePrefix = filePrefix,
            partitionOutputList = runParallelGenotyper.partitionOutput,
            referenceBundle = referenceBundle,
            diskSize = diskSize + 20,
            numPreempt = numPreempt
    }

    call mergeResults {
        input:
            gtRun = mergeGenotyping.gtRun,
            delOutput = mergeDiscovery.delOutput,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    output {
        File delOutput = mergeResults.mergedOutput
    }
}