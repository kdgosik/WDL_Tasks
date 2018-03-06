task parseBatchInfo {
    File batchInfo

    command {
        cut -f1 ${batchInfo} > "batch.list"
        cut -f1,2 ${batchInfo} > "batch_md_path.map"
        cut -f3 ${batchInfo} > "del_output_path.list"
        cut -f1,3 ${batchInfo} > "batch_del_output_path.map"
    }

    output {
        Array[String] batchList = read_lines("batch.list")
        Map[String, String] batchMdPathMap = read_map("batch_md_path.map")
        Array[String] delOutputPathList = read_lines("del_output_path.list")
        Map[String, String] batchDelOutputPathMap = read_map("batch_del_output_path.map")
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
    }
}

task extractBatchCalls {
    File delOutputPath
    File referenceBundle

    Int diskSize
    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/del/extract_batch_calls.sh ${delOutputPath} ${referenceBundle}
    }

    output {
        File selectedCallsFile = "SelectedCalls.dat"
    }

    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task mergeDelCalls {

    Array[File] selectedCallsFileList
    File referenceBundle
    
    Int diskSize
    Int numPreempt

    command {
        echo ${sep=' ' selectedCallsFileList} | sed 's/ /\n/g' > selected_calls_file.list
        $SV_DIR/scripts/firecloud/del/merge_del_calls.sh selected_calls_file.list ${referenceBundle}
    }

    output {
        File mergedSitesVcf = "gs_dels.sites.vcf.gz"
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
    String batch
    File partitionFile
    String partitionName
    File mdPath
    File delOutputPath
    File referenceBundle
    File credentialsKeyFile

    Int memory
    Int diskSize
    Int numPreempt
   
    String bamFileList = "del_output/gs_dels.pairs.bam"

    Int adjustedDiskSize = diskSize + round(size(mdPath, "G")) + round(2.2 * size(delOutputPath, "G"))

    command {
        cat ${partitionFile} | awk -v partitionName=${partitionName} '$1 == partitionName' | cut -f 2 > partition.arg

        tar -xvzf ${delOutputPath} || exit 1
        rm ${delOutputPath}

        $SV_DIR/scripts/firecloud/genotyping/run_parallel_genotyper.sh ${vcfFile} ${partitionName} "$(cat partition.arg)" ${mdPath} ${bamFileList} ${referenceBundle} false ${credentialsKeyFile}
    }

    output {
        String batchList = "${batch}"
        File partitionVcf = "gtrun/${partitionName}.genotypes.vcf.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        memory: "${memory}GB"
        disks: "local-disk ${adjustedDiskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task computeBatchVcfList {
    String batch
    Array[String] batchList
    Array[String] partitionVcfList

    Int numPreempt

    command {
        echo ${sep=' ' batchList} | sed 's/ /\n/g' > batch.list
        echo ${sep=' ' partitionVcfList} | sed 's/ /\n/g' > vcf.list
        paste batch.list vcf.list | awk -v batch=${batch} '$1 == batch' | cut -f 2 > batch_vcf.list
    }

    output {
        Array[String] batchVcfList = read_lines("batch_vcf.list")
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        preemptible: "${numPreempt}"
    }
}

task mergeBatchPartitions {
    Array[File] partitionVcfList
    File referenceBundle

    Int diskSize
    Int numPreempt
    
    command {
        echo ${sep=' ' partitionVcfList} | sed 's/ /\n/g' > vcf.list

        source $SV_DIR/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

        java -cp $SV_CLASSPATH -Xmx4g \
            org.broadinstitute.sv.apps.VCFMerge \
            -R $referenceFile \
            -vcf vcf.list \
            -includeInfoTag END \
            -includeInfoTag GSELENGTH \
            -includeInfoTag SVTYPE \
            -O gs_del.genotypes.vcf.gz \
            || exit 1
    }

    output {
        File batchVcf = "gs_del.genotypes.vcf.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task mergeBatches {
    Array[File] batchVcfList
    File referenceBundle

    Int cpu
    Int memory
    Int diskSize
    
    command {
        echo ${sep=' ' batchVcfList} | sed 's/ /\n/g' > vcf.list

        source $SV_DIR/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

        java -cp $SV_CLASSPATH -Xmx4g \
            org.broadinstitute.sv.apps.VCFMerge \
            -R $referenceFile \
            -vcf vcf.list \
            -includeInfoTag END \
            -includeInfoTag GSELENGTH \
            -includeInfoTag SVTYPE \
            -O gs_del.genotypes.vcf.gz \
            || exit 1
    }

    output {
        File mergedGenotypesVcf = "gs_del.genotypes.vcf.gz"
        File mergedGenotypesVcfIndex = "gs_del.genotypes.vcf.gz.tbi"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        cpu: "${cpu}"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
    }
}

task setupRedundancyFiltering {

    File vcfFile
    File referenceBundle

    Int windowSize

    Int diskSize
    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/common/create_genome_partitions.sh ${vcfFile} ${windowSize} intervals.list ${referenceBundle}
    }

    output {
        Array[String] intervalsList = read_lines("intervals.list")
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task filterRedundantSites {
    File vcfFile
    File vcfFileIndex
    String interval
    File referenceBundle

    Int cpu
    Int memory
    Int diskSize
    Int numPreempt

    Int adjustedDiskSize = diskSize + round(size(vcfFile, "G"))
    
    command {
        source $SV_DIR/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

        java -cp $SV_CLASSPATH -Xmx4g \
            org.broadinstitute.sv.apps.FilterRedundantSites \
            -R $referenceFile \
            -vcf ${vcfFile} \
            -L ${interval} \
            -O gs_del.genotypes.vcf.gz \
            || exit 1
    }

    output {
        File dedupedVcf = "gs_del.genotypes.vcf.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        cpu: "${cpu}"
        memory: "${memory}GB"
        disks: "local-disk ${adjustedDiskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task createFinalCallset {
    Array[File] vcfFileList
    File referenceBundle

    Int diskSize

    command {
        echo ${sep=' ' vcfFileList} | sed 's/ /\n/g' > vcf_files.list
        $SV_DIR/scripts/firecloud/del/create_final_callset.sh vcf_files.list ${referenceBundle}
    }

    output {
        File delCallset = "del_callset.tar.gz"
    }

    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
    }
}

workflow gs_create_del_callset_wf {
    File batchInfo
    Int genotypingParallelRecords
    File referenceBundle
    File credentialsKeyFile

    Int memory
    Int diskSize
    Int numThreads
    Int numPreempt

    call parseBatchInfo {
        input:
            batchInfo = batchInfo
    }

    scatter(batch in parseBatchInfo.batchList) {
        call extractBatchCalls {
            input:
                delOutputPath = parseBatchInfo.batchDelOutputPathMap[batch],
                referenceBundle = referenceBundle,
                diskSize = diskSize,
                numPreempt = numPreempt
        }
    }

    call mergeDelCalls {
        input:
            selectedCallsFileList = extractBatchCalls.selectedCallsFile,
            referenceBundle = referenceBundle,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    call setupGenotyping {
        input:
            vcfFile = mergeDelCalls.mergedSitesVcf,
            parallelRecords = genotypingParallelRecords,
            numPreempt = numPreempt
    }

    Array[Pair[String, String]] batchPartitionPairs = cross(parseBatchInfo.batchList, setupGenotyping.partitionList)
    scatter(pair in batchPartitionPairs) {
        String batch = pair.left
        call runParallelGenotyper {
            input: 
                vcfFile = mergeDelCalls.mergedSitesVcf,
                batch = batch,
                partitionFile = setupGenotyping.partitionFile,
                partitionName = pair.right,
                mdPath = parseBatchInfo.batchMdPathMap[batch],
                delOutputPath = parseBatchInfo.batchDelOutputPathMap[batch],
                referenceBundle = referenceBundle,
                credentialsKeyFile = credentialsKeyFile,
                memory = memory,
                diskSize = diskSize,
                numPreempt = numPreempt
        }
    }

    scatter(batch in parseBatchInfo.batchList) {
        call computeBatchVcfList {
            input:
                batch = batch,
                batchList = runParallelGenotyper.batchList,
                partitionVcfList = runParallelGenotyper.partitionVcf,
                numPreempt = numPreempt
        }
    }
    
    scatter(batchVcfList in computeBatchVcfList.batchVcfList) {
        call mergeBatchPartitions {
            input:
                partitionVcfList = batchVcfList,
                referenceBundle = referenceBundle,
                diskSize = diskSize + 2,
                numPreempt = numPreempt
        }
    }

    call mergeBatches {
        input:
            batchVcfList = mergeBatchPartitions.batchVcf,
            referenceBundle = referenceBundle
    }

    call setupRedundancyFiltering {
        input:
            vcfFile = mergeDelCalls.mergedSitesVcf,
            referenceBundle = referenceBundle,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    scatter(interval in setupRedundancyFiltering.intervalsList) {
        call filterRedundantSites {
            input:
                vcfFile = mergeBatches.mergedGenotypesVcf,
                vcfFileIndex = mergeBatches.mergedGenotypesVcfIndex,
                interval = interval,
                referenceBundle = referenceBundle,
                diskSize = diskSize,
                numPreempt = numPreempt
        }
    }

    call createFinalCallset {
        input:
            vcfFileList = filterRedundantSites.dedupedVcf,
            referenceBundle = referenceBundle
    } 

    output {
        File delCallset = createFinalCallset.delCallset
    }
}
