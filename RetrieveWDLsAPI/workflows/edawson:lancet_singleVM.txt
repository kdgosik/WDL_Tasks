task runLancet{
    Int threads
    File tumorBam
    File tumorIndex
    File normalBam
    File normalIndex
    String sampleName
    File region_file
    File ref

    command {
        ## Try to bail on repeat areas by ignoring regions w/ greater than 500X coverage
        lancet --tumor ${tumorBam} --normal ${normalBam} --ref ${ref} --bed ${region_file} --num-threads ${threads} --max-avg-cov 500 > ${sampleName}.lancet.vcf
    }

    runtime {
        docker : "erictdawson/svdocker"
        cpu : "${threads}"
        memory : "118 GB"
        disks : "local-disk 1000 HDD"
    }

    output{
        File calls = "${sampleName}.lancet.vcf"
    }
}

workflow lancetFULL{
    File tumor
    File tumorIND
    File normal
    File normalIND
    Int threads
    String name
    File regionFile
    File reference


    call runLancet{
            input : sampleName = name,
                threads=threads,
                region_file = regionFile,
                tumorBam = tumor,
                tumorIndex = tumorIND,
                normalBam = normal,
                normalIndex= normalIND,
                ref = reference
        }

}