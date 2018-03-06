task runLancet{
    Int threads
    File tumorBam
    File tumorIndex
    File normalBam
    File normalIndex
    String sampleName
    String region
    File ref

    String pre_region_str = sub(region, ":", "_")
    String fin_region_str = sub(pre_region_str, "-", "_")

    command {
        ## Try to bail on repeat areas by ignoring regions w/ greater than 1000X coverage
        lancet --tumor ${tumorBam} --normal ${normalBam} --ref ${ref} --reg ${region} --num-threads ${threads} --max-avg-cov 1000 > ${sampleName}.${fin_region_str}.lancet.vcf
    }

    runtime {
        docker : "erictdawson/svdocker"
        cpu : "${threads}"
        memory : "118 GB"
        disks : "local-disk 1000 HDD"
    }

    output{
        
        File calls = "${sampleName}.${fin_region_str}.lancet.vcf"
        #File calls = "${sampleName}.lancet.vcf"
    }
}

task mergeVCF{
    Array[File] subVCFs
    String sampleName

    runtime{
        docker : "erictdawson/svdocker"
        cpu : "1"
        memory : "24 GB"
        disks : "local-disk 250 HDD"
    }

    command {
        vcfcat ${sep=" " subVCFs} > ${sampleName}.lancet.vcf
    }

    output{
        File vcf = "${sampleName}.lancet.vcf"
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

    Array[String] regions = read_lines(regionFile)

    scatter (reg in regions){
        call runLancet{
            input : sampleName = name,
                threads=threads,
                region = reg,
                tumorBam = tumor,
                tumorIndex = tumorIND,
                normalBam = normal,
                normalIndex= normalIND,
                ref = reference
        }
    }

    call mergeVCF{
        input: 
            sampleName = name,
            subVCFs = runLancet.calls
    }
}