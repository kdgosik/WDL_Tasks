task runWhamg{
    File bamFile
    File bamIndex
    File reference
    File referenceIndex
    Int threads
    String sampleName

    String fin_str = "${sampleName}"
    
    runtime{
        docker : "erictdawson/svdocker"
        cpu : "${threads}"
        memory : "55 GB"
        disks : "local-disk 1000 HDD"
    }

    command {
        whamg -f ${bamFile} -a ${reference} -x ${threads} > ${fin_str}.wham.vcf
    }

    output {
        File calls = "${fin_str}.wham.vcf"
    }
}

workflow whamFULL{
    File bamFile
    File bamIndex
    File reference
    File referenceIndex
    String sampleName
    Int threads

    call runWhamg{
        input :
            bamFile = bamFile,
            bamIndex = bamIndex,
            sampleName = sampleName,
            reference = reference,
            referenceIndex = referenceIndex,
            threads = threads
    }
}