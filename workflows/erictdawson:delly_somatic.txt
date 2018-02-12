task dellyCall{
    File tumorBAM
    File normalBAM
    Int threads
    File reference
    File index
    String type
    String sampleName

    command{
        export OMP_NUM_THREADS=${threads} && delly call --type ${type} -g ${reference} -o ${sampleName}.somatic.${type}.bcf ${tumorBAM} ${normalBAM}
    }

    runtime{
        docker : "erictdawson/svdocker"
        memory : "40 GB"
        cpu : "${threads}"
        disks : "local-disk 1000 HDD"
    }

    output{
        File xbcf = "${sampleName}.somatic.${type}.bcf"
    }
}

task vcflibMerge{
    File insBCF
    File invBCF
    File delBCF
    String sampleName

    command {
        bcftools view ${insBCF} > ${sampleName}.delly.somatic.ins.vcf
        bcftools view ${invBCF} > ${sampleName}.delly.somatic.inv.vcf
        bcftools view ${delBCF} > ${sampleName}.delly.somatic.del.vcf
        vcfcombine ${sampleName}.delly.somatic.ins.vcf ${sampleName}.delly.somatic.inv.vcf ${sampleName}.delly.somatic.del.vcf > ${sampleName}.inv.ins.del.delly.somatic.vcf
    }

    runtime {
        docker : "erictdawson/svdocker"
        memory : "16GB"
        cpu : "1"
        disks : "local-disk 1000 HDD"
    }
    output{
        File merged = "${sampleName}.merged.delly.vcf"
    }
}




workflow dellyAll{
    File tumorBAM
    File normalBAM
    File index
    File reference
    Int threads
    String name
    
    call dellyCall as insCall{
        input:
           tumorBAM=tumorBAM,
           normalBAM=normalBAM,
           reference=reference,
           index=index,
           threads=threads,
           type="INS",
           sampleName=name
    }

    call dellyCall as invCall{
        input:
           tumorBAM=tumorBAM,
           normalBAM=normalBAM,
           reference=reference,
           index=index,
           threads=threads,
           type="INV",
           sampleName=name
    }

    call dellyCall as delCall{
        input:
           tumorBAM=tumorBAM,
           normalBAM=normalBAM,
           reference=reference,
           index=index,
           type="DEL",
           threads=threads,
           sampleName=name
    }

    call vcflibMerge{
        input:
            insBCF=insCall.xbcf,
            delBCF=delCall.xbcf,
            invBCF=invCall.xbcf,
            sampleName=name
        }
    }