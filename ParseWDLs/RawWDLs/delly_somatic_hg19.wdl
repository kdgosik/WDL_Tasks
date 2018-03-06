task dellyCall{
    File tumorBAM
    File normalBAM
    Int threads
    File reference
    File tumorIndex
    File normalIndex
    String type
    String sampleName

    command{
        export OMP_NUM_THREADS=${threads} && delly call --type ${type}  -x /app/delly/excludeTemplates/human.hg19.excl.tsv -g ${reference} -o ${sampleName}.somatic.${type}.bcf ${tumorBAM} ${normalBAM}
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
    File bndBCF
    File dupBCF
    String sampleName

    command {
        bcftools view ${insBCF} > ${sampleName}.delly.somatic.ins.vcf
        bcftools view ${invBCF} > ${sampleName}.delly.somatic.inv.vcf
        bcftools view ${delBCF} > ${sampleName}.delly.somatic.del.vcf
        bcftools view ${bndBCF} > ${sampleName}.delly.somatic.bnd.vcf
        bcftools view ${dupBCF} > ${sampleName}.delly.somatic.dup.vcf
        vcfcombine ${sampleName}.delly.somatic.ins.vcf ${sampleName}.delly.somatic.inv.vcf ${sampleName}.delly.somatic.del.vcf ${sampleName}.delly.somatic.bnd.vcf ${sampleName}.delly.somatic.dup.vcf > ${sampleName}.inv.ins.del.bnd.dup.delly.somatic.vcf
    }

    runtime {
        docker : "erictdawson/svdocker"
        memory : "16GB"
        cpu : "1"
        disks : "local-disk 1000 HDD"
    }
    output{
        File merged = "${sampleName}.inv.ins.del.bnd.dup.delly.somatic.vcf"
    }
}




workflow dellyAll{
    File tumorBAM
    File normalBAM
    File tumorIndex
    File normalIndex
    File reference
    Int threads
    String name
    
    call dellyCall as insCall{
        input:
           tumorBAM=tumorBAM,
           normalBAM=normalBAM,
           reference=reference,
           tumorIndex=tumorIndex,
           normalIndex=normalIndex,
           threads=threads,
           type="INS",
           sampleName=name
    }

    call dellyCall as invCall{
        input:
           tumorBAM=tumorBAM,
           normalBAM=normalBAM,
           reference=reference,
           tumorIndex=tumorIndex,
           normalIndex=normalIndex,
           threads=threads,
           type="INV",
           sampleName=name
    }

    call dellyCall as delCall{
        input:
           tumorBAM=tumorBAM,
           normalBAM=normalBAM,
           reference=reference,
           tumorIndex=tumorIndex,
           normalIndex=normalIndex,
           type="DEL",
           threads=threads,
           sampleName=name
    }

    call dellyCall as bndCall{
        input:
           tumorBAM=tumorBAM,
           normalBAM=normalBAM,
           reference=reference,
           tumorIndex=tumorIndex,
           normalIndex=normalIndex,
           type="BND",
           threads=threads,
           sampleName=name
    }

    call dellyCall as dupCall{
        input:
           tumorBAM=tumorBAM,
           normalBAM=normalBAM,
           reference=reference,
           tumorIndex=tumorIndex,
           normalIndex=normalIndex,
           type="DUP",
           threads=threads,
           sampleName=name
    }



    call vcflibMerge{
        input:
            insBCF=insCall.xbcf,
            delBCF=delCall.xbcf,
            invBCF=invCall.xbcf,
            bndBCF=bndCall.xbcf,
            dupBCF=dupCall.xbcf,
            sampleName=name
        }
    }