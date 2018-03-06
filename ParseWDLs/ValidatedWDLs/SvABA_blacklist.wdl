task svabaCall{
    File tumorBAM
    File tumorIndex
    File normalBAM
    File normalIndex
    File reference
    File refFAIndex
    File refBWTIndex
    File refSAIndex
    File refANNIndex
    File refAMBIndex
    File refPACIndex
    String id
    Int threads
    File redacted_regions
    File dbSNPVCF

    runtime{
        docker : "erictdawson/svdocker:latest"
        memory : "64 GB"
        cpu : "${threads}"
        disks : "local-disk 1000 HDD"
    }

    command{
        svaba run -p ${threads} -t ${tumorBAM} -n ${normalBAM} -G ${reference} -B ${redacted_regions} -a ${id} -D ${dbSNPVCF}
    }

    output{
        File SvABA_Somatic_Indel_VCF = "${id}.svaba.somatic.indel.vcf"
        File SvABA_Somatic_SV_VCF = "${id}.svaba.somatic.sv.vcf"
        File SvABA_Somatic_Unfiltered_indel_VCF = "${id}.svaba.unfiltered.somatic.indel.vcf"
        File SvABA_Somatic_Unfiltered_SV_VCF = "${id}.svaba.unfiltered.somatic.sv.vcf"
        File SvABA_Germline_Indel_VCF = "${id}.svaba.germline.indel.vcf"
        File SvABA_Germline_SV_VCF = "${id}.svaba.germline.sv.vcf"

    }

}

workflow svabaSomatic{
    File tumorBAM
    File tumorIndex

    File normalBAM
    File normalIndex

    File reference
    File refFAIndex
    File refBWTIndex
    File refSAIndex
    File refANNIndex
    File refAMBIndex
    File refPACIndex

    String id
    Int threads
    File redacted_regions
    File dbSNPVCF

    call svabaCall{
        input:
            tumorBAM=tumorBAM,
            tumorIndex=tumorIndex,
            normalBAM=normalBAM,
            normalIndex=normalIndex,
            reference=reference,
            refFAIndex=refFAIndex,
            refBWTIndex=refBWTIndex,
            refSAIndex=refSAIndex,
            refANNIndex=refANNIndex,
            refAMBIndex=refAMBIndex,
            refPACIndex=refPACIndex,
            id=id,
            threads=threads,
            redacted_regions=redacted_regions,
            dbSNPVCF=dbSNPVCF
    }

}