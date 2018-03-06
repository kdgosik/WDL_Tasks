task getDiscordants{
    File inputBam
    File inputBamIndex
    String sample
    String tag
    Int threads

    String fin_sample_str = sub(sample, "-", "_")

        #sambamba view -h -f bam --num-filter /1294 -o ${sample}.${tag}.bambatmp.bam ${inputBam}  &&  sambamba sort -t ${threads} -o ${sample}.${tag}.discords.bam ${sample}.${tag}.bambatmp.bam

    runtime{
        bootDiskSizeGb: 550
        docker : "erictdawson/svdocker:latest"
	    cpu : "${threads}"
	    memory : "64 GB"
	    disks : "local-disk 1500 HDD"
    }

    command {
        find . && samtools view -h -F 1294 -b -@ ${threads} ${inputBam} > ${sample}.${tag}.discords.unsrt.bam && \
        samtools sort -@ ${threads} -m 2G ${sample}.${tag}.discords.unsrt.bam > ${fin_sample_str}.${tag}.discords.bam && find .
    }

    output {
        File discordsBam="${fin_sample_str}.${tag}.discords.bam"
    }
}

task getSplits{
    File bamToSplits
    File inputBamIndex
    Int threads
    String sample
    String tag

        #samtools view -h -@ ${threads} ${bamToSplits} | \
        #/app/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin > tmp.discord.bam && \
        #sambamba view -S -f bam -h -l 5 -t ${threads} -o /dev/stdout /dev/stdin | \
        #sambamba sort -t ${threads} -o ${sample}.${tag}.splits.bam /dev/stdin; find .
 
    String fin_sample_str = sub(sample, "-", "_")

    runtime{
        bootDiskSizeGb: 550
        docker : "erictdawson/svdocker:latest"
	    cpu : "${threads}"
	    memory : "64 GB"
	    disks : "local-disk 1500 HDD"
    }

    command {
        find . && samtools view -h ${bamToSplits} \
            | /app/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
            | samtools view -b - > ${sample}.${tag}.splits.unsrt.bam && \
            samtools sort -m 2G -@ ${threads} ${sample}.${tag}.splits.unsrt.bam > ${fin_sample_str}.${tag}.splits.bam && find .
          }

    output {
        File splitsBam="${fin_sample_str}.${tag}.splits.bam"
    }
}

task lumpyexpress{
    File tumorBam
    File normalBam
    File tumorBamIndex
    File normalBamIndex
    File tumorSplits
    File tumorDiscords
    File normalSplits
    File normalDiscords
    Int threads
    String sampleName

    String fin_sample_str = sub(sampleName, "-", "_")

    command {
        lumpyexpress -B ${tumorBam},${normalBam} -t ${threads} -S ${tumorSplits},${normalSplits} -D ${tumorDiscords},${normalDiscords} -o ${fin_sample_str}.tn.lumpy.vcf && sleep 1 && find .
    }
    runtime {
        bootDiskSizeGb: 550
        docker : "erictdawson/svdocker:latest"
	    cpu : "${threads}"
	    memory : "120 GB"
	    disks : "local-disk 1500 HDD"
    }
    output {
        File outVCF = "${fin_sample_str}.tn.lumpy.vcf"
    }
}


workflow lumpyexpressFULL {
    File tumorBam
    File normalBam
    File tumorBamIndex
    File normalBamIndex
    Int lumpyThreads
    Int preThreads
    String name

    call getDiscordants as tumorDiscord{
        input:
            inputBam=tumorBam,
            inputBamIndex=tumorBamIndex,
            threads=preThreads,
            sample=name,
            tag="tumor"
    }

    call getDiscordants as normalDiscord{
        input:
            inputBam=normalBam,
            inputBamIndex=normalBamIndex,
            threads=preThreads,
            sample=name,
            tag="normal"
    }

    call getSplits as tumorSplit{
        input:
            bamToSplits=normalBam,
            inputBamIndex=tumorBamIndex,
            threads=preThreads,
            sample=name,
            tag="tumor"
    }

 
    call getSplits as normalSplit{
        input:
            bamToSplits=normalBam,
            inputBamIndex=normalBamIndex,
            threads=preThreads,
            sample=name,
            tag="normal"
    }

   

    call lumpyexpress{
        input:
            tumorBam=tumorBam,
            normalBam=normalBam,
            tumorBamIndex=tumorBamIndex,
            normalBamIndex=normalBamIndex,
            threads=lumpyThreads,
            tumorSplits=tumorSplit.splitsBam,
            normalSplits=normalSplit.splitsBam,
            tumorDiscords=tumorDiscord.discordsBam,
            normalDiscords=normalDiscord.discordsBam,
            sampleName=name
    }

    output {
        File call_vcf = lumpyexpress.outVCF
    }
}