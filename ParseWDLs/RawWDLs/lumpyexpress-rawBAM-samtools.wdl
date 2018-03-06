
task getDiscordants{
    File inputBam
    Int threads
    command{
        samtools view -h -@ ${threads} -F 1294 -u -b -h  ${inputBam} > ${inputBam}.temp && \
        samtools sort -@ ${threads} -m 3G ${inputBam}.temp > discords.bam && rm ${inputBam}.temp
    }
    runtime{
        docker : "erictdawson/lumpy-sv"
    }
    output {
        File discordsBam="discords.bam"
    }
}

task getSplits{
    File bamToSplits
    Int threads
    command {
        samtools view -h -@ ${threads} ${bamToSplits} | \
        /app/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | \
        samtools view -@ ${threads} -b -u - > ${bamToSplits}.temp && \
        samtools sort -@ ${threads} -m 3G ${bamToSplits}.temp > splits.bam && rm ${bamToSplits}.temp
    }
    runtime{
        docker : "erictdawson/lumpy-sv"
    }
    output {
        File splitsBam="splits.bam"
    }
}

task lumpyexpress{
    File inputBam
    File bamSplits
    File bamDiscords
    Int threads

    command {
        lumpyexpress -B ${inputBam} -t ${threads} -S ${bamSplits} -D ${bamDiscords} -o calls.vcf
    }
    runtime {
        docker : "erictdawson/lumpy-sv"
    }
    output{
        File outVCF="calls.vcf"
    }
}

workflow lumpyexpressFULL {
    File inputBam
    Int threads
   
    call getDiscordants{
        input:
            inputBam=inputBam,
            threads=threads
    }

    call getSplits{
        input:
            bamToSplits=inputBam,
            threads=threads
    }

    call lumpyexpress{
        input:
            inputBam=inputBam,
            threads=threads,
            bamSplits=getSplits.splitsBam,
            bamDiscords=getDiscordants.discordsBam
    }
}
