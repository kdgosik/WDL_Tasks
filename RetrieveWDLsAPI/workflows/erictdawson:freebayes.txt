task pFBTask{
        File bamFile
        File bamIndex
        Int threads
        String id
        File reference
        File referenceIndex

# cd /app/freebayes/scripts/ && ./freebayes-parallel <(python fasta_generate_regions.py ${referenceIndex} 100000) ${threads} -f ${reference} ${bamFile} > ${id}.freebayes.vcf
        String dollar = "$"
        command <<<
            
            if [ -e jobfile.txt ]
                then rm -f jobfile.txt
            fi
            
            for i in `python /app/freebayes/scripts/fasta_generate_regions.py ${referenceIndex} 100000`;
            do
                echo "freebayes -f ${reference} --region ${dollar}{i} ${bamFile} > ${id}.${dollar}{i}.fb.vcf" >> jobfile.txt
            done
            python /app/LaunChair/launcher.py -i jobfile.txt -c 1 -n ${threads} && \
            cat ${id}.*.fb.vcf | /usr/bin/vcffirstheader | /usr/bin/vcfstreamsort -w 1000 | /usr/bin/vcfuniq > ${id}.freebayes.vcf
        >>>

    runtime {
        docker : "erictdawson/svdocker"
        memory : "28GB"
        cpu : "${threads}"
        disks : "local-disk 1000 HDD"
    }

    output{
        File FBvcf = "${id}.freebayes.vcf"
    }
}

workflow FreeBayes{
    File bamFile
        File bamIndex
        Int threads
        String id
        File reference
        File referenceIndex

        call pFBTask{
            input:
            bamFile=bamFile,
                bamIndex=bamIndex,
                threads=threads,
                id=id,
                reference=reference,
                referenceIndex=referenceIndex
        }

}