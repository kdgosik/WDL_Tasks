task MantaTumorNormal {
    String sample_name
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index
    File reference_fasta
    File reference_fasta_index
    String sequencing_type
    String extra_args

    command {
        /usr/manta-1.1.1/bin/configManta.py --tumorBam ${tumor_bam} \
                                            --normalBam ${normal_bam} \
                                            --referenceFasta ${reference_fasta} \
                                            ${extra_args} \
                                            --runDir . &&
        ./runWorkflow.py --mode local \
                         --jobs 32 \
                         --memGb 96 &&

        mv results/variants/somaticSV.vcf.gz ${sample_name}.${sequencing_type}.somaticSV.vcf.gz &&
        mv results/variants/somaticSV.vcf.gz.tbi ${sample_name}.${sequencing_type}.somaticSV.vcf.gz.tbi
    }
    runtime {
        docker: "jnktsj/manta-1.1.1:0.1.0"
        memory: "100 GB"
        cpu: "32"
        disks: "local-disk 500 HDD"
    }
    output {
        File output_vcf = "${sample_name}.${sequencing_type}.somaticSV.vcf.gz"
        File output_vcf_index = "${sample_name}.${sequencing_type}.somaticSV.vcf.gz.tbi"
    }
}


task ConvertToBedTabix {
    File? interval
    String output_bed = "./interval.bed"

    command <<<
        if [[ "${interval}" == *.interval_list ]]
        then
            # Convert Picard-style intervals to BED
            grep -v "^@" ${interval} | awk '{$2-=1; print $1,$2,$3,$5}' OFS="\t" > ${output_bed}
        else
            cp ${interval} ${output_bed}
        fi

        ./bgzip ${output_bed}
        ./tabix -p bed ${output_bed}.gz
    >>>
    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.2.5-1486412288"
    }
    output {
        File out_interval = "${output_bed}.gz"
        File out_interval_index = "${output_bed}.gz.tbi"
    }
}


workflow MantaTumorNormalWorkflow {
    String sample_name
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index
    File reference_fasta
    File reference_fasta_index
    File? interval

    Boolean is_exome = defined(interval)

    if (is_exome) {
        call ConvertToBedTabix {
            input: interval = interval
        }
        call MantaTumorNormal as MantaTumorNormalWES {
            input: sample_name = sample_name,
                   tumor_bam = tumor_bam,
                   tumor_bam_index = tumor_bam_index,
                   normal_bam = normal_bam,
                   normal_bam_index = normal_bam_index,
                   reference_fasta = reference_fasta,
                   reference_fasta_index = reference_fasta_index,
                   sequencing_type = "WES",
                   extra_args = "--exome --callRegions " + ConvertToBedTabix.out_interval
        }
        output {
            MantaTumorNormalWES.*
        }
    }
    if (!is_exome) {
        call MantaTumorNormal as MantaTumorNormalWGS {
            input: sample_name = sample_name,
                   tumor_bam = tumor_bam,
                   tumor_bam_index = tumor_bam_index,
                   normal_bam = normal_bam,
                   normal_bam_index = normal_bam_index,
                   reference_fasta = reference_fasta,
                   reference_fasta_index = reference_fasta_index,
                   sequencing_type = "WGS",
                   extra_args = ""
        }
        output {
            MantaTumorNormalWGS.*
        }
    }
}