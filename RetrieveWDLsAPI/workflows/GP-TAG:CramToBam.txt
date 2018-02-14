task CramToBam {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File cram_file
    String sample_id

    String docker_image
    Int? mem_size
    Int? disk_size
    Int? preemptible_attempts

    command <<<
        set -e
        set -o pipefail

        ln -vs ${ref_fasta} reference.fasta
        ln -vs ${ref_fasta_index} reference.fasta.fai
        ln -vs ${ref_dict} reference.dict

        samtools view -h -T reference.fasta ${cram_file} |
        samtools view -b -o ${sample_id}.bam -
        samtools index -b ${sample_id}.bam
        mv ${sample_id}.bam.bai ${sample_id}.bai
    >>>
    runtime {
        docker: "${docker_image}"
        memory: select_first([mem_size, 4]) + " GB"
        disks: "local-disk " + select_first([disk_size, 200]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }
    output {
        File output_bam = "${sample_id}.bam"
        File output_bam_index = "${sample_id}.bai"
    }
}


task ValidateBam {
    File bam_file
    String sample_id
    String picard_jar

    String docker_image
    Int? mem_size
    Int? disk_size
    Int? preemptible_attempts

    command {
        java -Xmx${default=3 mem_size}g -jar ${picard_jar} \
            ValidateSamFile \
            INPUT=${bam_file} \
            OUTPUT="${sample_id}_validation_output.txt" \
            MODE=SUMMARY \
            IS_BISULFITE_SEQUENCED=false
    }
    runtime {
        docker: "${docker_image}"
        memory: select_first([mem_size, 4]) + " GB"
        disks: "local-disk " + select_first([disk_size, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }
    output {
        File validation_output = "${sample_id}_validation_output.txt"
    }
}


workflow CramToBamWorkflow {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File cram_file
    String sample_id

    String docker_image
    String picard_jar

    call CramToBam {
        input: ref_fasta = ref_fasta,
               ref_fasta_index = ref_fasta_index,
               ref_dict = ref_dict,
               cram_file = cram_file,
               sample_id = sample_id,
               docker_image = docker_image
    }

    call ValidateBam {
        input: bam_file = CramToBam.output_bam,
               sample_id = sample_id,
               picard_jar = picard_jar,
               docker_image = docker_image
    }

    output {
        File bam_file = CramToBam.output_bam
        File bam_index_file = CramToBam.output_bam_index
        File validation_output = ValidateBam.validation_output
    }
}
