task Manta {
    String sample_name
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index

    String manta_docker
    String config_manta

    Int? disk_size
    Int? mem_size
    Int? cpu_num
    Int? preemptible_attempts

    command {
        # sample links
        ln -vs ${input_bam} ${sample_name}.bam
        ln -vs ${input_bam_index} ${sample_name}.bai

        # reference links
        ln -vs ${ref_fasta} reference.fasta
        ln -vs ${ref_fasta_index} reference.fasta.fai

        ${config_manta} --bam ${sample_name}.bam \
                        --referenceFasta reference.fasta \
                        --runDir .

        ./runWorkflow.py --mode local \
                         --jobs ${default=32 cpu_num} \
                         --memGb ${default=100 mem_size}

        # change the default names with sample prefix
        mv results/variants/diploidSV.vcf.gz ${sample_name}.diploidSV.vcf.gz
        mv results/variants/diploidSV.vcf.gz.tbi ${sample_name}.diploidSV.vcf.gz.tbi
        mv results/variants/candidateSV.vcf.gz ${sample_name}.candidateSV.vcf.gz
        mv results/variants/candidateSV.vcf.gz.tbi ${sample_name}.candidateSV.vcf.gz.tbi
        mv results/variants/candidateSmallIndels.vcf.gz ${sample_name}.candidateSmallIndels.vcf.gz
        mv results/variants/candidateSmallIndels.vcf.gz.tbi ${sample_name}.candidateSmallIndels.vcf.gz.tbi
    }
    runtime {
        docker: "${manta_docker}"
        memory: select_first([mem_size, 100]) + " GB"
        cpu: select_first([cpu_num, 32])
        disks: "local-disk " + select_first([disk_size, 200]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }
    output {
        File output_sv_vcf = "${sample_name}.diploidSV.vcf.gz"
        File output_sv_vcf_index = "${sample_name}.diploidSV.vcf.gz.tbi"
        File candidate_sv_vcf = "${sample_name}.candidateSV.vcf.gz"
        File candidate_sv_vcf_index = "${sample_name}.candidateSV.vcf.gz.tbi"
        File candidate_indel_vcf = "${sample_name}.candidateSmallIndels.vcf.gz"
        File candidate_indel_vcf_index = "${sample_name}.candidateSmallIndels.vcf.gz.tbi"
    }
}

workflow MantaGermlineSV {
    String sample_name
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index

    String manta_docker
    String config_manta

    call Manta {
        input: sample_name = sample_name,
               input_bam = input_bam,
               input_bam_index = input_bam_index,
               ref_fasta = ref_fasta,
               ref_fasta_index = ref_fasta_index,
               manta_docker = manta_docker,
               config_manta = config_manta
    }
    output {
        Manta.*
    }
}
