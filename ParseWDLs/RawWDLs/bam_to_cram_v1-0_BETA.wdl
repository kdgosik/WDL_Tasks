task bam_to_cram {

    File bam_file
    File reference_fasta
    File reference_fasta_index

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Converting BAM to CRAM")
        prefix=$(basename ${bam_file} ".bam")
        cram_file=$prefix".cram"
        samtools view -T ${reference_fasta} -C -o $cram_file ${bam_file}
        samtools index $cram_file
        echo $(date +"[%b %d %H:%M:%S] done")
    }

    output {
        File cram_file=glob("*.cram")[0]
        File cram_index=glob("*.cram.crai")[0]
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow bam_to_cram_workflow {
    call bam_to_cram
}
