task samtools_count_reads {

    File bam_file

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] samtools: counting reads")
        samtools view -c ${bam_file} > num_reads.txt
        echo $(date +"[%b %d %H:%M:%S] done.")
    }

    output {
        String num_reads = read_string("num_reads.txt")
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


workflow samtools_count_reads_workflow {
    call samtools_count_reads
}
