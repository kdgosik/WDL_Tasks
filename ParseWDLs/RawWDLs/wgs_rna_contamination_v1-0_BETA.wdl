task wgs_rna_contamination {

    File bam_file
    File bam_index
    File top_genes_bed
    File annotation_gtf
    File hisat_genome_index
    File genome_fasta
    File genome_fasta_index
    String sample_id

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        touch ${bam_index}
        touch ${genome_fasta_index}
        mkdir hisat_index
        tar -zxf ${hisat_genome_index} -C hisat_index --strip-components=1
        /src/rnacont.sh ${bam_file} ${sample_id} ${top_genes_bed} ${annotation_gtf} hisat_index ${genome_fasta}
        tail -1 ${sample_id}.rnaContam.tsv | cut -f9 > pct_contam.txt
    }

    output {
        File rna_contam = "${sample_id}.rnaContam.tsv"
        File split_reads = "${sample_id}.rnacont.tar.gz"
        String pct_rna_contam = read_string("pct_contam.txt")
    }

    runtime {
        docker: "broadinstitute/gtex_genotype_qc"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Xiao Li"
    }
}

workflow wgs_rna_contamination_workflow {
    call wgs_rna_contamination
}
