task Picard {
    File input_bam
    String bam_index_prefix
    
    command {
        java -jar /usr/gitc/picard.jar BuildBamIndex I=${input_bam} O=${bam_index_prefix}.bai
    }
    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.2.5-1486412288"
        memory: "4G"
        disks: "local-disk 500 HDD"
    }
    output {
        File bam_index_file = "${bam_index_prefix}.bai"
    }
}

workflow BuildBamIndex {
    File input_bam
    String bam_index_prefix
    call Picard {
        input: input_bam = input_bam,
               bam_index_prefix = bam_index_prefix
    }
}