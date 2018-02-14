workflow GenotypeConcordanceWF {
	File call_vcf
    File call_index
    String call_sample
    File truth_vcf
    File truth_index
    String truth_sample
    String output_name
    Array[File] intervals
    File jar
    Int disk_size

	call GenotypeConcordance{
    	input:
        call_vcf=call_vcf,
        call_index=call_index,
        call_sample=call_sample,
        truth_vcf=truth_vcf,
        truth_index=truth_index,
        truth_sample=truth_sample,
        output_name=output_name,
        intervals=intervals,
        jar=jar,
        disk_size=disk_size
    }
}

task GenotypeConcordance {
    File call_vcf
    File call_index
    String call_sample
    File truth_vcf
    File truth_index
    String truth_sample
    String output_name
    Array[File] intervals
    File jar
    Int disk_size

    command {
        java -jar ${jar} GenotypeConcordance CALL_VCF=${call_vcf} CALL_SAMPLE=${call_sample} TRUTH_VCF=${truth_vcf} \
        TRUTH_SAMPLE=${truth_sample} O=${output_name} INTERVALS=${sep=' INTERVALS=' intervals} INTERSECT_INTERVALS=true OUTPUT_VCF=true
    }

    runtime {
         docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:7bc64948a0a9f50ea55edb8b30c710943e44bd861c46a229feaf121d345e68ed"
         memory: "3 GB"
         cpu: "1"
         disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File summary_metrics = "${output_name}.genotype_concordance_summary_metrics"
        File detail_metrics = "${output_name}.genotype_concordance_detail_metrics"
        File contingency_metrics = "${output_name}.genotype_concordance_contingency_metrics"
        File output_vcf = "${output_name}.genotype_concordance.vcf.gz"
        File output_vcf_index = "${output_name}.genotype_concordance.vcf.gz.tbi"
    }
}