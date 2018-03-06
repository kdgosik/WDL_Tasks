import "https://api.firecloud.org/ga4gh/v1/tools/bgranger:GCHaplotypeCaller/versions/2/plain-WDL/descriptor" as hc
import "https://api.firecloud.org/ga4gh/v1/tools/GPTAG:GenotypeConcordance/versions/1/plain-WDL/descriptor" as gc

workflow DownsampleSubworkflow {

	File ref_fasta = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
	File ref_fasta_index = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
	File ref_dict = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.dict"

	File bam
	File bam_index
	Float bam_coverage
	Array[Float] downsample_coverages
	File picard_jar
	String sample_name
	String gc_sample_name
	Int disk_size

	File nist_vcf
	File nist_vcf_index
    Array[File] intervals

	File gatk_jar

	scatter (coverage in downsample_coverages) {

		call DownsampleSam {
			input:
				bam = bam,
				bam_index = bam_index,
				coverage = bam_coverage,
				desired_coverage = coverage,
				jar = picard_jar,
				output_name = sample_name + "_" + coverage,
				disk_size = disk_size
		}

		call MarkDuplicates {
			input:
				bam = DownsampleSam.output_bam,
				bam_index = DownsampleSam.output_bam_index,
				jar = picard_jar,
				output_name = sample_name + "_" + coverage,
				disk_size = disk_size
		}

		call CollectWgsMetrics {
			input:
				bam = MarkDuplicates.output_bam,
				bam_index = MarkDuplicates.output_bam_index,
				ref = ref_fasta,
				ref_index =  ref_fasta_index,
				jar = picard_jar,
				output_name = sample_name + "_" + coverage,
				disk_size = disk_size
		}

		call hc.HaplotypeCallerWf as HaplotypeCaller {
		    input:
                input_bam = MarkDuplicates.output_bam,
                input_bam_index = MarkDuplicates.output_bam_index,
                sample_name = sample_name + "_" + coverage,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                disk_size = disk_size
		}

        call GenotypeGVCFs {
            input:
                gvcf = HaplotypeCaller.gvcf,
                gvcf_index = HaplotypeCaller.gvcf_index,
        	    ref = ref_fasta,
                ref_index = ref_fasta_index,
        	    ref_dict = ref_dict,
                output_name = sample_name + "_" + coverage,
        	    jar = gatk_jar,
                disk_size = disk_size
        }

		call gc.GenotypeConcordance as GenotypeConcordance{
            input:
                call_vcf = GenotypeGVCFs.output_vcf,
                call_index = GenotypeGVCFs.output_vcf_index,
                call_sample = gc_sample_name,
                truth_vcf = nist_vcf,
                truth_index = nist_vcf_index,
                truth_sample = "HG001",
                intervals = intervals,
                output_name = sample_name + "_" + coverage,
                jar = picard_jar,
                disk_size = disk_size
        }
	}

	output {
        Array[File] gc_summary_metrics = GenotypeConcordance.summary_metrics
        Array[File] wgs_metrics = CollectWgsMetrics.metrics
	}
}

task DownsampleSam {
	File bam
	File bam_index
	Float coverage
	Float desired_coverage
	File jar # picard jar
	String output_name
	Int disk_size

	command <<<
		downsampleFraction="$(python -c "print ${desired_coverage}/${coverage}")"
		java -jar ${jar} PositionBasedDownsampleSam INPUT=${bam} OUTPUT=${output_name}.bam F=$downsampleFraction CREATE_INDEX=true
	>>>

	runtime {
         docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:7bc64948a0a9f50ea55edb8b30c710943e44bd861c46a229feaf121d345e68ed"
         memory: "3 GB"
         cpu: "1"
         disks: "local-disk " + disk_size + " HDD"
    }

	output {
		File output_bam = "${output_name}.bam"
		File output_bam_index = "${output_name}.bai"
	}
}

task MarkDuplicates {

	File bam
	File bam_index
	File jar # picard jar
	String output_name
	Int disk_size

	command {
		java -jar ${jar} MarkDuplicates I=${bam} O=${output_name}.bam CREATE_INDEX=true METRICS_FILE=${output_name}.duplication_metrics \
		VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CLEAR_DT="false" ADD_PG_TAG_TO_READS=false
	}

	runtime {
         docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:7bc64948a0a9f50ea55edb8b30c710943e44bd861c46a229feaf121d345e68ed"
         memory: "3 GB"
         cpu: "1"
         disks: "local-disk " + disk_size + " HDD"
    }

	output {
		File output_bam = "${output_name}.bam"
		File output_bam_index = "${output_name}.bai"
		File metrics = "${output_name}.duplication_metrics"
	}
}

task CollectWgsMetrics {
	File bam
	File bam_index
	File ref
	File ref_index
	File jar # picard jar
	String output_name
	Int disk_size

	command {
		java -jar ${jar} CollectWgsMetrics I=${bam} O=${output_name}.wgs_metrics R=${ref}
	}

	runtime {
         docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:7bc64948a0a9f50ea55edb8b30c710943e44bd861c46a229feaf121d345e68ed"
         memory: "3 GB"
         cpu: "1"
         disks: "local-disk " + disk_size + " HDD"
    }

	output {
		File metrics = "${output_name}.wgs_metrics"
	}
}

task GenotypeGVCFs {
	File gvcf
	File gvcf_index
	File ref
	File ref_index
	File ref_dict
    String output_name
	File jar
    Int disk_size

	command {
		java -jar ${jar} -T GenotypeGVCFs -R ${ref} --variant ${gvcf} -o ${output_name}.vcf
	}

	runtime {
         docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:7bc64948a0a9f50ea55edb8b30c710943e44bd861c46a229feaf121d345e68ed"
         memory: "3 GB"
         cpu: "1"
         disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_vcf = "${output_name}.vcf"
        File output_vcf_index = "${output_name}.vcf.idx"
    }
}