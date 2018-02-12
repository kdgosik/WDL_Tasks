workflow runThisFirst {
	String capture_platform
	File ref_fasta
	File ref_fasta_index
	File ref_fasta_dict
	File input_bam
	File input_bam_index
	File normal_bam
	File normal_bam_index
	String sample_name

	# mimicking PicardTargetMapper
	Array[String] ice_files = [
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/CRSP_ICE_hg19_wex_illumina_v1.no_X_Y_MT.bed",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/ICE_Combined_EEW_EWS_v1.pon",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/ICE_406_Normal_samples_PoN.vcf",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/whole_exome_illumina_coding_v1_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/ice_targets.bed",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/ice_rcs_eval.v1.pd250.spark.pon"
	]

	Array[String] agilent_files = [
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/whole_exome_agilent_1.1_refseq_plus_3_boosters_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets_germline_copy_number_variants_X_Y_removed.bed",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/BI_with_FFPE_Agilent_collapsed_v1.pon",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/gaf_20111020+broad_wex_1.1_hg19.bed",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/agilent_targets.bed",
	"gs://fc-d9e360d2-1eb7-49d0-8007-356220e59f1a/param_files/agilent_rcs_eval_FFPE.v1.pd250.spark.pon"
	]

	# choose the appropriate set of reference files
	Array[String] files_to_use = if capture_platform == "ICE" then ice_files else agilent_files

	call depthOfCov {
		input: refFasta=ref_fasta,
			refFastaIndex=ref_fasta_index,
			refFastaDict=ref_fasta_dict,
			inputBam=input_bam,
			inputBamIndex=input_bam_index,
			sampleName=sample_name
	}

	call contestForCaptureArrayFree {
		input: refFasta=ref_fasta,
			targetsIntervalList=files_to_use[0],
			refFastaIndex=ref_fasta_index,
			refFastaDict=ref_fasta_dict,
			inputBam=input_bam,
			inputBamIndex=input_bam_index,
			normalBam=normal_bam,
			normalBamIndex=normal_bam_index,
			sampleName=sample_name
	}

	output {
		# mimicking ExtractBamIDs
		String bamfile_id = basename(input_bam, ".bam")

		# mimicking PicardTargetMapper
		File contest_capture_targets_interval_list = files_to_use[0]
		File ReCapSeg_target_bed = files_to_use[1]
		File ReCapSeg_PON = files_to_use[2]
		File mutect_panel_of_normals_vcf_capture = files_to_use[3]
		File mutect_capture_targets_interval_list = files_to_use[4]
		File gatk4cnv_target_bed_capture = files_to_use[5]
		File gatk4cnv_pon_capture = files_to_use[6]

		# DepthOfCoverage output files
		File sample_gene_summary = depthOfCov.sampleGeneSummary
		File sample_summary = depthOfCov.sampleSummary
		File sample_statistics = depthOfCov.sampleStatistics
		File sample_interval_summary = depthOfCov.sampleIntervalSummary
		File sample_interval_statistics = depthOfCov.sampleIntervalStatistics
		File sample_cumulative_coverage_proportions = depthOfCov.sampleCumulativeCoverageProportions
		File sample_cumulative_coverage_counts = depthOfCov.sampleCumulativeCoverageCounts

		# ContEst for Capture Array-Free output files
		File contamination = contestForCaptureArrayFree.contamination
		File contest_base_report = contestForCaptureArrayFree.baseReport
	}
}

task depthOfCov {
	File inputBam
	File inputBamIndex
	File geneList
	Int memoryGb
	Int disksGb
	Int minBaseQuality
	Int minMappingQuality
	String sampleName
	File refFasta
	File intervalList
	File refFastaDict
	File refFastaIndex

	command <<<
		java -Xmx15g -jar /usr/GenomeAnalysisTK.jar \
		-R ${refFasta} \
		-T DepthOfCoverage \
		-o ${sampleName} \
		-omitBaseOutput \
		-pt sample \
		-geneList ${geneList} \
		-I ${inputBam} \
		-L ${intervalList} \
		--minBaseQuality ${minBaseQuality} \
		--minMappingQuality ${minMappingQuality}
	>>>

	output {
		File sampleGeneSummary = "${sampleName}.sample_gene_summary"
		File sampleSummary = "${sampleName}.sample_summary"
		File sampleStatistics = "${sampleName}.sample_statistics"
		File sampleIntervalSummary = "${sampleName}.sample_interval_summary"
		File sampleIntervalStatistics = "${sampleName}.sample_interval_statistics"
		File sampleCumulativeCoverageProportions = "${sampleName}.sample_cumulative_coverage_proportions"
		File sampleCumulativeCoverageCounts = "${sampleName}.sample_cumulative_coverage_counts"
	}

	runtime {
		docker: "broadinstitute/gatk3:3.7-0"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${disksGb} HDD"
		
	}
}

task contestForCaptureArrayFree {
	File inputBam
	File inputBamIndex
	File normalBam
	File normalBamIndex
	File refFasta
	File refFastaIndex
	File refFastaDict
	File popVcf
	File dbSnpSixIntervalList
	File targetsIntervalList
	String sampleName
	Int memoryGb
	Int diskSpaceGb
	String inputBamName = basename(inputBam, ".bam")
	String normalBamName = basename(normalBam, ".bam")

	command <<<
	if [ ${inputBamName} == ${normalBamName} ]
	then
		java -Xmx4096m -jar /usr/local/bin/GenomeAnalysisTK.jar \
		-T ContEst \
	    -R ${refFasta} \
	    -I:eval,genotype ${inputBam} \
	    -l INFO \
	    -pf ${popVcf} \
	    -o ${sampleName}.contamination.txt \
	    -L ${dbSnpSixIntervalList} \
	    -L ${targetsIntervalList} \
	    -isr INTERSECTION \
	    --trim_fraction 0.03  \
		--beta_threshold 0.05 \
		-br ${sampleName}.contamination.txt.base_report.txt \
		-mbc 100  \
		--min_genotype_depth 30  \
		--min_genotype_ratio 0.8
	else
		java -Xmx4096m -jar /usr/local/bin/GenomeAnalysisTK.jar \
		-T ContEst \
	    -R ${refFasta} \
	    -I:eval ${inputBam} \
	    -I:genotype ${normalBam} \
	    -l INFO \
	    -pf ${popVcf} \
	    -o ${sampleName}.contamination.txt \
	    -L ${dbSnpSixIntervalList} \
	    -L ${targetsIntervalList} \
	    -isr INTERSECTION \
	    --trim_fraction 0.03  \
		--beta_threshold 0.05 \
		-br ${sampleName}.contamination.txt.base_report.txt \
		-mbc 100  \
		--min_genotype_depth 30  \
		--min_genotype_ratio 0.8
	fi
	>>>

	output {
		File contamination = "${sampleName}.contamination.txt"
		File baseReport = "${sampleName}.contamination.txt.base_report.txt"
	}

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
	
}