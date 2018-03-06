workflow contestForCaptureArrayFreeWorkflow {
	File ref_fasta
	File ref_fasta_index
	File ref_fasta_dict
	File pop_vcf
	File dbsnp_six_interval_list

	call contestForCaptureArrayFree {
		input: refFasta=ref_fasta,
		refFastaIndex=ref_fasta_index,
		refFastaDict=ref_fasta_dict,
		popVcf=pop_vcf,
		dbSnpSixIntervalList=dbsnp_six_interval_list
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