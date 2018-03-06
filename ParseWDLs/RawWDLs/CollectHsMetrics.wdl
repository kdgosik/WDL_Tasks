workflow collectHsMetricsTest {
	File ref_fasta
	File ref_fasta_index
	File ref_fasta_dict
	File target_interval_list
	File bait_interval_list
	
	call collectHsMetrics {
		input: refFasta=ref_fasta,
		targetIntervalList=target_interval_list,
		baitIntervalList=bait_interval_list,
		refFastaDict=ref_fasta_dict,
		refFastaIndex=ref_fasta_index
	}
}

task collectHsMetrics {
	File inputBam
	File refFasta
	File refFastaIndex
	File refFastaDict
	File targetIntervalList
	File baitIntervalList
	String sampleName
	Int memory
	Int disks


	command <<<
		java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
		-Xmx1500m -jar /usr/gitc/picard.jar CollectHsMetrics \
		TMP_DIR=. \
		INPUT=${inputBam} \
		OUTPUT=${sampleName}.hybrid_selection_metrics \
		REFERENCE_SEQUENCE=${refFasta} \
		TARGET_INTERVALS=${targetIntervalList} \
		BAIT_INTERVALS=${baitIntervalList}
	>>>

	output {
		File HsMetrics = "${sampleName}.hybrid_selection_metrics"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memory} GB"
		cpu: "1"
		disks: "local-disk ${disks} HDD"
	}
}