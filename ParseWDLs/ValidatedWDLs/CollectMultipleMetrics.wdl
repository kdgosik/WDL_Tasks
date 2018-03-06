workflow collectMultipleMetricsTest {
	call collectMultipleMetrics
}

task collectMultipleMetrics {
	File inputBam
	File refFasta
	String sampleName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
		-Xmx4000m -jar /usr/gitc/picard.jar CollectMultipleMetrics \
		TMP_DIR=. \
		METRIC_ACCUMULATION_LEVEL=null \
		METRIC_ACCUMULATION_LEVEL=ALL_READS \
		PROGRAM=null \
		PROGRAM=MeanQualityByCycle \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=CollectInsertSizeMetrics \
		PROGRAM=CollectAlignmentSummaryMetrics \
		INPUT=${inputBam} \
		REFERENCE_SEQUENCE=${refFasta} \
		ASSUME_SORTED=true \
		OUTPUT=${sampleName}.multiple_metrics.bam
	>>>

	output {
		File cycleMetrics = "${sampleName}.multiple_metrics.bam.quality_by_cycle_metrics"
		File distributionMetrics = "${sampleName}.multiple_metrics.bam.quality_distribution_metrics"
		File insertSizeMetrics = "${sampleName}.multiple_metrics.bam.insert_size_metrics"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"	
	}
}