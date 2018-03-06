workflow collectQualityYieldMetricsWorkflow {
	call collectQualityYieldMetrics
}

task collectQualityYieldMetrics {
	File inputBam
	String sampleName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx128m -jar /usr/gitc/picard.jar CollectQualityYieldMetrics \
		TMP_DIR=. \
		I=${inputBam} \
		O=${sampleName}_quality_yield_metrics.txt
	>>>

	output {
		File metricsFile = "${sampleName}_quality_yield_metrics.txt"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}