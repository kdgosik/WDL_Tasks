workflow collectAlignmentSummaryMetricsTest {

	call collectAlignmentSummaryMetrics
}

task collectAlignmentSummaryMetrics {
	File refFasta
	File inputBam
	String sampleName
	Boolean isBisulfateSequenced
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
		-Xmx1000m -jar /usr/gitc/picard.jar CollectAlignmentSummaryMetrics \
		TMP_DIR=. \
		INPUT=${inputBam} \
		OUTPUT=${sampleName}.alignment_summary_metrics \
		REFERENCE_SEQUENCE=${refFasta} \
		IS_BISULFITE_SEQUENCED=${isBisulfateSequenced}
	>>>

	output {
		File alignmentSummaryMetrics = "${sampleName}.alignment_summary_metrics"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}