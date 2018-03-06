workflow markDuplicatesTest {
	call markDuplicates
}

task markDuplicates {
	File inputBam
	String sampleName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar /usr/gitc/picard.jar MarkDuplicates \
		TMP_DIR=. \
		CREATE_INDEX=true \
		CREATE_MD5_FILE=true \
		INPUT=${inputBam} \
		OUTPUT=${sampleName}.aligned.duplicates_marked.bam \
		METRICS_FILE=${sampleName}.duplicate_metrics \
		VALIDATION_STRINGENCY=SILENT \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
		READ_ONE_BARCODE_TAG=RX
	>>>

	output {
		File markDuplicatesMetrics = "${sampleName}.duplicate_metrics"
		File markedDuplicatesBam = "${sampleName}.aligned.duplicates_marked.bam"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}

}