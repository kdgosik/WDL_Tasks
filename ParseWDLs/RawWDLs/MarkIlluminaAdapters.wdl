workflow markIlluminaAdaptersTest {
	call markIlluminaAdapters
}

task markIlluminaAdapters {
	File inputBam
	String sampleName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 \
		-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m \
		-jar /usr/gitc/picard.jar MarkIlluminaAdapters \
		TMP_DIR=. \
		INPUT=${inputBam} \
		OUTPUT=${sampleName}.unmapped.bam \
		M=${sampleName}.adapter_metrics \
		PE=true \
		ADAPTERS=PAIRED_END
	>>>

	output {
		File unmappedBamFile = "${sampleName}.unmapped.bam"
		File adapterMetricsFile = "${sampleName}.adapter_metrics"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
		
	}
}