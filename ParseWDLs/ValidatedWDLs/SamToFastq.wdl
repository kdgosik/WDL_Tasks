workflow samToFastqTest {
	call samToFastq
}

task samToFastq {
	File inputBam
	String sampleName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Dsamjdk.buffer_size=131072 -Dsamjdk.use_async_io=true \
		-Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx256m \
		-jar /usr/gitc/picard.jar SamToFastq \
		TMP_DIR=. \
		INPUT=${inputBam} \
		FASTQ=${sampleName}.1.fastq.gz \
		INTERLEAVE=false \
		SECOND_END_FASTQ=${sampleName}.2.fastq.gz \
		INCLUDE_NON_PF_READS=true \
		CLIPPING_ATTRIBUTE=XT \
		CLIPPING_ACTION=2 \
		UNPAIRED_FASTQ=${sampleName}.unpaired.fastq.gz
	>>>

	output {
		File firstEndFastq = "${sampleName}.1.fastq.gz"
		File secondEndFastq = "${sampleName}.2.fastq.gz"
		File unpairedFastq = "${sampleName}.unpaired.fastq.gz"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
		
	}
}