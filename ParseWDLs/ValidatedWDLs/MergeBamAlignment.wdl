workflow mergeBamAlignmentTest {
	call mergeBamAlignment
}

task mergeBamAlignment {
	File refFasta
	File refFastaDict
	File refFastaIndex
	File unmappedBam
	File alignedBam
	File unpairedAlignedBam
	String sampleName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:+UseStringCache -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx5000m -jar /usr/gitc/picard.jar MergeBamAlignment \
		TMP_DIR=. \
		VALIDATION_STRINGENCY=SILENT \
		CREATE_INDEX=true \
		ALIGNED_BAM=${alignedBam} \
		ALIGNED_BAM=${unpairedAlignedBam} \
		EXPECTED_ORIENTATIONS=FR \
		ATTRIBUTES_TO_RETAIN=X0 \
		UNMAPPED_BAM=${unmappedBam} \
		OUTPUT=${sampleName}.merged.aligned.bam \
		REFERENCE_SEQUENCE=${refFasta} \
		PAIRED_RUN=true \
		IS_BISULFITE_SEQUENCE=false \
		ALIGNED_READS_ONLY=false \
		CLIP_ADAPTERS=false \
		MAX_RECORDS_IN_RAM=2000000 \
		ADD_MATE_CIGAR=true
	>>>

	output {
		File mergedAlignmentBam = "${sampleName}.merged.aligned.bam"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
		
}