workflow printReadsWorkflow {
	call printReads
}

task printReads {
	File inputBam
	File inputBamIndex
	String sampleName
	File refFasta
	File refFastaIndex
	File refFastaDict
	File recalData
	String nameBase = basename(inputBam, ".bam")
	Int memoryGb
	Int diskSpaceGb

	
	command <<<
		mv ${inputBam} ${nameBase}.bam
		mv ${inputBamIndex} ${nameBase}.bai

		java -Djava.io.tmpdir=/directory -Dsamjdk.use_async_io=true -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8000m -jar /usr/GenomeAnalysisTK.jar -T PrintReads \
		--disable_auto_index_creation_and_locking_when_reading_rods \
		--generate_md5 -U \
		-R ${refFasta} \
		-I ${nameBase}.bam \
		--useOriginalQualities \
		-o ${sampleName}.bqsr.bam \
		--disable_indel_quals \
		-BQSR ${recalData} \
		--emit_original_quals
	>>>

	output {
		File printReadsBam = "${sampleName}.bqsr.bam"
	}

	runtime {
		docker: "broadinstitute/gatk3:3.7-0"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}

}