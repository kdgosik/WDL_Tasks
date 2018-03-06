workflow validateSamWorkflow {
	call validateSamFile
}

task validateSamFile {
	File inputBam
	File inputBamIndex
	String sampleName
	File refFasta
	File refFastaIndex
	File refFastaDict
	Boolean isBisulfateSequenced
	Int memoryGb
	Int diskSpaceGb


	command <<<
		java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar /usr/gitc/picard.jar ValidateSamFile \
		TMP_DIR=. \
		CREATE_MD5_FILE=false \
		I=${inputBam} \
		O=${sampleName}.validation_metrics \
		REFERENCE_SEQUENCE=${refFasta} \
		MODE=SUMMARY \
		IS_BISULFITE_SEQUENCED=${isBisulfateSequenced} \
		INDEX_VALIDATION_STRINGENCY=LESS_EXHAUSTIVE
	>>>

	output {
		File validationMetrics = "${sampleName}.validation_metrics"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"	
	}
}