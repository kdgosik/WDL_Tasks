workflow AddOrReplaceReadGroupsWorkflow {
	call AddOrReplaceReadGroups
}

task AddOrReplaceReadGroups {
	File inputBam
	String sampleName
	Int memoryGb
	Int diskSpaceGb
	
	command <<<
		java -jar /usr/gitc/picard.jar AddOrReplaceReadGroups \
		I=${inputBam} \
		O=${sampleName}.readgroupadded.bam \
		RGID=4 \
		RGLB=lib1 \
		RGPL=illumina \
		RGPU=unit1 \
		RGSM=20
	>>>

	output {
		File bamWithReadGroupAdded = "${sampleName}.readgroupadded.bam"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}