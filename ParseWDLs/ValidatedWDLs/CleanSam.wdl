workflow cleanSamWorkflow {
	call cleanSam
}

task cleanSam {
	File inputBam
	String sampleName
	Int diskSpaceGb
	Int memoryGb

	command <<<
	  java -jar /usr/gitc/picard.jar CleanSam \
	  I=${inputBam} \
	  O=${sampleName}.cleaned.bam
	>>>

	output {
		File cleanedBam = "${sampleName}.cleaned.bam"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}