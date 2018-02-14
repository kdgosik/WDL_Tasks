workflow sortSamWorkflow {
	call sortSam
}

task sortSam {
	File inputBam
	String sampleName
	Int diskSpaceGb
	Int memoryGb

	command <<<
	  java -jar /usr/gitc/picard.jar SortSam \
	  I=${inputBam} \
	  O=${sampleName}.sorted.bam \
	  SORT_ORDER=coordinate
	>>>

	output {
		File sortedBam = "${sampleName}.sorted.bam"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}

}