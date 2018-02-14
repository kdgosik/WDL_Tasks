workflow BedToIntervalListWorkflow {
	call BedToIntervalList
}

task BedToIntervalList {
	File bedFile
	String bedFileName = basename(bedFile, ".bed")
	File refFastaDict
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -jar /usr/gitc/picard.jar BedToIntervalList \
		I=${bedFile} \
		O=${bedFileName}.interval_list \
		SD=${refFastaDict}
	>>>

	output {
		File intervalList = "${bedFileName}.interval_list"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}