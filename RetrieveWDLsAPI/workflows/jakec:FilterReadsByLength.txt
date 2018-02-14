workflow trimLengthZeroReads {
	call trimReads
}

task trimReads {
	File inputBam
	String bamName = basename(inputBam, ".bam")
	Int readLength
	Int memoryGb

	Float bamSize = size(inputBam, "GB")
	Float diskSpace = (bamSize * 3)
	String float_to_int_diskSpace = sub(diskSpace, "\\..*", "")
	Int diskSpaceGb = float_to_int_diskSpace


	command <<<
	samtools
	samtools view -h ${inputBam} | awk 'length($10) > ${readLength} || $1 ~ /^@/' | samtools view -bS - > ${bamName}.filtered.bam
	>>>

	output {
		File filteredBam = "${bamName}.filtered.bam"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}