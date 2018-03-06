workflow ConvertBedToTargetFileWorkflow {
	call ConvertBedToTargetFile
}

task ConvertBedToTargetFile {
	File targetBed
	String targetBedName = basename(targetBed, ".bed")
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Xmx4g -jar /gatk/gatk.jar ConvertBedToTargetFile \
		--input ${targetBed} \
		--output ${targetBedName}.tsv
	>>>

	output {
		File targetsTsv = "${targetBedName}.tsv"
	}

	runtime {
		docker: "broadinstitute/gatk:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}