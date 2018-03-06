workflow gatk4CNVpadTargetsForCaptureWorkflow {
	call gatk4CNVpadTargetsForCapture
}

task gatk4CNVpadTargetsForCapture {
	File targetsTsv
	Int padding
	String sampleName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Xmx1g -jar /gatk/gatk.jar PadTargets \
		--targets ${targetsTsv} \
		--output ${sampleName}.capture.padded.tsv \
		--padding ${padding}
	>>>

	output {
		File paddedTsv = "${sampleName}.capture.padded.tsv"
	}

	runtime {
		docker: "broadinstitute/gatk:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
	
}