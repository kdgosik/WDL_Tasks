workflow gatk4CNVproportionalCoverageForCaptureWorkflow {
	call gatk4CNVproportionalCoverageForCapture
}

task gatk4CNVproportionalCoverageForCapture {
	File refFasta
	File refFastaIndex
	File refFastaDict
	File paddedTargetTsv
	File inputBam
	String sampleName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -jar /gatk/gatk.jar CalculateTargetCoverage \
		--output ${sampleName}.pcov \
		--groupBy SAMPLE \
		--transform PCOV \
		--targets ${paddedTargetTsv} \
		--targetInformationColumns FULL \
		--input ${inputBam} \
		--reference ${refFasta} \
		--cohortName "<ALL>"
	>>>

	output {
		File pcov = "${sampleName}.pcov"
	}

	runtime {
		docker: "broadinstitute/gatk:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}

}