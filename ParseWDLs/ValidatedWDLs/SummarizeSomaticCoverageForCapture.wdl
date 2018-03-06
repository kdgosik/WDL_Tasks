workflow SummarizeSomaticCoverageForCaptureWorkflow {
	call SummarizeSomaticCoverageForCapture
}

task SummarizeSomaticCoverageForCapture {
	File coverageWigFile
	String pairName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		perl /usr/gitc/summarizeWigFile.pl ${coverageWigFile} ${pairName}.somatic_coverage_summary.txt
	>>>

	output {
		File coverageSummaryFile = "${pairName}.somatic_coverage_summary.txt"
	}

	runtime {
		docker: "jakeconway/somatic_var_calling:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}