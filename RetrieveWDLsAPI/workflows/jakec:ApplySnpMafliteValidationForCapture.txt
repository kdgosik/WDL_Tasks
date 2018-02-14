workflow ApplySnpMafliteValidationForCaptureWorkflow {
	call ApplySnpMafliteValidationForCapture
}

task ApplySnpMafliteValidationForCapture {
	File inputMaf
	String matchmode
	String tumorSampleName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Xmx1g -jar /usr/gitc/ApplyMAFValidation.jar M=${inputMaf} \
		OUTPUT_MAF=${tumorSampleName}.maf.annotated MATCH_MODE=${matchmode} \
		V=/usr/gitc/validation
	>>>

	output {
		File validatedMafliteFile = "${tumorSampleName}.maf.annotated"
	}

	runtime {
		docker: "jakeconway/somatic_var_calling:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}