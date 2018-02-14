workflow gatk4CNVtangentNormalizationForCaptureWorkflow {
	call gatk4CNVtangentNormalizationForCapture
}

task gatk4CNVtangentNormalizationForCapture {
	File pCovFile
	File paddedTargetTsv
	File ponFile
	String sampleName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Xmx4g -jar /gatk/gatk.jar NormalizeSomaticReadCounts \
		--input ${pCovFile} \
		--targets ${paddedTargetTsv} \
		--panelOfNormals ${ponFile} \
		--tangentNormalized ${sampleName}.tn.tsv \
		--factorNormalizedOutput ${sampleName}.fnt.tsv \
		--betaHatsOutput ${sampleName}.beta_hats.tsv \
		--preTangentNormalized ${sampleName}.pre_tn.tsv
	>>>

	output {
		File tangentNormalized = "${sampleName}.tn.tsv"
		File factorNormalized = "${sampleName}.fnt.tsv"
		File betaHats = "${sampleName}.beta_hats.tsv"
		File preTangentNormalized = "${sampleName}.pre_tn.tsv"
	}

	runtime {
		docker: "broadinstitute/gatk:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}