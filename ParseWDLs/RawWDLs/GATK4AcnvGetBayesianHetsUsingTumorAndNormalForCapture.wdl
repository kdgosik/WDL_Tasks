workflow gatk4AcnvGetBayesianHetsUsingTumorAndNormalForCaptureWorkflow {
	call gatk4AcnvGetBayesianHetsUsingTumorAndNormalForCapture
}

task gatk4AcnvGetBayesianHetsUsingTumorAndNormalForCapture {
	File normalBam
	File normalBamIndex
	File tumorBam
	File tumorBamIndex
	File refFasta
	File refFastaIndex
	File refFastaDict
	File snpIntervalList
	String pairName
	Int readDepthThreshold
	Int minMq
	Int minBq
	Float hetStringency
	Float minAbnormalFraction
	Float maxAbnormalFraction
	Int maxCopyNumber
	Int quadratureOrder
	Float errorAdjustment
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -jar /gatk/gatk.jar GetBayesianHetCoverage \
		--normal ${normalBam} \
		--tumor ${tumorBam} \
		--reference ${refFasta} \
		--snpIntervals ${snpIntervalList} \
		--normalHets ${pairName}.normal.hets.tsv \
		--tumorHets ${pairName}.tumor.hets.tsv \
		--readDepthThreshold ${readDepthThreshold} \
		--minimumMappingQuality ${minMq} \
		--minimumBaseQuality ${minBq} \
		--hetCallingStringency ${hetStringency} \
		--minimumAbnormalFraction ${minAbnormalFraction} \
		--maximumAbnormalFraction ${maxAbnormalFraction} \
		--maximumCopyNumber ${maxCopyNumber} \
		--quadratureOrder ${quadratureOrder} \
		--errorAdjustmentFactor ${errorAdjustment}

	>>>

	output {
		File normalHets = "${pairName}.normal.hets.tsv"
		File tumorHets = "${pairName}.tumor.hets.tsv"
	}

	runtime {
		docker: "broadinstitute/gatk:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}