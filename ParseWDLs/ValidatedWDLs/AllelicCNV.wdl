workflow gatk4AcnvAllelicCnvForCaptureWorkflow {
	call gatk4AcnvAllelicCnvForCapture
}

task gatk4AcnvAllelicCnvForCapture {
	File tumorSegFile
	File tumorHetFile
	File tumorTangentNormalizedFile
	String caseSampleName
	Int smallSegThreshold
	Int numOfSamplesCopyRatio
	Int numBurninCopyRatio
	Int numSamplesAlleleFraction
	Int numBurninAlleleFraction
	Float intervalThresholdCopyRatio
	Float intervalThresholdAlleleFraction
	Int maxNumIterSimSeg
	Int maxNumIterSnpSeg
	Boolean useAllCopyRatioSegs
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -jar /gatk/gatk.jar AllelicCNV \
		--tumorHets ${tumorHetFile} \
		--tangentNormalized ${tumorTangentNormalizedFile} \
		--segments ${tumorSegFile} \
		--outputPrefix ${caseSampleName} \
		--smallSegmentThreshold ${smallSegThreshold} \
		--numSamplesCopyRatio ${numOfSamplesCopyRatio} \
		--numBurnInCopyRatio ${numBurninCopyRatio} \
		--numSamplesAlleleFraction ${numSamplesAlleleFraction} \
		--numBurnInAlleleFraction ${numBurninAlleleFraction} \
		--intervalThresholdCopyRatio ${intervalThresholdCopyRatio} \
		--intervalThresholdAlleleFraction ${intervalThresholdAlleleFraction} \
		--maxNumIterationsSimSeg ${maxNumIterSimSeg} \
		--maxNumIterationsSNPSeg ${maxNumIterSnpSeg} \
		--useAllCopyRatioSegments ${useAllCopyRatioSegs}
	>>>

	output {
		File tumorSigFinalSeg = "${caseSampleName}-sim-final.seg"
	}

	runtime {
		docker: "broadinstitute/gatk:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}