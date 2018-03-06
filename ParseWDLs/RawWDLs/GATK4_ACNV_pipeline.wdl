workflow GATK4_ACNV {
	String tumor_sample_name

	call gatk4AcnvGetBayesianHetsUsingTumorAndNormalForCapture

	call gatk4AcnvAllelicCnvForCapture {
		input: caseSampleName=tumor_sample_name,
		tumorHetFile=gatk4AcnvGetBayesianHetsUsingTumorAndNormalForCapture.tumorHets
	}

	call GATK4AcnvCallCnLohAndBalancedSegmentsForCapture {
		input: tumorSampleName=tumor_sample_name,
		hetFile=gatk4AcnvGetBayesianHetsUsingTumorAndNormalForCapture.tumorHets,
		segFile=gatk4AcnvAllelicCnvForCapture.tumorSigFinalSeg
	}

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
		docker: "broadinstitute/gatk:4.beta.6"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
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
		docker: "broadinstitute/gatk:4.beta.6"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}

task GATK4AcnvCallCnLohAndBalancedSegmentsForCapture {
	String tumorSampleName
	File hetFile
	File segFile
	File tnFile
	Float rhoThreshold
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Xmx4g -jar /gatk/gatk-protected.jar CallCNLoHAndSplits \
		--segments ${segFile} \
		--tumorHets ${hetFile} \
		--tangentNormalized ${tnFile} \
		--rhoThreshold ${rhoThreshold} \
		--outputDir .

		ls
	>>>

	output {
		File finalAcsSeg = "${tumorSampleName}-sim-final.acs.seg"
		File finalCnbSeg = "${tumorSampleName}-sim-final.cnb_called.seg"
		File finalCnvSeg = "${tumorSampleName}-sim-final.cnv.seg"
		File finalTitanHet = "${tumorSampleName}-sim-final.titan.het.tsv"
		File finalTitanTn = "${tumorSampleName}-sim-final.titan.tn.tsv"

	}

	runtime {
		docker: "jakeconway/gatk-protected:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
	
}