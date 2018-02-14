task plotSegmentedCopyRatio {
	File tangentNormalizedFile
	File segmentFile
	File refFastaDict
	String sampleName
	File preTangentNormalizedFile
	Boolean log2Input
	Int memoryGb
	Int diskSpaceGb

	command <<<
		mkdir ${sampleName}

		java -jar /gatk/gatk.jar PlotSegmentedCopyRatio \
		--tangentNormalized ${tangentNormalizedFile} \
		--segments ${segmentFile} \
		--output ${sampleName} \
		--outputPrefix ${sampleName} \
		--preTangentNormalized ${preTangentNormalizedFile} \
		--sequenceDictionaryFile ${refFastaDict} \
		--log2Input ${log2Input}
	>>>

	runtime {
		docker: "broadinstitute/gatk:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}

task callSegments {
	File tangentNormalizedFile
	File segmentFile
	String sampleName
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -jar /gatk/gatk.jar CallSegments \
		--tangentNormalized ${tangentNormalizedFile} \
		--segments ${segmentFile} \
		--output ${sampleName}.called
	>>>


	output {
		File calledSegFile = "${sampleName}.called"
	}

	runtime {
		docker: "broadinstitute/gatk:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}

task performSegmentation {
	File tangentNormalizedFile
	String sampleName
	Float alpha
	Float eta
	Float trim
	Float undoPrune
	Int nperm
	Int minWidth
	Int kmax
	Int nmin
	Int undoSD
	Boolean log2Input
	Int memoryGb
	Int diskSpaceGb
	
	command <<<
		java -Xmx4g -jar /gatk/gatk.jar PerformSegmentation \
		--tangentNormalized ${tangentNormalizedFile} \
		--output ${sampleName}.seg \
		--alpha ${alpha} \
		--nperm ${nperm} \
		--pmethod HYBRID \
		--minWidth ${minWidth} \
		--kmax ${kmax} \
		--nmin ${nmin} \
		--eta ${eta} \
		--trim 	${trim} \
		--undoSplits NONE \
		--undoPrune ${undoPrune} \
		--undoSD ${undoSD} \
		--log2Input ${log2Input}
	>>>

	output {
		File segFile = "${sampleName}.seg"
	}

	runtime {
		docker: "broadinstitute/gatk:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}

workflow gatk4CNVsegmentThenCallThenPlotForCaptureWorkflow {
	String sample_name
	File tangent_normalized_file

	call performSegmentation {
		input:sampleName=sample_name,
		tangentNormalizedFile=tangent_normalized_file
	}
	call callSegments {
		input: segmentFile=performSegmentation.segFile,
		sampleName=sample_name,
		tangentNormalizedFile=tangent_normalized_file
	}
	call plotSegmentedCopyRatio {
		input: sampleName=sample_name,
		tangentNormalizedFile=tangent_normalized_file,
		segmentFile=performSegmentation.segFile
	}
}