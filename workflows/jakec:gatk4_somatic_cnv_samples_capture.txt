workflow GATK4SomaticCnvToolchainCaptureForSamples {
	String sample_name
	File ref_fasta
	File ref_fasta_index
	File ref_fasta_dict
	File input_bam
	File input_bam_idx

	call gatk4CNVpadTargetsForCapture {
		input: sampleName=sample_name
	}

	call gatk4CNVproportionalCoverageForCapture  {
		input: sampleName=sample_name,
		refFasta=ref_fasta,
		refFastaIndex=ref_fasta_index,
		refFastaDict=ref_fasta_dict,
		inputBam=input_bam,
		paddedTargetTsv=gatk4CNVpadTargetsForCapture.paddedTsv
	}

	call gatk4CNVtangentNormalizationForCapture {
		input:sampleName=sample_name,
		pCovFile=gatk4CNVproportionalCoverageForCapture.pcov,
		paddedTargetTsv=gatk4CNVpadTargetsForCapture.paddedTsv
	}

	call performSegmentation {
		input:sampleName=sample_name,
		tangentNormalizedFile=gatk4CNVtangentNormalizationForCapture.tangentNormalized
	}

	call callSegments {
		input:tangentNormalizedFile=gatk4CNVtangentNormalizationForCapture.tangentNormalized,
		segmentFile=performSegmentation.segFile,
		sampleName=sample_name
	}

	call plotSegmentedCopyRatio {
		input:tangentNormalizedFile=gatk4CNVtangentNormalizationForCapture.tangentNormalized,
		segmentFile=performSegmentation.segFile,
		refFastaDict=ref_fasta_dict,
		sampleName=sample_name,
		preTangentNormalizedFile=gatk4CNVtangentNormalizationForCapture.preTangentNormalized
	}

	call oncotateGatkCnvFileForCapture {
		input:calledSegmentsFile=callSegments.calledSegFile,
		sampleName=sample_name
	}

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

task oncotateGatkCnvFileForCapture {
	String txMode
	String genomeBuild
	String sampleName
	String inputFormat
	String outputFormat
	File dataSourceTarGzFile
	String dataSourceDir = basename(dataSourceTarGzFile, ".tar.gz")
	File defaultConfig
	File calledSegmentsFile
	Int memoryGb
	Int diskSpaceGb

	command <<<
		
		cp -r /root/oncotator_venv/ $PWD
		ls -a

		tar -xzvf ${dataSourceTarGzFile}

		./oncotator_venv/bin/oncotator -v -i ${inputFormat} \
		-o ${outputFormat} \
		--db-dir ${dataSourceDir}/ \
		--no-multicore \
		--default_config ${defaultConfig} ${calledSegmentsFile} ${sampleName}.called.seg.annotated ${genomeBuild} \
		--tx-mode ${txMode}

	>>>

	output {
		File annotatedSegFile = "${sampleName}.called.seg.annotated"
	}

	runtime {
		docker: "broadinstitute/oncotator:1.9.3.0"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}