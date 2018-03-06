workflow depthOfCoverageTest {
	File ref_fasta
	File gene_list
	File interval_list
	File ref_fasta_index
	File ref_fasta_dict

	call depthOfCov {
		input: refFasta=ref_fasta,
			geneList=gene_list,
			intervalList=interval_list,
			refFastaIndex=ref_fasta_index,
			refFastaDict=ref_fasta_dict
	}
}

task depthOfCov {
	File inputBam
	File inputBamIndex
	Int minBaseQuality
	Int minMappingQuality
	String sampleName
	File intervalList
	File geneList
	File refFasta
	File refFastaDict
	File refFastaIndex
	Int diskSpaceGb
	Int memoryGb

	String bamName = basename(inputBam)
	String baiName = basename(inputBamIndex)

	command <<<
		mv ${inputBam} $PWD
		mv ${inputBamIndex} $PWD

		bamFile=$PWD/${bamName}
		baiFile=$PWD/${baiName}



		java -Xmx15g -jar /usr/GenomeAnalysisTK.jar \
		-R ${refFasta} \
		-T DepthOfCoverage \
		-o ${sampleName} \
		-omitBaseOutput \
		-pt sample \
		-geneList ${geneList} \
		-I $bamFile \
		-L ${intervalList} \
		--minBaseQuality ${minBaseQuality} \
		--minMappingQuality ${minMappingQuality}
	>>>

	output {
		File sampleGeneSummary = "${sampleName}.sample_gene_summary"
		File sampleSummary = "${sampleName}.sample_summary"
		File sampleStatistics = "${sampleName}.sample_statistics"
		File sampleIntervalSummary = "${sampleName}.sample_interval_summary"
		File sampleIntervalStatistics = "${sampleName}.sample_interval_statistics"
		File sampleCumulativeCoverageProportions = "${sampleName}.sample_cumulative_coverage_proportions"
		File sampleCumulativeCoverageCounts = "${sampleName}.sample_cumulative_coverage_counts"
	}

	runtime {
		docker: "broadinstitute/gatk3:3.7-0"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
		
	}
}
