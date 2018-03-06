task aggregate_clusters {
	#Inputs defined here
	File? miRseqMaturecnmf
	File? miRseqMaturehc
	File? mrnacnmf
	File? mrnahc
	File? mirnacnmf
	File? mirnahc
	File? copynumber
	File? methylationhc
	File? methylationcnmf
	File? PredictionU133aTCGA
	File? PredictionU133aPhilips
	File? PredictionUNCTCGA
	File? mRNAseqcnmf
	File? mRNAseqhc
	File? miRseqcnmf
	File? miRseqhc
	File? CustomEvents

	String mergedfile



	command {
		/R/RunR.sh -f main /src/aggregate_clusters.R ${"-w" + miRseqMaturecnmf} ${"-y" + miRseqMaturehc} ${"-c" + mrnacnmf} ${"-d" + mrnahc} ${"-e" + mirnacnmf} ${"-f" + mirnahc} ${"-g" + copynumber} ${"-h" + methylationhc} ${"-a" + methylationcnmf} ${"-i" + mergedfile} ${"-b" + PredictionU133aTCGA} ${"-j" + PredictionU133aPhilips} ${"-k" + PredictionUNCTCGA} ${"-l" + mRNAseqcnmf} ${"-m" + mRNAseqhc} ${"-n" + miRseqcnmf} ${"-p" + miRseqhc} ${"-x" + CustomEvents}
	}

	output {
		#Outputs defined here
		File aggregate_molecular_clusters="${mergedfile}.mergedcluster.txt"
		File aggregate_molecular_clusters_transposed="${mergedfile}.transferedmergedcluster.txt"
	}

	runtime {
		docker : "broadgdac/aggregate_clusters:3"
	}

	meta {
		author : "timdef"
		email : "timdef@broadinstitute.org"
	}

}

workflow aggregate_clusters_workflow {
	call aggregate_clusters
}
