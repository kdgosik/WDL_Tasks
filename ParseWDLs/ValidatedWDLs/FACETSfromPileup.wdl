workflow facetsFromPileupWorkflow {
	call facetsFromPileup
}

task facetsFromPileup {
	File pileupFile
	String pairName
	Int diskSpaceGb
	Int memoryGb


	command <<<
		# Remember starting directory b/c FireCloud doesn't let you put outputs of files w/ abs paths
		home_dir=$PWD

		# mv for R script
		mv ${pileupFile} /usr/gitc/pileup/${pairName}

		# Gotta add a little vanilla
		Rscript --vanilla /usr/gitc/facets.R

		cd /usr/gitc/pileup

		# move to home directory for output
		mv /usr/gitc/pileup/Facets_output.txt $home_dir/Facets_output.txt
		mv /usr/gitc/pileup/Facets_iterations.txt $home_dir/Facets_iterations.txt

		cd $home_dir
	>>>

	output {
		File facetsOutput = "Facets_output.txt"
		File facetsIterations = "Facets_iterations.txt"
	}

	runtime {
		docker: "jakeconway/facets:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}

}