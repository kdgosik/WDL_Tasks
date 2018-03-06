workflow facetsWorkflow {
	call facets
}

task facets {
	File tumorBam
	File tumorBamIndex
	File normalBam
	File normalBamIndex
	File refFasta
	String pairName
	Int diskSpaceGb
	Int memoryGb


	command <<<
		# move fasta for perl script
		mv ${refFasta} /usr/gitc/Homo_sapiens_assembly19.fasta

		# Remember starting directory b/c FireCloud doesn't let you put outputs of files w/ abs paths
		home_dir=$PWD

		cd /usr/gitc

		# Samtools mpileup wrapper
		perl /usr/gitc/scripts_snp_pileup/snp-pileup.pl ${normalBam} ${tumorBam} ${pairName}

		# mv for R script
		mv /usr/gitc/${pairName} /usr/gitc/pileup/${pairName}

		# Gotta add a little vanilla
		Rscript --vanilla /usr/gitc/facets.R

		cd /usr/gitc/pileup

		# move to home directory for output
		mv /usr/gitc/pileup/Facets_output.txt $home_dir/Facets_output.txt
		mv /usr/gitc/pileup/Facets_iterations.txt $home_dir/Facets_iterations.txt
		mv /usr/gitc/pileup/${pairName} $home_dir/${pairName}_pileup.txt

		cd $home_dir
	>>>

	output {
		File facetsOutput = "Facets_output.txt"
		File facetsIterations = "Facets_iterations.txt"
		File pileup = "${pairName}_pileup.txt"
	}

	runtime {
		docker: "jakeconway/facets:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}

}