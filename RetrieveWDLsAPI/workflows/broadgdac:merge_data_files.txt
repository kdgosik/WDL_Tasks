task merge_data_files {
	#Inputs defined here
	String outfilename
	String filetype
	String row_column
	Int num_headers

	Array[File] sampledatafiles
	Array[String] sampleids

	command {
    python <<CODE
keys = '${sep="," sampleids}'.split(",")
values = '${sep="," sampledatafiles}'.split(",")
fout = open("datapathsfile.tsv", "w")
for i in range(len(keys)):
    fout.write(keys[i] + "\t" + values[i] + "\n")
fout.close()
CODE

/src/merge_txt_files "datapathsfile.tsv" ${outfilename}.${filetype}.txt ${row_column} ${num_headers}

	}

	output {
		File merged = "${outfilename}.${filetype}.txt"
	}

	runtime {
		docker : "broadgdac/merge_data_files:31"
	}

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	}

}

workflow merge_data_files_workflow {
	call merge_data_files
}