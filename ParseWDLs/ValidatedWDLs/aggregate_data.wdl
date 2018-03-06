task aggregate_data {
	#Inputs defined here
	String outfilename
	String filetype
	String row_column
	Int num_headers

	Array[File] sampledatafiles
	Array[String] sampleids

	command {
    python <<CODE
NULL_SENTINEL = "GDAC_FC_NULL"

keys = '${sep="," sampleids}'.split(",")
values = '${sep="," sampledatafiles}'.split(",")
fout = open("datapathsfile.tsv", "w")
for i in range(len(keys)):
    if not values[i].endswith(NULL_SENTINEL):
        fout.write(keys[i] + "\t" + values[i] + "\n")
fout.close()
CODE

/src/merge_txt_files "datapathsfile.tsv" ${outfilename}.${filetype}.txt ${row_column} ${num_headers}

	}

	output {
		File merged = "${outfilename}.${filetype}.txt"
	}

	runtime {
		docker : "broadgdac/aggregate_data:32"
	}

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	}

}

workflow aggregate_data_workflow {
	call aggregate_data
}
