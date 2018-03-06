task maf_aggregator_task
	{

	String memGB
	String diskGB
	String pSetID
    String cpus
	Array[File]+ mafs

	command <<<

	#increase verbosity
	set -x

       python <<CODE
NULL_SENTINEL = "GDAC_FC_NULL"

mafpaths = '${sep="," mafs}'.split(",")
with open("datapathsfile.tsv", "w") as fout:
    for mafpath in mafpaths:
        if not mafpath.endswith(NULL_SENTINEL):
            fout.write("NA\t%s\n"%mafpath)
CODE


	#run catters
	python /usr/local/bin/tsvConcatListFile.py datapathsfile.tsv ${pSetID}.aggregated.maf

	>>>

	output {
		File aggregated_maf="${pSetID}.aggregated.maf"
		}

	runtime {
		docker: "broadinstitute/broadmutationcalling_filtering_beta@sha256:d2df6d9d705e90d3ee926b72a3edec5034dd5bdd64c0cf1cabd9fc8462702c79"
		memory: "${memGB} GB"
        cpu: "${cpus}"
		disks: "local-disk ${diskGB} HDD"
		}

	}

workflow maf_aggregator_workflow {

	call maf_aggregator_task

}