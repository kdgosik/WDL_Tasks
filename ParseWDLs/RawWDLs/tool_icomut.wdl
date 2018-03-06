task tool_icomut {
	File sample_feature_table
	String cohort
	Float primary_cutoff

	File? clinical_features

	command {
		Rscript /src/generate_comut.R -s /src -f ${sample_feature_table} -t ${cohort} ${"-c" + clinical_features} -x ${primary_cutoff}
	}

	output {
		#Outputs defined here
		File iCoMut_table="${cohort}.coMut_table.txt"
	}

	runtime {
		docker : "broadgdac/tool_icomut:1"
	}

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	}

}

workflow tool_icomut_workflow {
	call tool_icomut
}
