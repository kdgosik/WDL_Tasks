task preprocess_clinical {
  File stage2_file
  Int node_column
  String output_prefix

	command {
    Rscript /src/clinicalDataPicker.R --libdir /src/ --stage2File ${stage2_file} --nodeColumn ${node_column} --outPrefix ${output_prefix}
  }

  output {
    File tier1_clinical_data="${output_prefix}.clin.merged.picked.txt"
	}

	runtime {
		docker : "broadgdac/preprocess_clinical:113"
	}

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	}

}

workflow preprocess_clinical_workflow {
	call preprocess_clinical
}
