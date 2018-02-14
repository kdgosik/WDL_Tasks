task tool_gistic2_postprocess {
  File cnfile
  String outputprefix

	command {
    Rscript --vanilla /src/CN_postprocess.R -c ${cnfile} -o ${outputprefix}
	}

	output {
    File gistic_output_all_lesions_realcn="${outputprefix}.cn_focal_gistic.txt"
    File gistic_output_all_lesions_threshold="${outputprefix}.cn_focal_gistic_threshold.txt"
	}

	runtime {
		docker : "broadgdac/tool_gistic2_postprocess:13"
	}

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	}

}

workflow tool_gistic2_postprocess_workflow {
	call tool_gistic2_postprocess
}
