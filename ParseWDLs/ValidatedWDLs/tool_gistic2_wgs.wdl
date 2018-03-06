task tool_gistic2_wgs {
	File seg_file
	String markers_file
	File refgene_file
	File cnv_files

	Float amp_thresh
	Float del_thresh
	Float qv_thresh
	Float cap
	Float broad_length_cutoff
	Int remove_X
	Float conf_level
	Int join_segment_size
	Int arm_peel
	Int max_sample_segs
	Int do_gene_gistic
	String gene_collapse_method
	Int max_window_var

command {
	/src/call_gistic2 . 4 ${seg_file} ${markers_file} ${refgene_file} ${cnv_files} ${amp_thresh} ${del_thresh} ${qv_thresh } ${cap} ${broad_length_cutoff} ${remove_X} ${conf_level} ${join_segment_size} ${arm_peel} ${max_sample_segs} ${do_gene_gistic} ${gene_collapse_method} ${max_window_var} ./version.txt && /src/link_conf_wrapper.sh
}

output {
	File amp_genes="amp_genes.txt"
	File del_genes="del_genes.txt"
	File amp_qplot_png="amp_qplot.png"
	File del_qplot_png="del_qplot.png"
	File segmented_copy_number_png="raw_copy_number.png"
	File gistic_version="gisticVersion.txt"
	File arraylistfile="arraylistfile.txt"
	File broad_significance_results="broad_significance_results.txt"
	File gistic_inputs="gisticInputs.txt"
	File all_lesions="all_lesions.txt"
	File broad_values_by_arm="broad_values_by_arm.txt"
	File all_thresholded_by_genes="all_thresholded.by_genes.txt"
	File all_data_by_genes="all_data_by_genes.txt"
	File gistic_scores="scores.gistic"
}

	runtime {
		docker : "ggao/tool_gistic2_wgs:latest"
	}

	meta {
		author : "Galen Gao"
		email : "galengao@broadinstitute.org"
	}

}

workflow tool_gistic2_wgs_workflow {
	call tool_gistic2_wgs
}
