task tool_gistic2 {
	#Inputs defined here
	File seg_file 
	File markers_file 
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

	command {
		/src/call_gistic2 . 4 ${seg_file} ${markers_file} ${refgene_file} ${cnv_files} ${amp_thresh} ${del_thresh} ${qv_thresh } ${cap} ${broad_length_cutoff} ${remove_X} ${conf_level} ${join_segment_size} ${arm_peel} ${max_sample_segs} ${do_gene_gistic} ${gene_collapse_method} ./version.txt && /src/link_conf_wrapper.sh
	}

	output {
		File amp="amp_genes.txt"
		File del="del_genes.txt"
		File amp_qplot_png="amp_qplot.png"
		File del_qplot_png="del_qplot.png"
		File segmented_copy_number_png="raw_copy_number.png"
		File gistic_version="gisticVersion.txt"
		File arraylistfile="arraylistfile.txt"
		File broad_sig_res="broad_significance_results.txt"
		File gistic_inputs="gisticInputs.txt"
	}

	runtime {
		docker : "broadgdac/tool_gistic2:141"
	}

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	}

}

workflow tool_gistic2_workflow {
	call tool_gistic2
}
