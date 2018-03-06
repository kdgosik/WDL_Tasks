task report_gistic2 {
	#Inputs defined here
	String? type

	File? amp
	File? del
	File amp_qplot_png
	File del_qplot_png
	File segmented_copy_number_png
	File gistic_version
	File arraylistfile_txt
	File broad_sig_res
	File gistic_inputs

	command <<<
		ln -s ${amp_qplot_png} amp_qplot.png
		ln -s ${del_qplot_png} del_qplot.png
		ln -s ${segmented_copy_number_png} segmented_copy_number.png
        /R/RunR.sh -f main /src/Gistic2Report.R main -g/src/geneList.txt ${"-a" + amp} ${"-d" + del} -Aamp_qplot.png -Ddel_qplot.png -Ssegmented_copy_number.png -V${gistic_version} -T${arraylistfile_txt} -b${broad_sig_res} -I${gistic_inputs}
        TS=`date` ; NEW_DIV="<div class=\"sectionbody\">$TS</div>"  ; 
        awk -v NEW_DIV="$NEW_DIV" '{if($0~/class=.title/) { print NEW_DIV "\n" $0 } else { print $0}}' nozzle.html > nozzle_date.html ;
        mv -vf nozzle_date.html nozzle.html
	    >>>

	output {
		#Outputs defined here
		File nozzle_html="nozzle.html"
		File nozzle_Rdata="nozzle.RData"
		File amp_qplot_png_ignore="amp_qplot.png"
		File del_qplot_png_ignore="del_qplot.png"
		File segmented_copy_number_png_ignore="segmented_copy_number.png"

	    }

	runtime {
		docker : "broadgdac/report_gistic2:67"
	    }

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	    }

}


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
	Int memoryGB
	Int preemptible

	parameter_meta {

		seg_file: "A six-column tab-delimited file containing the segmented data for all tumor/normal pairs in the pair set."
		markers_file: "A three-column tab-delimited file identifying the names and positions all markers."
		refgene_file: "Contains information about the location of genes and cytobands on a given build of the genome.  These files are created in MATLAB."
		cnv_files: ""
		amp_thresh: "Threshold for copy number amplifications.  Regions with a log2 ratio above this value are considered amplified. (Recommended: 0.1)"
		del_thresh: "Threshold for copy number deletions.  REgions with a log2 ratio below the negative of this value are considered deletions.(Recommended: 0.1)"
		qv_thresh: "Significance threshold for Q-values.  Regions with Q-values below this number are considered signficant. (Recommended: 0.1)"
		cap: "Minimum and maximum cap values on analyzed data. Regions with a log2 ratio greater than the cap are set to the cap value; regions with a log2 ratio less than -cap value are set to -cap. (DEFAULT=1.5)"
		broad_length_cutoff: "Threshold used to distinguish broad from focal events, given in units of fraction of chromsome arm.  (Recommended: 0.7)"
		remove_X: "0/1 flag indicating whether to remove data from the X chromosome before analysis.  (Recommended: 0)"
		conf_level: "Confidence level used to calculate region containing the driver.  (Recommended: 0.99)"
		join_segment_size: "Smallest number of markers to allow in segements from the segemented data.  Segements that contain a number of markers less than or equal to this number are joined to the adjacent segement, closest in copy number. (Recommended: 4)"
		arm_peel: "0/1 flag indicating whether to perform arm-level peel off, wich helps spearte peaks and clean up noise.  (Recommended: 1)"
		max_sample_segs: "Maximum number of segements allowed for a smaple in the input data.    Samples with more segments than this are excluded from the analysis. (Recommended: 2000)"
		do_gene_gistic: "0/1 flag indicating tht the gene GISTIC aglrithm should be used to calculate signficance of deletions at the gene level instead of a marker level. (Recommended: 1)"
		gene_collapse_method: "Method for reducing marker-level copy number data to the gene-level copy number data in the gene tables. Markers contained in the gene are used when available, otherwise the flanking marker or markers are used. Allowed values are mean, median, min, max or extreme. The extreme method chooses whichever of min or max is furthest from diploid. (Recommended: extreme)"
		memoryGB: "Integer value specifying the minimum memory requirements (in GB) for the virtual machine running the GISTIC2 task."
		preemptible: "Integer value specifying the maximum number of times Cromwell should request a preemptible machine for this task before defaulting back to a non-preemptible one."
	}

	command {

		#write runtime info to file
		echo -n "Number of Processors: " > runtime_info.txt; 
		nproc >> runtime_info.txt
		cat /proc/meminfo | grep MemTotal >> runtime_info.txt
		df >> runtime_info.txt

	    /src/call_gistic2 . 4 ${seg_file} ${markers_file} ${refgene_file} ${cnv_files} ${amp_thresh} ${del_thresh} ${qv_thresh } ${cap} ${broad_length_cutoff} ${remove_X} ${conf_level} ${join_segment_size} ${arm_peel} ${max_sample_segs} ${do_gene_gistic} ${gene_collapse_method} ./version.txt && /src/link_conf_wrapper.sh
	}

	output {
		File runtime_info="runtime_info.txt"
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
		memory: "${memoryGB} GB"
		preemptible: "${preemptible}"
	}

	meta {
	     author : "Tim DeFreitas"
	     email : "timdef@broadinstitute.org"
	}

}

workflow copy_number_gistic2 {
	call tool_gistic2
	call report_gistic2 {
        input:
        amp=tool_gistic2.amp,
        del=tool_gistic2.del,
        amp_qplot_png=tool_gistic2.amp_qplot_png,
        del_qplot_png=tool_gistic2.del_qplot_png,
        segmented_copy_number_png=tool_gistic2.segmented_copy_number_png,
        gistic_version=tool_gistic2.gistic_version,
        arraylistfile_txt=tool_gistic2.arraylistfile,
        broad_sig_res=tool_gistic2.broad_sig_res,
        gistic_inputs=tool_gistic2.gistic_inputs
	}

}