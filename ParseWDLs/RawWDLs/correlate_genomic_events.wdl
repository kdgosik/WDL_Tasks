task correlate_genomic_event_clinical {
    
    File  cluster_label_file
    File  clinical_info_file 
    String? features_used # optional
    String filter_cluster_smaller_than
    String figure_height
    String figure_width 
    String survival_legend_position 
    String pvalue_cutoff 
    String qvalue_cutoff
    String make_figure
    String output_name 
    String options 
    String? report_title # optional 
    String? job_description # optional

	command {
		/R/RunR.sh -f main /src/Correlate_Genomic_Events.R /src -iD${cluster_label_file} -iC${clinical_info_file} -fF${features_used} -fC${filter_cluster_smaller_than} -fH${figure_height} -fW${figure_width} -fP${survival_legend_position} -cP${pvalue_cutoff} -cQ${qvalue_cutoff} -MF${make_figure} -oT${output_name} -OP${options} -iT${report_title} -iX${job_description}
	    # zip results, add as an output.
        zip -r correlate_genomic_event_clinical.zip . -i *.pdf *.png
    }

	output {
        File R_file="analysis.results.RData"
        File zip_results="correlate_genomic_event_clinical.zip"
	}

	runtime {
		docker : "broadgdac/correlate_genomic_event_clinical:35"
	}

}

task correlate_genomic_events {
    
    File  cluster_label_file
    File  clinical_info_file 
    String? features_used # optional
    String filter_cluster_smaller_than
    String figure_height
    String figure_width 
    String survival_legend_position 
    String pvalue_cutoff 
    String qvalue_cutoff
    String make_figure
    String output_name 
    String options 
    String? report_title # optional 
    String? job_description # optional

	command {
		/R/RunR.sh -f main /src/Correlate_Genomic_Events.R /src -iD${cluster_label_file} -iC${clinical_info_file} -fF${features_used} -fC${filter_cluster_smaller_than} -fH${figure_height} -fW${figure_width} -fP${survival_legend_position} -cP${pvalue_cutoff} -cQ${qvalue_cutoff} -MF${make_figure} -oT${output_name} -OP${options} -iT${report_title} -iX${job_description}
	    # zip results, add as an output.
        zip -r correlate_genomic_events.zip . -i *.pdf *.png
    }

	output {
        File R_file="analysis.results.RData"
        File zip_results="correlate_genomic_events.zip"
	}

	runtime {
		docker : "broadgdac/correlate_genomic_events:35"
	}

}

task correlate_genomic_events_nozzle {
	#Inputs defined here
	File correlate_RData
	String? genevent_name

	command {
		/R/RunR.sh /src/GDAC_Correlate_Genomic_Events_NozzleReport.R /src  -iF${correlate_RData} ${"-gE" + genevent_name}
	}

	output {
		#Outputs defined here
		File nozzle_RData="nozzle.RData"
		File nozzle_html="nozzle.html"
	}

	runtime {
		docker : "broadgdac/correlate_genomic_events_nozzle:39"
	}

	meta {
		author : "Your name"
		email : "youremail@broadinstitute.org"
	}

}

task report_correlate_genomic_event_clinical {
	#Inputs defined here
	File correlate_RData
	String? genevent_name

	command {
		/R/RunR.sh /src/GDAC_Correlate_Genomic_Events_NozzleReport.R /src  -iF${correlate_RData} ${"-gE" + genevent_name}
	}

	output {
		#Outputs defined here
		File nozzle_RData="nozzle.RData"
		File nozzle_html="nozzle.html"
	}

	runtime {
		docker : "broadgdac/report_correlate_genomic_event_clinical:39"
	}

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	}

}

workflow correlate_genomic_events {
	call correlate_genomic_event_clinical 
	call report_correlate_genomic_event_clinical {
		input: correlate_RData=correlate_genomic_event_clinical.R_file
	}
}

