task correlate_expression_clinical {
	String? option_datatype
	File data_file_name
	File clinical_file_name
	String? genes_filtered_by_variance
	File? genes_filtered_by_list
	File? samples_filtered_by_list
	String? features_used
	String? output_file_name
	String? job_description

	command {
		Rscript /src/ClinicalAnalysisAllGenes.R --libdir /src ${"--option_datatype " + option_datatype} --data.file.name ${data_file_name} --clinical.file.name ${clinical_file_name} ${"--genes.filtered.by.variance " + genes_filtered_by_variance} ${"--genes.filtered.by.list " + genes_filtered_by_list} ${"--samples.filtered.by.list " + samples_filtered_by_list} ${"--features.used " + features_used} ${"--output.file.name " + output_file_name} ${"--description " + job_description}
	}

	output {
		File results_text="clinical-results.txt"
		File results_binary="clinical-results.RData"

	}

	runtime {
		docker : "broadgdac/correlate_expression_clinical:74"
	}

	meta {
		author : "Michael S. Noble"
		email : "mnoble@broadinstitute.org"
	}

}

task report_correlate_expression_clinical {
	#Inputs defined here
	File correlation_results		# Results from correlate_expression_clinical
	File rData
	Int max_table_rows
	Float cutoff_p
	Float cutoff_fdr
	String? report_type				# Mutation, miR, mRNA, RPPA, etc
	String? report_title
	String? job_description

	command {
		/R/RunR.sh /src/ClinicalAnalysisAllGenesReport_nozzle.R /src  -iR${correlation_results} -iI${rData} -rN${max_table_rows} -cP${cutoff_p} -cQ${cutoff_fdr} ${"-OP" + report_type} ${"-iT" + report_title} ${"-iX" + job_description}
		zip -r report.zip .
	}

	output {
		File ReportArchive="report.zip"

	}

	runtime {
		docker : "broadgdac/report_correlate_expression_clinical:99"
	}

	meta {
		author : "Michael S. Noble"
		email : "mnoble@broadinstitute.org"
	}

}

workflow correlate_clinical_vs_mrna {

	call correlate_expression_clinical
	call report_correlate_expression_clinical {
		input: correlation_results=correlate_expression_clinical.results_text,
		rData=correlate_expression_clinical.results_binary
	}

}

