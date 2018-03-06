task correlate_mir_mrna {
	#Inputs defined here
	File miR
	File mRNA
	String outputprefix
	Float mutual_info_threshold
	Float pearson_cor_threshold

	command {
		#Command goes here
		Rscript --slave --vanilla /src/CLR_all.R -m ${miR} -o ${outputprefix} -g ${mRNA} -t ${mutual_info_threshold} -p ${pearson_cor_threshold}
	}

	output {
		#Outputs defined here
		File correlation_significance="${outputprefix}.corr.sig.txt"
		File cg_summary="${outputprefix}.first.summary.gct"
		File summary="${outputprefix}.s.summary.gct"

	}

	runtime {
		docker : "broadgdac/correlate_mir_mrna:1"
	}

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	}

}

task postprocess_correlate_mir_mrna {
	#Inputs defined here
	String platform
	File mapfile
	File corrfile
	String outputprefix


	command {
		#Command goes here
		perl /src/GDAC_DirectTarget_v2.pl -a${platform} -p${mapfile} -c${corrfile} -o${outputprefix}
	}

	output {
		#Outputs defined here
		File mrna_topnodes="${outputprefix}.mrna.topnodes"
		File mrna_distr="${outputprefix}.mrna.distr"
		File corr3db="${outputprefix}.corr.3db.txt"
		File mirna_topnodes="${outputprefix}.mirna.topnodes"
		File mirna_sortednodes="${outputprefix}.mirna.sortednodes"
		File mrna_sortednodes="${outputprefix}.mrna.sortednodes"
		File corr="${outputprefix}.corr.txt"
		File summary="${outputprefix}.summary.gct"
		File mir_distr="${outputprefix}.mir.distr"
	}

	runtime {
		docker : "broadgdac/postprocess_correlate_mir_mrna:1"
	}

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	}

}

task report_mir_targets {
	#Inputs defined here
	File summary
	File cg_summary
	File dt_summary
	File mir_distr
	File mrna_distr
	File mir_topnodes
	File mrna_topnodes
	File dt_corr
	File corr_3db
	File mirna_sortednodes
	File mrna_sortednodes
	String platform


	command {
		#Command goes here
		Rscript /src/GDAC_miRDirectTargetsReport.R  -s ${summary} -r ${cg_summary} -d ${dt_summary} -m ${mir_distr} -g ${mrna_distr} -a ${mir_topnodes} -b ${mrna_topnodes} -w ${dt_corr} -x ${corr_3db} -c ${platform} -y ${mirna_sortednodes} -z ${mrna_sortednodes} && \
		zip -r nozzle_report.zip . -i *.PNG *.png *.html *.RData

	}

	output {
		# Nozzle output, as a zip file
		File NozzleZip="nozzle_report.zip"

	}

	runtime {
		docker : "broadgdac/report_mir_targets:1"
	}

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	}

}

workflow mir_find_direct_targets {
	String outputprefix
	String platform

	call correlate_mir_mrna {
		input: outputprefix=outputprefix
	}
	call postprocess_correlate_mir_mrna {
		input: corrfile=correlate_mir_mrna.correlation_significance, outputprefix=outputprefix, platform=platform
	}
	call report_mir_targets {
		input: summary=correlate_mir_mrna.summary, cg_summary=correlate_mir_mrna.cg_summary, dt_summary=postprocess_correlate_mir_mrna.summary,
		mir_distr=postprocess_correlate_mir_mrna.mir_distr, mrna_distr=postprocess_correlate_mir_mrna.mrna_distr,
		mir_topnodes=postprocess_correlate_mir_mrna.mirna_topnodes, mrna_topnodes=postprocess_correlate_mir_mrna.mrna_topnodes,
		dt_corr=postprocess_correlate_mir_mrna.corr, corr_3db=postprocess_correlate_mir_mrna.corr3db,
		mirna_sortednodes=postprocess_correlate_mir_mrna.mirna_sortednodes, mrna_sortednodes=postprocess_correlate_mir_mrna.mrna_sortednodes,
		platform=platform

	}

}
