workflow callStatsToMafLiteForCaptureWorkflow {
	call callStatsToMafLiteForCapture
}

task callStatsToMafLiteForCapture {
	File callStatsFile
	String pairName
	String triallelic_mode_KEEP_or_REJECT
	String mode
	Int f_threshold
	Int genomeBuild
	Int memoryGb
	Int diskSpaceGb

	command <<<
		perl /usr/gitc/call_stats_to_maflite.pl ${callStatsFile} ${genomeBuild} ${mode} ${f_threshold} \
		${triallelic_mode_KEEP_or_REJECT} ${pairName}.maf tumor_f,init_t_lod,t_lod_fstar,t_alt_count,t_ref_count,judgement
	>>>

	output {
		File outputMaf = "${pairName}.maf"
	}

	runtime {
		docker: "jakeconway/somatic_var_calling:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}