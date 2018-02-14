workflow OncotateSnpForCaptureWorkflow {
	call OncotateSnpForCapture
}

task OncotateSnpForCapture {
	String txMode
	String genomeBuild
	String pairName
	String inputFormat
	String outputFormat
	File dataSourceTarGzFile
	String dataSourceDir = basename(dataSourceTarGzFile, ".tar.gz")
	File defaultConfig
	File validatedMafliteFile
	Int memoryGb
	Int diskSpaceGb

	command <<<
		cp -r /root/oncotator_venv/ $PWD
		cp /root/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt $PWD
		ls -a

		tar -xzvf ${dataSourceTarGzFile}

		./oncotator_venv/bin/oncotator -v -i ${inputFormat} \
		-o ${outputFormat} \
		--db-dir ${dataSourceDir}/ \
		--no-multicore \
		--default_config ${defaultConfig} ${validatedMafliteFile} ${pairName}.snp.capture.maf.annotated ${genomeBuild} \
		--tx-mode ${txMode} \
		--log_name oncotator_firecloud.log --prepend --infer-onps -c tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt --collapse-number-annotations
	>>>

	output {
		File oncotatedMaf = "${pairName}.snp.capture.maf.annotated"
	}

	runtime {
		docker: "broadinstitute/oncotator:1.9.3.0"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}