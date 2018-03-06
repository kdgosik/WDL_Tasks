workflow oncotateGatkCnvFileForCaptureWorkflow {
	call oncotateGatkCnvFileForCapture
}

task oncotateGatkCnvFileForCapture {
	String txMode
	String genomeBuild
	String sampleName
	String inputFormat
	String outputFormat
	File dataSourceTarGzFile
	String dataSourceDir = basename(dataSourceTarGzFile, ".tar.gz")
	File defaultConfig
	File calledSegmentsFile
	Int memoryGb
	Int diskSpaceGb

	command <<<
		
		cp -r /root/oncotator_venv/ $PWD
		ls -a

		tar -xzvf ${dataSourceTarGzFile}

		./oncotator_venv/bin/oncotator -v -i ${inputFormat} \
		-o ${outputFormat} \
		--db-dir ${dataSourceDir}/ \
		--no-multicore \
		--default_config ${defaultConfig} ${calledSegmentsFile} ${sampleName}.called.seg.annotated ${genomeBuild} \
		--tx-mode ${txMode}

	>>>

	output {
		File annotatedSegFile = "${sampleName}.called.seg.annotated"
	}

	runtime {
		docker: "broadinstitute/oncotator:1.9.3.0"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}