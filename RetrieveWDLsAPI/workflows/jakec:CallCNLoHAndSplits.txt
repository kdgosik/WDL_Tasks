workflow GATK4AcnvCallCnLohAndBalancedSegmentsForCaptureWorkflow {
	call GATK4AcnvCallCnLohAndBalancedSegmentsForCapture
}

task GATK4AcnvCallCnLohAndBalancedSegmentsForCapture {
	String tumorSampleName
	File hetFile
	File segFile
	File tnFile
	Float rhoThreshold
	Int memoryGb
	Int diskSpaceGb

	command <<<
		java -Xmx4g -jar /gatk/gatk-protected.jar CallCNLoHAndSplits \
		--segments ${segFile} \
		--tumorHets ${hetFile} \
		--tangentNormalized ${tnFile} \
		--rhoThreshold ${rhoThreshold} \
		--outputDir .

		ls
	>>>

	output {
		File finalAcsSeg = "${tumorSampleName}-sim-final.acs.seg"
		File finalCnbSeg = "${tumorSampleName}-sim-final.cnb_called.seg"
		File finalCnvSeg = "${tumorSampleName}-sim-final.cnv.seg"
		File finalTitanHet = "${tumorSampleName}-sim-final.titan.het.tsv"
		File finalTitanTn = "${tumorSampleName}-sim-final.titan.tn.tsv"

	}

	runtime {
		docker: "jakeconway/gatk-protected:latest"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
	
}