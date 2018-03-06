workflow revertSamTest {
	call revertSam
}

task revertSam {
	File inputBam
	String sampleName
	Int memoryGb
    Int diskSpaceGb

	command <<<
	java -Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 \
	-XX:GCTimeLimit=50 \
	-XX:GCHeapFreeLimit=10 \
	-Xmx4000m -jar /usr/gitc/picard.jar RevertSam \
	TMP_DIR=. \
	VALIDATION_STRINGENCY=SILENT \
	OUTPUT=${sampleName}.reverted.bam \
	INPUT=${inputBam} \
	SORT_ORDER=queryname \
	RESTORE_ORIGINAL_QUALITIES=true \
	REMOVE_DUPLICATE_INFORMATION=true \
	REMOVE_ALIGNMENT_INFORMATION=true \
	SANITIZE=true \
	MAX_DISCARD_FRACTION=0.01 \
	ATTRIBUTE_TO_CLEAR=X0 \
	ATTRIBUTE_TO_CLEAR=X1 \
	ATTRIBUTE_TO_CLEAR=XA \
	ATTRIBUTE_TO_CLEAR=XC \
	ATTRIBUTE_TO_CLEAR=XG \
	ATTRIBUTE_TO_CLEAR=XM \
	ATTRIBUTE_TO_CLEAR=XN \
	ATTRIBUTE_TO_CLEAR=XO \
	ATTRIBUTE_TO_CLEAR=XT \
	ATTRIBUTE_TO_CLEAR=AM \
	ATTRIBUTE_TO_CLEAR=AS \
	ATTRIBUTE_TO_CLEAR=BQ \
	ATTRIBUTE_TO_CLEAR=CC \
	ATTRIBUTE_TO_CLEAR=CP \
	ATTRIBUTE_TO_CLEAR=E2 \
	ATTRIBUTE_TO_CLEAR=H0 \
	ATTRIBUTE_TO_CLEAR=H1 \
	ATTRIBUTE_TO_CLEAR=H2 \
	ATTRIBUTE_TO_CLEAR=HI \
	ATTRIBUTE_TO_CLEAR=IH \
	ATTRIBUTE_TO_CLEAR=MF \
	ATTRIBUTE_TO_CLEAR=NH \
	ATTRIBUTE_TO_CLEAR=OC \
	ATTRIBUTE_TO_CLEAR=OP \
	ATTRIBUTE_TO_CLEAR=PQ \
	ATTRIBUTE_TO_CLEAR=R2 \
	ATTRIBUTE_TO_CLEAR=S2 \
	ATTRIBUTE_TO_CLEAR=SM \
	ATTRIBUTE_TO_CLEAR=SQ \
	ATTRIBUTE_TO_CLEAR=U2 \
	ATTRIBUTE_TO_CLEAR=XQ
	>>>

	output {
		File revertedBam = "${sampleName}.reverted.bam"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}