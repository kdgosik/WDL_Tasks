task make_acnv_pon
	{

	String memGB
	String diskGB
	String pSetID
    String cpus
	Array[File] pcovs

	command <<<

	#increase verbosity
	set -x

	for F in ${sep =" " pcovs} do ;
		echo "writing $F to agg.pcovs"
		echo $F >> agg.pcovs
	done ;



	java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar /root/gatk-protected.jar  CombineReadCounts \
	 	--inputList agg.pcovs \
    	-O merged.pcovs -MOF 200 


	java -Xmx16g -Djava.library.path=/usr/lib/jni/ -jar /root/gatk-protected.jar  CreatePanelOfNormals \
	 	-I merged.pcovs  \
       	-O ${pSetID}.acnv.pon

	>>>

	output {
		File acnv_pon="${pSetID}.acnv.pon"
		}

	runtime {
    	disks: "local-disk ${diskGB} HDD"
    	memory: "${memGB} GB"
    	docker: "broadinstitute/hellbender_pon"
   	 	preemptible: 3
		}

	}

workflow make_acnv_pon_workflow {
	String memGB
	String diskGB
	String pSetID
    String cpus
	Array[File] pcovs
	
	call make_acnv_pon {
	input:
	diskGB=diskGB,
	memGB=memGB,
	pSetID=pSetID,
	pcovs=pcovs

}
}