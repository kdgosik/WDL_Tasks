task acnv_make_pon
	{

	File jar_fn
	String memGB
	String diskGB
	String pSetID
	String cpus
	Array[File] pcovs

	command <<<

	#increase verbosity
	set -x

	for F in ${sep =" " pcovs}; do 
		echo "writing $F to agg.pcovs" ;
		echo $F >> agg.pcovs ;
		done ;

	java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar ${jar_fn} CombineReadCounts --inputList agg.pcovs -O merged.pcovs -MOF 200 

	java -Xmx16g -Djava.library.path=/usr/lib/jni/ -jar ${jar_fn} CreatePanelOfNormals -I merged.pcovs  -O ${pSetID}.acnv.pon

	>>>

	output {
		File acnv_pon="${pSetID}.acnv.pon"
		}

	runtime {
		disks: "local-disk ${diskGB} HDD"
		memory: "${memGB} GB"
		docker: "dlivitzbroad/cnvvdev"
		preemptible: 3
		}

	}

workflow acnv_make_pon_workflow {
	
	call acnv_make_pon {}

}