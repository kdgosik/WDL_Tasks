task oncotator {
  File oncoDBTarBall
  File IN
  String OUT_TYPE
  String IN_TYPE
  String OTHER_FLAGS
  String id

  command {
	#obtain the name of the directory for oncodb
	ONCO_DB_DIR_NAME=`gunzip -c ${oncoDBTarBall} |tar -tf /dev/stdin|head -1` ; 

	#unpack the oncodir! (may take a while...)
	tar -xzf ${oncoDBTarBall}
	OUTNAME="${id}.${OUT_TYPE}"

	if [ "${OUT_TYPE}" = "MAF" ]
	then
	 OUT_T="TCGAMAF"
	else
	 OUT_T="${OUT_TYPE}"
	fi
	if [ "${IN_TYPE}" = "MAF" ]
	then
	 IN_T="TCGAMAF"
	else
	 IN_T="${IN_TYPE}"
	fi

	#Run the merged filtered VCF (from both mutects through Oncotator) 
	/root/oncotator_venv/bin/Oncotator -i $IN_T --db-dir `pwd`/$ONCO_DB_DIR_NAME -o $OUT_T ${OTHER_FLAGS} ${IN} $OUTNAME hg19 

}

  runtime {
    memory: "7 GB"
    docker: "gcr.io/broad-firecloud-itools/oncotator"
    disks: "local-disk 100 HDD"
    preemptible: 0
  }

  output {
    File oncotator_out = "${id}.${OUT_TYPE}"
  }
}

workflow oncotator_wkflw {
  call oncotator
}