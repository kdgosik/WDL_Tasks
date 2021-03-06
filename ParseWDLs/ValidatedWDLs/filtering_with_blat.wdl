
task calcDiskSizes {

	Float bamSize
    Float baiSize
    Float ponSize
    Float  mafSize
    
    command <<<
    echo -ne "(${bamSize}+${baiSize}+5000000000.0)/1000000000.0\n" | perl -ne 'print eval($_)."\n"' | grep -Po '^\d+'  > blat_size.dat
    echo -ne "(${bamSize}+${baiSize}+${ponSize}+5000000000.0)/1000000000.0\n" | perl -ne 'print eval($_)."\n"' | grep -Po '^\d+'> pon_size.dat
    echo -ne "(${bamSize}+${baiSize}+${ponSize}+15000000000.0)/1000000000.0\n" | perl -ne 'print eval($_)."\n"' | grep -Po '^\d+'> obf_size.dat
    echo -ne "(${mafSize}+1000000000.0)/1000000000.0\n" |  perl -ne 'print eval($_)."\n"' | grep -Po '^\d+'  > one_maf.dat
    echo -ne "(2.0*${mafSize}+1000000000.0)/1000000000.0\n" | perl -ne 'print eval($_)."\n"' | grep -Po '^\d+'  > two_maf.dat
    echo -ne "(3.0*${mafSize}+1000000000.0)/1000000000.0\n" | perl -ne 'print eval($_)."\n"' | grep -Po '^\d+'  > three_maf.dat
    echo -ne "(10.0*${mafSize}+1000000000.0)/1000000000.0\n" |  perl -ne 'print eval($_)."\n"' | grep -Po '^\d+'  > ten_maf.dat
    echo -ne "(25.0*${mafSize}+1000000000.0)/1000000000.0\n" |  perl -ne 'print eval($_)."\n"' | grep -Po '^\d+'  > 25_maf.dat

    >>>
    
    
  runtime {
    docker: "ubuntu:14.04"
    disks: "local-disk 1 HDD"
    preemptible: 1
  } 
  
  output {
  	Int ponDisk=read_int("pon_size.dat")
    Int blatDisk=read_int("blat_size.dat")
    Int obfDisk=read_int("obf_size.dat")
    Int oneMafGB=read_int("one_maf.dat")
    Int twoMafGB=read_int("two_maf.dat")
    Int threeMafGB=read_int("three_maf.dat")
    Int tenMafGB=read_int("ten_maf.dat")
    Int tfMafGB=read_int("25_maf.dat")
  }


	}



task blat {
  File bam
  File bai
  File maf
  String id
  Int diskGB

  command {
    python /opt/realign.py ${bam} ${maf} ${id}
}

  runtime {
    memory: "7 GB"
    disks: "local-disk ${diskGB} HDD" 
    docker: "gcr.io/broad-firecloud-itools/blat_filter_v2"
  }

  output {
    File blat_results = "${id}.blat.maf"
    File debug_results = "${id}.blat.rejected.maf"
    File blat_all_maf = "${id}.blat.all.maf"
  }
}



task MAFPonFilter {

	File MAFFile
	File PONFile
	File cytoBandFile
	File parameterFile
	Boolean useParameterFile
	String PairID
	String TOTNStr
	Int NMIN
	Float thresh
	Float WCUT
	String CODING_ONLY
	Int MIN_ALT_COUNT
    String filter_stub
    Int diskGB    

	command <<<

		#increase verbosity
		set -x

		echo "useParameterFile is ${useParameterFile}"

		#the cytoBand file should be in the expected location
		mkdir -pv /xchip/cga/reference/annotation/db/ucsc/hg19/
		cp -vf  ${cytoBandFile} /xchip/cga/reference/annotation/db/ucsc/hg19/cytoBand.txt

		#run the filter
		mkdir -v outdir
		if [ "${useParameterFile}" == "true" ] ; 
		then
			bash -c "source /matlab_source_file_2013a.sh && /usr/local/bin/maf_pon_filter \
			 ${MAFFile} ${PONFile} ${PairID} ${TOTNStr} ${parameterFile} ${NMIN} ${thresh} ${WCUT} \
			  ./outdir ${CODING_ONLY} ${MIN_ALT_COUNT} ${filter_stub} "  ; 
		else
			bash -c "source /matlab_source_file_2013a.sh && /usr/local/bin/maf_pon_filter \
			 ${MAFFile} ${PONFile} ${PairID} ${TOTNStr} . ${NMIN} ${thresh} ${WCUT} \
			  ./outdir ${CODING_ONLY} ${MIN_ALT_COUNT} ${filter_stub} "  ; 
		fi ;
        
		#capture output
		zip -r maf_pon_filter.out.zip ./outdir

	>>>


	output  {
			File mpOutzip="maf_pon_filter.out.zip"
		}

	runtime {
		docker: "broadinstitute/broadmutationcalling_filtering_beta:updated_maf_pon_filter_w_stub"
		memory: "24 GB"
		disks: "local-disk ${diskGB} SSD"
		}
	}

task OrientationBias_filter_Task {
	String ID
	File BAM
	File BAI
	File MAF
	File REFERENCE
	File DBSNP
	File DBSNPIDX
	String CONTEXT
	String ALTALLELE
	String STUB
	File detailMetrics
    Int diskGB
	
	command <<<
	#increase verbosity
	set -x

	#prep the MAF
    /usr/local/bin/maf_filter_prep.py ${MAF} > no_indels.maf


	#create output directory
	mkdir -v out

	METRICS_IN_LEN=`echo '${detailMetrics}'|tr -d "\n"|wc -c` ;
	echo "DM IS ${detailMetrics}"
	echo "METRICS_IN_LEN IS $METRICS_IN_LEN" ; 
	if [ "$METRICS_IN_LEN" -eq "0" ] ;
	then
		echo "Starting orientation bias filter from BAM ..." ;
		bash -x /usr/local/run_OrientationBias_filter.sh ${ID}  ${BAM} \
		 no_indels.maf ${REFERENCE} ${DBSNP} ${CONTEXT} ${ALTALLELE} ${STUB} /dev/null out ;
	else
		echo "Starting OrientationBiasFilter from passed-in pre_adapter_detail_metrics file ..." ;
		bash -x /usr/local/run_OrientationBias_filter.sh ${ID}  ${BAM} \
		 no_indels.maf ${REFERENCE} ${DBSNP} ${CONTEXT} ${ALTALLELE} ${STUB} ${detailMetrics} out ;
	fi ;

	#zip up all the output
	zip -r out_and_mat.zip ./out ./mat


	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_filtering_beta:1"
		memory: "7 GB"
		disks: "local-disk ${diskGB} HDD"
		}

	output  {
		#lane lists
		File detailMetricsOut="out/${ID}.pre_adapter_detail_metrics"
		File MatLabAndFilterDataZip="out_and_mat.zip"
		File filtered_maf="out/${ID}.OrientationBiasFilter.maf"
		File unfiltered_maf="out/${ID}.OrientationBiasFilter.unfiltered.maf"
		}
}



task maf_merge_task {

	#data files
	File oxog_all
	File oxog_pass
	File pon_one_zip
	File pon_two_zip
	File blat_maf_pass
	File blat_maf_all

	String id
	Int diskGB
	String memGB
	File inputMaf


	command <<<

	#increase verbosity
	set -x

	#unzip any zip files
	unzip -d mp1_out ${pon_one_zip}
	unzip -d mp2_out ${pon_two_zip}

	#remove any line starting with #
	for M in `find *_out oxog_unzip   ${inputMaf} ${blat_maf_pass} ${blat_maf_all} | grep -i '\.maf$'   `; do
		echo "Removing lines starting with # from $M" ;
		grep -Pv '^#' $M > tmp.maf ;
		mv -vf tmp.maf $M ;
	done ;

	#make UNION filtered from all data 
	#output contains every row of original input AND has additional columns with each additional column coming from one of the filters
	#1) BLAT
	#2) oxoG
	#3) MP1
	#4) MP2
	OXOG_UNFILTERED=${oxog_all}
	echo "OXOG_UNFILTERED IS $OXOG_UNFILTERED" ; 
	MP1_ALL=`find mp1_out |grep -Pi '\.maf$'|grep -Piv '\.pass\.'` ; 
	echo "MP1_ALL IS $MP1_ALL" ; 
	MP2_ALL=`find mp2_out |grep -Pi '\.maf$'|grep -Piv '\.pass\.'`
	echo "MP2_ALL IS $MP2_ALL" ; 
	python /usr/local/bin/maf_maf_merge.py ${inputMaf} ${blat_maf_all} with_blat.maf && \
	 python /usr/local/bin/maf_maf_merge.py with_blat.maf $OXOG_UNFILTERED blat.oxog.maf && \
	 python /usr/local/bin/maf_maf_merge.py blat.oxog.maf $MP1_ALL blat.oxog.mp1.maf && \
	 python /usr/local/bin/maf_maf_merge.py blat.oxog.mp1.maf $MP2_ALL blat.oxog.mp1.mp2.maf && \
	 mv -vf blat.oxog.mp1.mp2.maf ${id}.filter.union.maf
	FIRST_RES=$? ;


	if [ "$FIRST_RES" -eq "0" ] ; 
		then
		#make INTERSECTION filtered
		#output contains rows such that each row passes ALL the filters
		#make union_filtered from all data 
		#output contains every row of original input AND has additional columns with each additional column coming from one of the filters
		#1) BLAT
		#2) oxoG
		#3) MP1
		#4) MP2
		OXOG_FILTERED=${oxog_pass}
		echo "OXOG_FILTERED IS $OXOG_FILTERED" ;
		MP1_PASS=`find mp1_out | grep -Pi '\.maf$'|grep -Pi '\.pass\.'`
		echo "MP1_PASS IS $MP1_PASS" ; 
		MP2_PASS=`find mp2_out | grep -Pi '\.maf$'|grep -Pi '\.pass\.'`
		echo "MP2_PASS IS $MP2_PASS" ; 
		python /usr/local/bin/maf_maf_merge.py -i ${inputMaf} ${blat_maf_pass} with_blat.i.maf && \
		 python /usr/local/bin/maf_maf_merge.py -i with_blat.i.maf $OXOG_FILTERED blat.oxog.i.maf && \
		 python /usr/local/bin/maf_maf_merge.py -i blat.oxog.i.maf $MP1_PASS blat.oxog.mp1.i.maf && \
		 python /usr/local/bin/maf_maf_merge.py -i blat.oxog.mp1.i.maf $MP2_PASS blat.oxog.mp1.mp2.i.maf && \
		 mv -vf blat.oxog.mp1.mp2.i.maf ${id}.filter.intersection.maf
	else
		echo "Error !" ;
		/bin/bash -c "exit 2;" 
	fi ;


	>>>

	output 
		{
		File union_maf="${id}.filter.union.maf"
		File intersection_maf="${id}.filter.intersection.maf" 
		}

	runtime {
		docker : "broadinstitute/blat_filtering:eddie_end_filter_wf_merger"
		disks: "local-disk ${diskGB} HDD"
		memory: "${memGB}GB"
		}

	}


task mergeInIndelsToOBF
	{

	File obfFilteredMaf
	File obfUnfilteredMaf
	File origMaf
	Int diskSize
	String id

	command <<<

		#validate column headers
		#Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2
		VALID_HEADERS=`	grep -P '^Hugo'  ${origMaf} | head -1 | awk '{ if($11=="Reference_Allele" && $12=="Tumor_Seq_Allele1" && $13=="Tumor_Seq_Allele2") { print "go"  }}'`
		if [ "$VALID_HEADERS" == "go" ] ; then 
			#extract indels from input maf
			grep -P '^Hugo' ${origMaf} > ${id}.indels.maf ;
			cat ${origMaf} |grep -Pv '^#'|grep -Pv '^Hugo'  |awk '{if($11~/\-/ || $12~/\-/ || $13~/\-/) { print $0 } }' >> ${id}.indels.maf ; 

			#merge back indels into OBF output
			python /usr/local/bin/tsvConcatFiles.py  ${id}.indels.maf ${obfFilteredMaf}   --outputFilename=${id}.filt.indels.maf ;
			python /usr/local/bin/tsvConcatFiles.py  ${id}.indels.maf ${obfUnfilteredMaf} --outputFilename=${id}.unfilt.indels.maf ;
		fi ;


		>>>

	runtime {
		docker : "eddiebroad/tsvconcatfiles"
		disks: "local-disk ${diskSize} HDD"
		memory: "0.25 GB"
		}

	output
		{
		File indels="${id}.indels.maf"
		File filtIndels="${id}.filt.indels.maf"
		File unFiltIndels="${id}.unfilt.indels.maf"
		}
	}




workflow FilterWorkflow  {

	Int preemptible
	File oncoDBTarBall
	File inMAF
	File cytoBandFile
	File parameterFile
	String PairID
	String TOTNStr
	Int NMIN
	Float thresh
	Float WCUT
	String CODING_ONLY
	Int MIN_ALT_COUNT
	Boolean useParameterFile
	String TUM_ID
	File REFERENCE
	File DBSNP
	File DBSNPIDX	
	File tumorBam
	File tumorBamIdx
    File metricsFile
    File ponFileOne
    File ponFileTwo

    call  calcDiskSizes {
    input:

	 bamSize=size(tumorBam),
     baiSize=size(tumorBamIdx),
     ponSize=size(ponFileOne),
     mafSize=size(inMAF)
    }

	call blat {
		input:
			bam=tumorBam,
			bai=tumorBamIdx,
			maf=inMAF,
			id=TUM_ID,
            diskGB=calcDiskSizes.blatDisk
		}



	call MAFPonFilter as mp1 {
		input : 
        	PONFile=ponFileOne,
			MAFFile=inMAF,
			cytoBandFile=cytoBandFile,
			useParameterFile=useParameterFile,
			parameterFile=parameterFile,
			PairID=PairID,
			TOTNStr=TOTNStr,
			NMIN=NMIN,
			thresh=thresh,
			WCUT=WCUT,
			CODING_ONLY=CODING_ONLY,
			MIN_ALT_COUNT=MIN_ALT_COUNT,
            diskGB=calcDiskSizes.ponDisk
		}

	call MAFPonFilter as mp2 {
		input : 
           	PONFile=ponFileTwo,
			MAFFile=inMAF,
			cytoBandFile=cytoBandFile,
			useParameterFile=useParameterFile,
			parameterFile=parameterFile,
			PairID=PairID,
			TOTNStr=TOTNStr,
			NMIN=NMIN,
			thresh=thresh,
			WCUT=WCUT,
			CODING_ONLY=CODING_ONLY,
			MIN_ALT_COUNT=MIN_ALT_COUNT,
            diskGB=calcDiskSizes.ponDisk
		}

	call OrientationBias_filter_Task as oxoGOBF {
		input :
			detailMetrics=metricsFile,
			ID=TUM_ID,
			BAM=tumorBam,
			BAI=tumorBamIdx,
			MAF=inMAF,
			REFERENCE=REFERENCE,
			DBSNP=DBSNP,
			DBSNPIDX=DBSNPIDX,
            diskGB=calcDiskSizes.obfDisk
		}

	call mergeInIndelsToOBF {
    	input:
              obfFilteredMaf=oxoGOBF.filtered_maf,
              obfUnfilteredMaf=oxoGOBF.unfiltered_maf,
              origMaf=inMAF,
              diskSize=calcDiskSizes.threeMafGB,
              id=TUM_ID    
	    	}

	call maf_merge_task {
		input:
			oxog_all=mergeInIndelsToOBF.unFiltIndels,
			oxog_pass=mergeInIndelsToOBF.filtIndels,
			pon_one_zip=mp1.mpOutzip,
			pon_two_zip=mp2.mpOutzip,
			blat_maf_pass=blat.blat_results,
			blat_maf_all=blat.blat_all_maf,
			id=TUM_ID,
			diskGB=calcDiskSizes.tfMafGB,
			memGB=4,
			inputMaf=inMAF
		}


	output {
		maf_merge_task.union_maf
		maf_merge_task.intersection_maf
		}

	}