task MAFPonFilter {

	#Tool inputs
	File MAFFile
	File PONFile
	File cytoBandFile
	String PairID
	String TOTNStr
	Int NMIN
	Float thresh
	Float WCUT
	Int CODING_ONLY
	Int MIN_ALT_COUNT
	Int diskGB

	parameter_meta {
	  MAFFile : "input MAF for PON filter analysis"
	  PONFile : "formatted panel-of-normals file"
	  cytoBandFile : "TSV file of chromsomal annotation ; chr, start, end, band, stain"
	  PairID : "used to name the output files and title other output"
	  TOTNStr : "indicating pair status : can be 'TO' or 'TN'"
	  NMIN : "min number of supporting reads"
	  thresh : "threshold for pass/fail with log likelihood"
	  WCUT : "threshold for pass/fail with pon-computed weight"
	  CODING_ONLY : "analyze coding regions only?"
	  MIN_ALT_COUNT : "minimum alt allele count from PoN to consider for filtering"
	}



	command <<<

		#increase verbosity
		set -x

		#the cytoBand file should be in the expected location
		mkdir -pv /xchip/cga/reference/annotation/db/ucsc/hg19/
		cp -v  ${cytoBandFile} /xchip/cga/reference/annotation/db/ucsc/hg19/cytoBand.txt

		#run the filter
		#Note "." is used for the parameter file whose usage is *not* exposed/functional in this WDL
		mkdir -v outdir
		bash -c "source /matlab_source_file_2013a.sh && /usr/local/bin/maf_pon_filter \
		 ${MAFFile} ${PONFile} ${PairID} ${TOTNStr} . ${NMIN} ${thresh} ${WCUT} \
		  ./outdir ${CODING_ONLY} ${MIN_ALT_COUNT}" ; 
        
		#capture output
		zip -r maf_pon_filter.out.zip ./outdir



	>>>


	output  {
			File mpOutzip="maf_pon_filter.out.zip"
			File allMaf="outdir/${PairID}.pon_annotated.maf"
			File passMaf="outdir/${PairID}.pon_annotated.pass.maf"
		}

	runtime {
		docker: "broadinstitute/broadmutationcalling_filtering_beta@sha256:d2df6d9d705e90d3ee926b72a3edec5034dd5bdd64c0cf1cabd9fc8462702c79"
		memory: "7 GB"
		disks: "local-disk ${diskGB} HDD"
		}
	}

task OrientationBias_filter_Task {

	#tool inputs
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
	File? detailMetrics
	Int diskGB


  parameter_meta {
    ID : "String prefix of the output"
    BAM: "the tumor bam file"
    BAM_IDX : "its .bai index file"
    detailMetrics : "output from the generateMetrics_task from picard CollectSequencingArtifactMetrics run with the settings here ; this allows for passthrough instead of recomputing"
    DBSNP : "VCF file of DB SNP variants"
    DBSNPIDX :"its index file"
    REFERENCE : "fasta used by the BAM and the DBSNP"
    CONTEXT : "three bases before, of, and after the effect"
    ALTALLELE : "the alternate allel of the effect"
    STUB : "string used to indicate in the output the effect name"
    detailMetrics : "output from the generateMetrics_task from picard CollectSequencingArtifactMetrics run with the settings here ; this allows for passthrough instead of recomputing"
	}

	
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

    #add indels back 
    VALID_HEADERS=`	grep -P '^Hugo'  ${MAF} | head -1 | awk '{ if($11=="Reference_Allele" && $12=="Tumor_Seq_Allele1" && $13=="Tumor_Seq_Allele2") { print "go"  }}'`
		if [ "$VALID_HEADERS" == "go" ] ; then 
			#extract indels from input maf
			grep -P '^Hugo' ${MAF} > ${ID}.indels.maf ;
			cat ${MAF} |grep -Pv '^#'|grep -Pv '^Hugo'  |awk '{if($11~/\-/ || $12~/\-/ || $13~/\-/) { print $0 } }' >> ${ID}.indels.maf ; 

			#merge back indels into OBF output
			python /usr/local/bin/tsvConcatFiles.py  ${ID}.indels.maf out/${ID}.OrientationBiasFilter.maf            --outputFilename=${ID}.filt.indels.maf   && mv -vf ${ID}.filt.indels.maf   out/${ID}.OrientationBiasFilter.maf 
			python /usr/local/bin/tsvConcatFiles.py  ${ID}.indels.maf out/${ID}.OrientationBiasFilter.unfiltered.maf --outputFilename=${ID}.unfilt.indels.maf && mv -vf ${ID}.unfilt.indels.maf out/${ID}.OrientationBiasFilter.unfiltered.maf ;
		fi ;



	#zip up all the output
	zip -r out_and_mat.zip ./out ./mat

	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_filtering_beta@sha256:d2df6d9d705e90d3ee926b72a3edec5034dd5bdd64c0cf1cabd9fc8462702c79"
		memory: "7 GB"
		disks: "local-disk ${diskGB} HDD"
		}

	output  {
		#lane lists
		Float q_val=read_float("out/${ID}.orientation_BiasQ.txt")
		File detailMetricsOut="out/${ID}.pre_adapter_detail_metrics"
		File MatLabAndFilterDataZip="out_and_mat.zip"
		File filtered_maf="out/${ID}.OrientationBiasFilter.maf"
		File unfiltered_maf="out/${ID}.OrientationBiasFilter.unfiltered.maf"
		}
}





task generateMetrics_task {

	#tool inputs
	String ID
	File BAM
	File BAM_IDX
	File? detailMetrics
	File DBSNP
	File DBSNPIDX
	File REFERENCE
	Int diskGB


  parameter_meta {
    ID : "String prefix of the output"
    BAM: "the tumor bam file"
    BAM_IDX : "its .bai index file"
    detailMetrics : "output from the generateMetrics_task from picard CollectSequencingArtifactMetrics run with the settings here ; this allows for passthrough instead of recomputing"
    DBSNP : "VCF file of DB SNP variants"
    DBSNPIDX :"its index file"
    REFERENCE : "fasta used by the BAM and the DBSNP"
	}


	command <<<

	#increase verbosity
	set -x

	#see if a metrics file was passed in and create one if not
	METRICS_IN_LEN=`echo '${detailMetrics}'|tr -d "\n"|wc -c` ;
	if [ "$METRICS_IN_LEN" -eq "0" ] ;
	then 
		#file not found, create it!
		/usr/local/bin/java  -Xmx3600M  -jar /usr/local/CollectSequencingArtifactMetrics/picard.1.895.jar CollectSequencingArtifactMetrics \
		  DB_SNP=${DBSNP} INPUT=${BAM} OUTPUT=${ID}  REFERENCE_SEQUENCE=${REFERENCE} \
		  MINIMUM_QUALITY_SCORE=20 MINIMUM_MAPPING_QUALITY=30 MINIMUM_INSERT_SIZE=60 MAXIMUM_INSERT_SIZE=600  \
		  INCLUDE_UNPAIRED=false TANDEM_READS=false USE_OQ=true CONTEXT_SIZE=1 ASSUME_SORTED=true \
		  STOP_AFTER=100000000 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 \
		  MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json
	else
		#file found, pass it on!
		echo "Passing the metrics file along!"
		ln -v ${detailMetrics} ${ID}.pre_adapter_detail_metrics
	fi ;

	>>>

	output {
		File metricsFile="${ID}.pre_adapter_detail_metrics"
		}

	runtime {
		docker : "broadinstitute/broadmutationcalling_filtering_beta@sha256:d2df6d9d705e90d3ee926b72a3edec5034dd5bdd64c0cf1cabd9fc8462702c79"
		memory: "3 GB"
		disks: "local-disk ${diskGB} HDD"
	}
}


task VCF_to_MAF_task {

	#tool inputs
	File inputVCF
	File oncoDBTarBall
	String pairName
	String caseName
	String ctrlName
	Int diskGB


  parameter_meta {
    inputVCF : "un-onco-tated VCF to be converted to MAF by oncotator"
    oncoDBTarBall : ".tar or .tgz or .tar.gz file containing an oncotator DB directory with just hg19"
	}



	command <<<

		#increase verbosity
		set -x

		#find TARBALL type
		TYPE=`echo 'if("${oncoDBTarBall}"=~m/z$/) { print "GZ" ; } else { print "TAR" ; } '| perl` ; 

		#obtain the name of the directory for oncodb
		#and unpack based on the TYPE
		if [ "$TYPE" == "GZ" ] ; then
			ONCO_DB_DIR_NAME=`gunzip -c ${oncoDBTarBall} |tar -tf /dev/stdin|head -1` ; 
			echo "ONCO_DB_DIR_NAME is $ONCO_DB_DIR_NAME" ; 
			tar -xzf ${oncoDBTarBall}
		else
			ONCO_DB_DIR_NAME=`tar -tf ${oncoDBTarBall} |head -1` ;
			tar -xf ${oncoDBTarBall} ;
		fi ;
		
		#Run the VCF through oncotator for VCF->MAF conversion
		/root/oncotator_venv/bin/Oncotator -i VCF   --skip-no-alt  \
            --db-dir `pwd`/$ONCO_DB_DIR_NAME ${inputVCF} output.maf hg19 ;

#add in the tumor/normal sampleIDs
python_prog_str="
import sys
def main():
	input_filename=str(sys.argv[1])
	output_filename=str(sys.argv[2])
	case_name=str(sys.argv[3])
	ctrl_name=str(sys.argv[4])
	reader=open(input_filename,'r')
	writer=open(output_filename,'w')
	for line in reader:
		if(line.startswith('#')):
			writer.write(line)
		elif(line.startswith('Hugo_Symbol')):
			writer.write(line)
		else:
			pieces=line.split('\t')
			#columns 16 and 17 (tumor/normal barcodes) are at indices 15 and 16
			pieces[15]=case_name
			pieces[16]=ctrl_name
			writer.write('\t'.join(pieces))
	writer.close()
	reader.close()
if __name__ == '__main__':
	main()
"
		/root/oncotator_venv/bin/python -c "$python_prog_str" output.maf  output.with_bcs.maf ${caseName} ${ctrlName}


        #add the pair name to the output
        cat output.with_bcs.maf |grep -P '^#' > ${pairName}.maf #header
        cat output.with_bcs.maf |awk '{if($1=="Hugo_Symbol" ) { print $0 "\tPairID" }}' >> ${pairName}.maf
        cat output.with_bcs.maf |grep -Pv '^#'|grep -Pv '^Hugo_Symbol'|awk '{print $0 "\t${pairName}"}' >> ${pairName}.maf



		>>>

	runtime {
		memory: "3 GB"
		docker : "broadinstitute/oncotator@sha256:434b232606f5f7df02a3a64be46ce2214308c9cdad577759524b9e68c4d70979"
		disks: "local-disk ${diskGB} HDD"
		}

	output {
		File rawOnco="output.maf"
		File outputMAF="output.with_bcs.maf"
		File outputPairMaf="${pairName}.maf"
		}


	}



task pf_counter_task {

	File allMaf
	File passMaf
	Int diskGB


	command <<<

	#increase verbosity
	set -x

	#get counts
	cat ${allMaf}  | grep -Pv '^#'|grep -Piv '^Hugo'|wc -l > all_count.dat
	cat ${passMaf} | grep -Pv '^#'|grep -Piv '^Hugo'|wc -l > pass_count.dat

	>>>

	runtime {
		docker : "broadinstitute/broadmutationcalling_filtering_beta@sha256:d2df6d9d705e90d3ee926b72a3edec5034dd5bdd64c0cf1cabd9fc8462702c79"
		disks  : "local-disk ${diskGB} HDD"
		memory: "0.1 GB"
		}

	output {
		Int all_count=read_int("all_count.dat")
		Int pass_count=read_int("pass_count.dat")
		}

}



task maf_merge_task {

	#data files
	File oxog_zip
	File ffpe_zip
	File pon_zip
	File in_vcf
	File in_maf
	String PairID
	Int diskGB

  parameter_meta {
    oxog_zip : "output zip from oxoG orientation bias filter"
    ffpe_zip : "output zip from ffpe orientation bias filter"
    pon_zip : "output from MAF pon filter"
    in_vcf : "the same input to VCF_to_MAF_task which will be annotated with filter information then output"
	}


	command <<<

	#increase verbosity
	set -x

	#unzip each file and obtain inputs for merge script
	unzip -d oxog_out ${oxog_zip}
	OXOG_MAF=`find $PWD/oxog_out | grep -P 'OrientationBiasFilter.unfiltered.maf$'` ;
	OXOG_MAF_PASS=`find $PWD/oxog_out | grep -P 'OrientationBiasFilter.maf$'`
	unzip -d ffpe_out ${ffpe_zip}
	FFPE_MAF=`find $PWD/ffpe_out | grep -P 'OrientationBiasFilter.unfiltered.maf$'` ; 
	FFPE_MAF_PASS=`find $PWD/ffpe_out | grep -P 'OrientationBiasFilter.maf$'` ; 
	unzip -d pon_out ${pon_zip}
	PON_MAF=`find $PWD/pon_out | grep -P 'pon_annotated.maf$'`
	PON_PASS_MAF=`find $PWD/pon_out | grep -P 'pon_annotated.pass.maf$'` ; 

	#define output
	OUTPUT_FILE=${PairID}.OrientationBiasAndPON_filtered.vcf ; 
	
	#run the merger
	/usr/local/bin/maf_filters_vcf_merge.py ${in_vcf} $OXOG_MAF $FFPE_MAF $PON_MAF $PON_PASS_MAF $OUTPUT_FILE

	#merge as both union and intersection using maf_maf_merge.py
 
	#make union MAF from filters (note that '${in_maf}'' is passed second during first call because
	#when conflicts arise at a cell, the value from the first file passed to maf_maf_merge is kept')
	python /usr/local/bin/maf_maf_merge.py $OXOG_MAF ${in_maf}  all_oxog.maf && \
	python /usr/local/bin/maf_maf_merge.py all_oxog.maf $FFPE_MAF all_oxog_ffpe.maf && \
	python /usr/local/bin/maf_maf_merge.py all_oxog_ffpe.maf $PON_MAF ${PairID}.union.maf ; 

	#make intersection MAF from filters
	python /usr/local/bin/maf_maf_merge.py -intersection  $OXOG_MAF_PASS  ${in_maf} i_all_oxog.maf
	python /usr/local/bin/maf_maf_merge.py -intersection i_all_oxog.maf $FFPE_MAF_PASS i_all_oxog_ffpe.maf
	python /usr/local/bin/maf_maf_merge.py -intersection i_all_oxog_ffpe.maf $PON_PASS_MAF ${PairID}.intersection.maf


	>>>


	runtime {
		docker : "broadinstitute/broadmutationcalling_filtering_beta@sha256:d2df6d9d705e90d3ee926b72a3edec5034dd5bdd64c0cf1cabd9fc8462702c79"
		disks: "local-disk ${diskGB} HDD"
		}

	output {
		#output VCF
		File merged_filtered_vcf="${PairID}.OrientationBiasAndPON_filtered.vcf"
		File union_maf="${PairID}.union.maf"
		File intersection_maf="${PairID}.intersection.maf"
		}


	}


task oncotation_task {

	File oncotationTarBall
	File inVCF
	String PairID
	Int diskGB

	command <<<

	#increase verbosity
	set -x

	#find TARBALL type
	TYPE=`echo 'if("${oncotationTarBall}"=~m/z$/) { print "GZ" ; } else { print "TAR" ; } '| perl` ; 

	#obtain the name of the directory for oncodb
	#and unpack based on the TYPE
	if [ "$TYPE" == "GZ" ] ; then
		ONCO_DB_DIR_NAME=`gunzip -c ${oncotationTarBall} |tar -tf /dev/stdin|head -1` ; 
		echo "ONCO_DB_DIR_NAME is $ONCO_DB_DIR_NAME" ; 
		tar -xzf ${oncotationTarBall}
	else
		ONCO_DB_DIR_NAME=`tar -tf ${oncotationTarBall} |head -1` ;
		tar -xf ${oncotationTarBall} ;
	fi ;

	#oncotate!
	/root/oncotator_venv/bin/Oncotator -i VCF  \
		--db-dir `pwd`/$ONCO_DB_DIR_NAME -o VCF  ${inVCF} ${PairID}.filtered.oncotated.vcf hg19 ;

	>>>	

	runtime {
		docker : "broadinstitute/oncotator@sha256:434b232606f5f7df02a3a64be46ce2214308c9cdad577759524b9e68c4d70979"
		disks: "local-disk ${diskGB} HDD"
		}

	output {
		File annotatedFilteredVCF="${PairID}.filtered.oncotated.vcf"
		}

	}



task sizesCalc {

	Float bamSize
	Float baiSize
	Float dbSNPSize
	Float dbSNPIDXSize
	Float VCFSize
	Float vcfMafFactor
	Float oncoDBSize
	Float oncoDBFactor
	Float fastaSize
	Float ponSize
	Float mafFactor



	command <<<

	set -x
	#generate metrics requires BAM, BAI, DBSNP
	echo -ne "((${bamSize}+${baiSize}+${dbSNPSize}+${dbSNPIDXSize})/1000000000)+5.0" | perl -ne 'print eval($_);'|grep -Po '^\d+' > met_size.dat ;
	#vcf to maf needs to hold VCF and vcf_size and db and unpacked DB
	echo -ne "((${VCFSize}+${VCFSize}*${vcfMafFactor}+${oncoDBSize}+${oncoDBSize}*${oncoDBFactor})/1000000000)+5.0" | perl -ne 'print eval($_);' | grep -Po '^\d+' > vcf_to_maf_size.dat
	#OBF filter neds bam, bai, reference (data, dict, idx), and DBSNP
	echo -ne "((${bamSize}+${baiSize}+${dbSNPSize}+${dbSNPIDXSize}+${fastaSize})/1000000000)+5.0"| perl -ne 'print eval($_);' | grep -Po '^\d+' > obf_disk.dat
	#maf-pon filter needs PON cytoband and mafffile.  Use vcfsize*mafFactor for size of MAF.  cytoband is small so ignore it
	echo -ne "(${ponSize}+${VCFSize}*${mafFactor})/1000000000+10.0" | perl -ne 'print eval($_);' | grep -Po '^\d+' > mp_disk.dat
	#maf merger needs room for MAFs and intermediate data
	echo -ne "(${VCFSize}*${mafFactor}*100)/1000000000+10"| perl -ne 'print eval($_);' | grep -Po '^\d+' > maf_disk.dat
	#oncotation task needs MAF, VCF, and oncodb
	cp -v vcf_to_maf_size.dat onco_disk.dat 


	>>>

	output {
		Int met_size=read_int("met_size.dat")
		Int vcf_to_maf_size=read_int("vcf_to_maf_size.dat")
		Int obf_disk=read_int("obf_disk.dat")
		Int mp_disk=read_int("mp_disk.dat")
		Int maf_disk=read_int("maf_disk.dat")
		Int onco_disk=read_int("onco_disk.dat")
		}



	runtime {
		docker : "broadinstitute/oncotator@sha256:434b232606f5f7df02a3a64be46ce2214308c9cdad577759524b9e68c4d70979"
		disks: "local-disk 1 HDD"
		memory: "0.1 GB"
		}


	}




workflow FilterWorkflow  {


	File PONFile
	File oncoDBTarBall
	File inVCF
	File cytoBandFile
	String PairID
	String TOTNStr
	Int NMIN
	Float thresh
	Float WCUT
	Int CODING_ONLY
	Int MIN_ALT_COUNT
	String TUM_ID
	File REFERENCE
	File DBSNP
	File DBSNPIDX	
	File tumorBam
	File tumorBamIdx
	File oncotationTarBall
	String caseName
	String ctrlName

	call  sizesCalc {
		input:
			bamSize=size(tumorBam),
			baiSize=size(tumorBamIdx),
			dbSNPSize=size(DBSNP),
			dbSNPIDXSize=size(DBSNPIDX),
			VCFSize=size(inVCF),
			oncoDBSize=size(oncotationTarBall),
			fastaSize=size(REFERENCE),
			ponSize=size(PONFile)
		}



	call generateMetrics_task {
		input:
			ID=TUM_ID,
			REFERENCE=REFERENCE,
			BAM=tumorBam,
			BAM_IDX=tumorBamIdx,
			DBSNP=DBSNP,
			DBSNPIDX=DBSNPIDX,
			diskGB=sizesCalc.met_size
		}


	call VCF_to_MAF_task {
		input:
			pairName=PairID,
			inputVCF=inVCF,
			oncoDBTarBall=oncoDBTarBall,
			diskGB=sizesCalc.vcf_to_maf_size,
			caseName=caseName,
			ctrlName=ctrlName
		}

	call MAFPonFilter  {
		input : 
			MAFFile=VCF_to_MAF_task.outputPairMaf,
			PONFile=PONFile,
			cytoBandFile=cytoBandFile,
			PairID=PairID,
			TOTNStr=TOTNStr,
			NMIN=NMIN,
			thresh=thresh,
			WCUT=WCUT,
			CODING_ONLY=CODING_ONLY,
			MIN_ALT_COUNT=MIN_ALT_COUNT,
			diskGB=sizesCalc.mp_disk
		}


	call OrientationBias_filter_Task as oxoGOBF {
		input :
			detailMetrics=generateMetrics_task.metricsFile,
			ID=TUM_ID,
			BAM=tumorBam,
			BAI=tumorBamIdx,
			MAF=VCF_to_MAF_task.outputPairMaf,
			REFERENCE=REFERENCE,
			DBSNP=DBSNP,
			DBSNPIDX=DBSNPIDX,
			diskGB=sizesCalc.obf_disk
		}

	call OrientationBias_filter_Task as ffpeOBF {
		input : 
			detailMetrics=generateMetrics_task.metricsFile,
			ID=TUM_ID,
			BAM=tumorBam,
			BAI=tumorBamIdx,
			MAF=VCF_to_MAF_task.outputPairMaf,
			REFERENCE=REFERENCE,
			DBSNP=DBSNP,
			DBSNPIDX=DBSNPIDX,
			diskGB=sizesCalc.obf_disk
		}

	call maf_merge_task {
		input:
			PairID=PairID,
			in_maf=VCF_to_MAF_task.outputPairMaf,
			oxog_zip=oxoGOBF.MatLabAndFilterDataZip,
			ffpe_zip=ffpeOBF.MatLabAndFilterDataZip,
			pon_zip=MAFPonFilter.mpOutzip,
			in_vcf=inVCF,
			diskGB=sizesCalc.maf_disk
		}

    call pf_counter_task as maf_pon_counter {
    	input:
    		allMaf=MAFPonFilter.allMaf,
    		passMaf=MAFPonFilter.passMaf,
    		diskGB=sizesCalc.maf_disk
    	}


    call pf_counter_task as obf_ffpe_counter {
    	input:
    		allMaf=ffpeOBF.unfiltered_maf,
    		passMaf=ffpeOBF.filtered_maf,
    		diskGB=sizesCalc.maf_disk
    	}


    call pf_counter_task as obf_oxoG_counter {
    	input:
    		allMaf=oxoGOBF.unfiltered_maf,
    		passMaf=oxoGOBF.filtered_maf,
    		diskGB=sizesCalc.maf_disk
    	}


	call oncotation_task {
		input:
			oncotationTarBall=oncotationTarBall,
			inVCF=maf_merge_task.merged_filtered_vcf,
			PairID=PairID,
			diskGB=sizesCalc.onco_disk
		}

	call pf_counter_task as filtered_counter {
		input:
			allMaf=maf_merge_task.union_maf,
			passMaf=maf_merge_task.intersection_maf,
			diskGB=sizesCalc.maf_disk
		}


	output {
		oxoGOBF.q_val
		ffpeOBF.q_val
		maf_merge_task.merged_filtered_vcf
		generateMetrics_task.metricsFile
		oncotation_task.annotatedFilteredVCF
		MAFPonFilter.mpOutzip
		ffpeOBF.MatLabAndFilterDataZip
		oxoGOBF.MatLabAndFilterDataZip
		maf_merge_task.union_maf
		maf_merge_task.intersection_maf
		filtered_counter.all_count
		filtered_counter.pass_count
		oxoGOBF.filtered_maf
		oxoGOBF.unfiltered_maf
		}

	}
