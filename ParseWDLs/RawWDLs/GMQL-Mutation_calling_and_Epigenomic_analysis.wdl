task CallSomaticMutations_131_Prepare_Task {
    #RUNTIME INPUT PARAMS
	Int preemptible

	#TASK INPUT PARAMS
	File mutectIntervals 
	File refFastaIdx
	File refFasta
	File refFastaDict
	File mutectForcecallIntervals
	Int nWay
	Float tBamSize
	Float tBaiSize
	Float nBamSize
	Float nBaiSize
	Float dbSNPVCFSize
	Float rgBLSize
	Float cosmicVCFSize
	Float normalPanelSize
	Float fcIntervalsSize
	Float fastaSize
	Float fastaDictSize
	Float fastaIdxSize
	Float mFactor
	Float mfcFactor
	Float oncoDBSize
	Float vepZipSize


	command <<<
		#increase verbosity
		set -x

		echo "TMPDIR IS $TMPDIR"
		mkdir -pv $TMPDIR

		#Calculate disk size for all shards of the mutect run
		SIZE_FILE=split_base_sizes_disk.dat
		awk 'BEGIN { print int((${tBamSize}+${tBaiSize}+${nBamSize}+${nBaiSize}+${dbSNPVCFSize}+${rgBLSize}+${cosmicVCFSize}+${normalPanelSize}+${fastaSize}+${fastaDictSize}+${fastaIdxSize}+2000000000)/1000000000+1) }' > $SIZE_FILE
		
		#Create list of indices for the scatter job
		seq 0 $((${nWay}-1)) > indices.dat

		#Run the prepare task that splits the .interval_list file into subfiles
		java -Xmx2g -jar /usr/local/bin/GatkScatterGatherPrepare.jar . ${nWay} --intervals ${mutectIntervals} --reference_sequence ${refFasta} 

		#sizes for disks for MutectFC, VEP, oncotator
		NUM_FC_BASES=`cat ${mutectForcecallIntervals}| awk '{print ($3-$2)+1}' |tr "\n" "+"|sed -r "s/\+$//g"|perl -ne 'print eval($_);'` ; 
		echo -ne "$NUM_FC_BASES" | awk '{ print ($1*${mfcFactor}+${tBamSize}+${tBaiSize}+${nBamSize}+${nBaiSize}+${dbSNPVCFSize}+${rgBLSize}+${cosmicVCFSize}+${normalPanelSize}+${fastaSize}+${fastaDictSize}+${fastaIdxSize}+2000000000)/1000000000 }' | grep -Po '^\d+' > fc_size.dat
		NUM_REG_BASES=`cat ${mutectIntervals} | awk '{print (($3-$2)+1)}'|tr "\n" "+"|sed -r "s/\+$//g"|perl -ne 'print eval($_);'` ; 
		echo -ne "$NUM_REG_BASES" |awk '{ print (($1*${mFactor}*10+${oncoDBSize}+${oncoDBSize}*4)+2000000000)/1000000000 }' | grep -Po '^\d+' > onco_disk_size.dat
		echo -ne "$NUM_REG_BASES" |awk '{ print (($1*${mFactor}*10+${vepZipSize}+${vepZipSize}*4)+2000000000)/1000000000 }' | grep -Po '^\d+' > vep_disk_size.dat

		>>>

	output  {
		Array[File] interval_files=glob("gatk-scatter.*")
		Int mutectDiskGB=read_int("split_base_sizes_disk.dat")
		Array[Int] scatterIndices=read_lines("indices.dat")
		Int mutectFCDisk=read_int("fc_size.dat")
		Int oncotatorDisk=read_int("onco_disk_size.dat")
		Int vepDisk=read_int("vep_disk_size.dat")
		}

	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta@sha256:2285809244dabf3aefe11f4e16e91fe3d541b894f8325e4f6c9c98a8dfc81deb"
		memory: "1 GB"
		}
	}

task Mutect1_Task {
        #RUNTIME INPUT PARAMS
	Int preemptible
	Int diskGB

	#TASK INPUT PARAMS
	File mutectIntervals 
	File tumorBam
	File normalBam
	File tumorBamIdx
	File normalBamIdx
	File refFasta
	File refFastaIdx
	File refFastaDict
	Float fracContam
	Float contamFloor
	File dbSNPVCF
	File cosmicVCF
	Int downsampleToCoverage
	File readGroupBlackList
	File normalPanel
	String pairName
	String caseName
	String ctrlName


	command <<<
	#increase verbosity
	set -x	

	#variable for normal panel
	NORMAL_PANEL_FLAG_AND_VAL=""
	if [ -s "${normalPanel}" ] ; then
		NORMAL_PANEL_FLAG_AND_VAL="--normal_panel ${normalPanel}" ;
	fi ;

	#compute/apply contamination floor for effective contamination
	EFFECTIVE_CONTAMINATION=`/usr/local/bin/python -c 'import sys;print sys.argv[1] if(float(sys.argv[1])>=float(sys.argv[2])) else sys.argv[2]'  ${fracContam} ${contamFloor}` ;
	echo "EFFECTIVE_CONTAMINATION computed to be $EFFECTIVE_CONTAMINATION"

	#mutect 1
	java -jar -Xmx4g /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
	 -L ${mutectIntervals}  --normal_sample_name ${ctrlName} -I:normal  ${normalBam}  \
	 --tumor_sample_name ${caseName} -I:tumor ${tumorBam}  \
	 --reference_sequence ${refFasta} \
	 --fraction_contamination $EFFECTIVE_CONTAMINATION --dbsnp ${dbSNPVCF} \
	 --cosmic ${cosmicVCF}  --read_group_black_list ${readGroupBlackList}  \
	 --out ${pairName}.MuTect1.call_stats.txt --coverage_file ${pairName}.MuTect1.coverage.wig.txt \
	 --power_file ${pairName}.MuTect1.power.wig.txt  --downsample_to_coverage ${downsampleToCoverage} \
	 $NORMAL_PANEL_FLAG_AND_VAL

	>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta@sha256:2285809244dabf3aefe11f4e16e91fe3d541b894f8325e4f6c9c98a8dfc81deb"
		memory: "3 GB"
		disks: "local-disk ${diskGB} HDD"
		}


	output {
		File mutect1_cs="${pairName}.MuTect1.call_stats.txt"
		File mutect1_pw="${pairName}.MuTect1.power.wig.txt"
		File mutect1_cw="${pairName}.MuTect1.coverage.wig.txt"
		}

	}

task Mutect2_Task {
        #RUNTIME INPUT PARAMS
	Int preemptible
	Int diskGB

	#TASK INPUT PARAMS
	File tumorBam
	File normalBam
	File tumorBamIdx
	File normalBamIdx
	File mutectIntervals
	File refFasta
	File refFastaIdx
	File refFastaDict
	Float fracContam
	Float contamFloor
	File dbSNPVCF
	Int downsampleToCoverage
	File readGroupBlackList
	File normalPanel
	String pairName
	String caseName
	String ctrlName
	File GATK4_JAR


	command <<<
	#increase verbosity
	set -x

	#variable for normal panel
	NORMAL_PANEL_FLAG_AND_VAL=""
	if [ -s "${normalPanel}" ] ; then
		BZ="${normalPanel}.gz"
		#bgzip the file and index it
		bgzip ${normalPanel}
		tabix $BZ 
		NORMAL_PANEL_FLAG_AND_VAL="--normal_panel $BZ" ;
	fi ;	

	#compute/apply contamination floor for effective contamination
	EFFECTIVE_CONTAMINATION=`/usr/local/bin/python -c 'import sys;print sys.argv[1] if(float(sys.argv[1])>=float(sys.argv[2])) else sys.argv[2]'  ${fracContam} ${contamFloor}` ;
	echo "EFFECTIVE_CONTAMINATION computed to be $EFFECTIVE_CONTAMINATION"

	#DBSNP indexing
	BZD="${dbSNPVCF}.gz"
	bgzip ${dbSNPVCF}
	tabix $BZD


	#Mutect2 wants names that match those in the BAMs so grab them from the BAMs
	BAM_NORMAL_NAME=`samtools view -H ${normalBam} |grep -Po 'SM:.*'|cut -f1|grep -Po ':.*'|perl -ne '$_=~s/^://g ; print $_'|head -1`
	echo "To use bam normal name $BAM_NORMAL_NAME" ;
	BAM_TUMOR_NAME=`samtools view -H ${tumorBam} |grep -Po 'SM:.*'|cut -f1|grep -Po ':.*'|perl -ne '$_=~s/^://g ; print $_'|head -1`
	echo "To use bam tumor name $BAM_TUMOR_NAME" ; 
	#mutect 2 ----- gatk4
	/usr/local/jre1.8.0_73/bin/java -jar -Xmx4g ${GATK4_JAR} Mutect2 \
	-R ${refFasta} \
	-I ${tumorBam} \
	-tumor "$BAM_TUMOR_NAME" \
	-I ${normalBam} \
	-normal "$BAM_NORMAL_NAME" \
	--dbsnp $BZD \
	--contamination_fraction_to_filter $EFFECTIVE_CONTAMINATION \
	$NORMAL_PANEL_FLAG_AND_VAL \
	-L ${mutectIntervals} \
	-O ${pairName}.MuTect2.call_stats.unfiltered.unaf.txt

	#filter the variants
	/usr/local/jre1.8.0_73/bin/java -jar -Xmx4g ${GATK4_JAR} FilterMutectCalls \
	 -O ${pairName}.MuTect2.call_stats.filtered.unaf.txt -V ${pairName}.MuTect2.call_stats.unfiltered.unaf.txt

process_af_py_str="
import re
from argparse import ArgumentParser

def processVCF(IN_VCF,OUT_VCF,TUMOR_STR,NORMAL_STR):
	reader=open(IN_VCF,'r')
	writer=open(OUT_VCF,'w')
	tumor_piece_idx=None
	for line in reader:
		sline=line.strip()
		if(sline.startswith('##')):
			#if begins with '#' then pass-thru
			writer.write(sline+'\n')
		elif(sline.startswith('#CHROM')):
			pieces=sline.split('\t')
			tmr_idx=pieces.index(TUMOR_STR)
			nrml_idx=pieces.index(NORMAL_STR)
			writer.write(sline+'\n')
		elif(sline.startswith('#')):
			writer.write(sline+'\n')
		else:
			pieces=sline.split('\t')
			info_piece=pieces[8]
			info_pieces=info_piece.split(':')
			AD_idx=info_pieces.index('AD')
			AF_idx=info_pieces.index('AF')
			tmr_piece=pieces[tmr_idx]
			nrml_piece=pieces[nrml_idx]
			##FORMAT=<ID=AD,Number=R,Type=Integer,Description='Allelic depths for the ref and alt alleles in the order listed>
			tmr_pieces=tmr_piece.split(':')
			nrml_pieces=nrml_piece.split(':')
			tmr_ads=tmr_pieces[AD_idx]
			nrml_ads=nrml_pieces[AD_idx]
			tmr_ads_vals=tmr_ads.split(',')
			nrml_ads_vals=nrml_ads.split(',')
			tmr_af=float(tmr_ads_vals[1])/(float(tmr_ads_vals[1])+float(tmr_ads_vals[0]))
			nrml_af=float(nrml_ads_vals[1])/(float(nrml_ads_vals[1])+float(nrml_ads_vals[0]))
			tmr_pieces[AF_idx]='{:.4f}'.format(tmr_af)
			nrml_pieces[AF_idx]='{:.4f}'.format(nrml_af)
			pieces[tmr_idx]=':'.join(tmr_pieces)
			pieces[nrml_idx]=':'.join(nrml_pieces)
			writer.write('\t'.join(pieces)+'\n')
	writer.close()
	reader.close()

if __name__ == '__main__':
	parser = ArgumentParser(description='computer AF from AD, only AF for indicated TUMOR_STR sample')
	parser.add_argument('IN_VCF',help='input VCF')
	parser.add_argument('OUT_VCF',help='output MAF')
	parser.add_argument('TUMOR_STR',help='name of tumor sample (for acquiring proper sample-column)')
	parser.add_argument('NORMAL_STR',help='name of the normal sample (for acquiring proper sample-column')
	args = parser.parse_args()
	if(args):
		print 'Picked up    IN_VCF  : '+str(args.IN_VCF)
		print 'Picked up   OUT_VCF  : '+str(args.OUT_VCF)
		print 'Picked up TUMOR_STR  : '+str(args.TUMOR_STR)
		print 'Picked up NORMAL_STR : '+str(args.NORMAL_STR)
		processVCF(args.IN_VCF,args.OUT_VCF,args.TUMOR_STR,args.NORMAL_STR)
	else:
		parser.print_help()
"
	/usr/local/bin/python -c "$process_af_py_str"  \
	 ${pairName}.MuTect2.call_stats.filtered.unaf.txt  \
	 ${pairName}.MuTect2.call_stats.txt  \
	 "$BAM_TUMOR_NAME" "$BAM_NORMAL_NAME"


	>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta@sha256:2285809244dabf3aefe11f4e16e91fe3d541b894f8325e4f6c9c98a8dfc81deb"
		memory: "3 GB"
		disks: "local-disk  ${diskGB}  HDD"
		}

	

	output {
		File mutect2_cs="${pairName}.MuTect2.call_stats.txt"
		}

	}

task MutectFC_Task {
        #RUNTIME INPUT PARAMS
	Int preemptible
	Int diskGB

	#TASK INPUT PARAMS
	File tumorBam
	File normalBam
	File tumorBamIdx
	File normalBamIdx
	File mutectIntervals
	File refFasta
	File refFastaIdx
	File refFastaDict
	Float fracContam
	Float contamFloor
	File dbSNPVCF
	File cosmicVCF
	Int downsampleToCoverage
	File readGroupBlackList
	File normalPanel
	String pairName
	String caseName
	String ctrlName

	command <<<
	#increase verbosity
	set -x

	#variable for normal panel
	NORMAL_PANEL_FLAG_AND_VAL=""
	if [ -s "${normalPanel}" ] ; then
		NORMAL_PANEL_FLAG_AND_VAL="--normal_panel ${normalPanel}" ;
	fi ;		

	#compute/apply contamination floor for effective contamination
	EFFECTIVE_CONTAMINATION=`/usr/local/bin/python -c 'import sys;print sys.argv[1] if(float(sys.argv[1])>=float(sys.argv[2])) else sys.argv[2]'  ${fracContam} ${contamFloor}` ;
	echo "EFFECTIVE_CONTAMINATION computed to be $EFFECTIVE_CONTAMINATION"




	java -jar -Xmx4g /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
	 -L ${mutectIntervals}  --normal_sample_name ${ctrlName} -I:normal  ${normalBam}  \
	 --tumor_sample_name ${caseName} -I:tumor ${tumorBam}  \
	 --reference_sequence ${refFasta} \
	 --fraction_contamination $EFFECTIVE_CONTAMINATION  --dbsnp ${dbSNPVCF} \
	 --cosmic ${cosmicVCF} \
	 --force_output   --read_group_black_list ${readGroupBlackList}    \
	 --out ${pairName}.MuTect1.call_stats.txt --coverage_file ${pairName}.MuTect1.coverage.wig.txt \
	 --power_file ${pairName}.MuTect1.power.wig.txt  --downsample_to_coverage ${downsampleToCoverage} \
	 $NORMAL_PANEL_FLAG_AND_VAL

	>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta@sha256:2285809244dabf3aefe11f4e16e91fe3d541b894f8325e4f6c9c98a8dfc81deb"
		memory: "3 GB"
		disks: "local-disk  ${diskGB}  HDD"
		}

	output {
		File mutectfc_cs="${pairName}.MuTect1.call_stats.txt"
		File mutectfc_pw="${pairName}.MuTect1.power.wig.txt"
		File mutectfc_cw="${pairName}.MuTect1.coverage.wig.txt"
		}

	}

task GatherAndOncotate_Task {
        #RUNTIME INPUT PARAMS
	Int preemptible
	Int memoryGB
	Int diskGB

	#TASK INPUT PARAMS
	Array[File] mutect1_cs
	Array[File] mutect2_cs
	File refFastaDict
	File oncoDBTarBall
	File VEP_File
	File refFasta
	File refFastaIdx
	String pairName

	command <<<
		#increase verbosity
		set -x

		#mutect1 call_stats merging
		MUTECT1_CS="${pairName}.MuTect1.call_stats.txt"
		head --lines=2 ${mutect1_cs[0]} > $MUTECT1_CS
		cat ${sep =' ' mutect1_cs} | grep -Pv '#'|grep -Pv '^contig' >> $MUTECT1_CS

		#mutect2 call_stats merging
		MUTECT2_CS="${pairName}.MuTect2.call_stats.txt"
		cat ${mutect2_cs[0]} |grep -P '^#' > $MUTECT2_CS ;
		cat ${sep=' ' mutect2_cs} |grep -Pv '^#' >> $MUTECT2_CS ;

		#convert them to VCFs
		MUTECT1_VCF="${pairName}.MuTect1.call_stats.vcf"
		MUTECT2_VCF="${pairName}.MuTect2.call_stats.vcf"
		python /usr/local/bin/callstats_to_vcf.py "${pairName}.MuTect1.call_stats" $MUTECT1_CS ;
		ln -vs $MUTECT2_CS $MUTECT2_VCF ;
		
		#split MuTect2 outputs into indels and non-indels
		MUTECT2_INDELS="${pairName}.MuTect2.call_stats.indels.vcf"
		MUTECT2_OTHER="${pairName}.MuTect2.call_stats.other.vcf"
		python /usr/local/bin/vcf_partition.py $MUTECT2_VCF $MUTECT2_INDELS $MUTECT2_OTHER

		#merge the M1 and M2 together 
		MUTECT_MERGED_RAW="${pairName}.MuTect.1.2.call_stats.M1All.M2IndelsOnly.raw.vcf"
		MUTECT_MERGED_FILTERED="${pairName}.MuTect.1.2.call_stats.M1All.M2IndelsOnly.filtered.vcf"
		#merge and output raw data
		java -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK.jar  -T CombineVariants -R ${refFasta} \
		--variant $MUTECT1_VCF --variant $MUTECT2_INDELS -o $MUTECT_MERGED_RAW \
		--assumeIdenticalSamples  -U ALLOW_SEQ_DICT_INCOMPATIBILITY
		#merge and output filtered data
		#Step 1 is get headers
		cat $MUTECT_MERGED_RAW|grep -P '^#' > $MUTECT_MERGED_FILTERED
		#Step 2 is get pass-only lines
		cat $MUTECT_MERGED_RAW|grep -Pv '^#' | awk '{if($7~/PASS/) print $0}' >> $MUTECT_MERGED_FILTERED

		#find TARBALL type
		TYPE=`echo 'if("${oncoDBTarBall}"=~m/z$/) { print "GZ" ; } else { print "TAR" ; } '| perl` ; 

		#obtain the name of the directory for oncodb
		#and unpack based on the TYPE
		if [ "$TYPE" == "GZ" ] ; then
			ONCO_DB_DIR_NAME=`gunzip -c ${oncoDBTarBall} |tar -tf /dev/stdin|head -1` ; 
			tar -xzf ${oncoDBTarBall}
		else
			ONCO_DB_DIR_NAME=`tar -tf ${oncoDBTarBall} |head -1` ;
			tar -xf ${oncoDBTarBall} ;
		fi ;
		echo "ONCO_DB_DIR_NAME is $ONCO_DB_DIR_NAME" ; 

		#Run the merged filtered VCF (from both mutects through Oncotator) 
		/usr/local/lib/python2.7/site-packages/Oncotator-1.8.0.0-py2.7.egg/oncotator/Oncotator.py\
			 -i VCF --db-dir `pwd`/$ONCO_DB_DIR_NAME -o VCF $MUTECT_MERGED_FILTERED $MUTECT_MERGED_FILTERED.annotated.vcf hg19

		>>>


	output {
		#Array[File] all_outs=glob("*")
		File WXS_Mutation_M1M2_noindels_filetered_not_oncotated="${pairName}.MuTect.1.2.call_stats.M1All.M2IndelsOnly.filtered.vcf"
		File oncotator_log="oncotator.log"
		File WXS_Mutation_M1M2_noindels_filetered_oncotated="${pairName}.MuTect.1.2.call_stats.M1All.M2IndelsOnly.filtered.vcf.annotated.vcf"
		File WXS_Mutation_M1="${pairName}.MuTect1.call_stats.txt"
		File WXS_Mutation_M2="${pairName}.MuTect2.call_stats.txt"
		File VEP_VCF="${pairName}.MuTect.1.2.call_stats.M1All.M2IndelsOnly.raw.vcf"
		}


	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta@sha256:2285809244dabf3aefe11f4e16e91fe3d541b894f8325e4f6c9c98a8dfc81deb"
		memory: "${memoryGB} GB"
		disks: "local-disk  ${diskGB}  HDD"
		}
	}
    
    
task gmql {
  File inFile				# the file given by the OncotatorTask containine the called mutations
  String repoPath			# Path to the ChiP-Seq experiment
  String annotationPath		# Path to the dataset of genome annotations

  command <<<
  set -x
  /usr/src/myapp/vcf2gdm.sh ${inFile} ds
  
  query="refSeqGenes = SELECT(annotation_type == 'gene' AND provider == 'RefSeq') %s; \
         myExp = SELECT() %s ; \
         myConfirmedExp = COVER(2, ANY) myExp;\
         myExpOverGenes  = JOIN(DIST < 0; output: RIGHT_DISTINCT) refSeqGenes myConfirmedExp;
         myMut  = SELECT(parser: bedscoreparser)  /usr/src/myapp/ds/; \
         myMutOverExp = MAP() myExpOverGenes myMut;\
         myFilteredExp = SELECT(region: count_myExpOverGenes_myMut > 0) myMutOverExp;\
         MATERIALIZE myFilteredExp INTO file:///usr/src/myapp/results/;"

  printf "$query" ${annotationPath} ${repoPath} > /usr/src/myapp/query.gmql
  
  cat  /usr/src/myapp/query.gmql

  java -jar /usr/src/myapp/uber-GMQL-Cli-2.0.jar -scriptpath /usr/src/myapp/query.gmql
  
  tar -zcvf out.tar.gz /usr/src/myapp/results/exp/
  
  >>>
  output {
    File response = stdout()
    File out = "out.tar.gz"
  }
  
  runtime {
		docker: "pp86/gmql"
        memory: "24 GB"
    	disks: "local-disk 100 SSD"
    }
}
    
    
    

task VEP_Task {

	#input data
	File mutectMergedRawVCF
	File VEP_File

	#runtime parameters
	Int preemptible
	Int memoryGB
	Int diskGB	

	command <<<

		#increase verbosity
		set -x

		#make a link from the home directory to the current directory to avoid running out of disk space on the boot disk
		mkdir -v vep_data_dir
		#delete the existing directory first to make a successful link
		rm -rf ~/.vep
		ln -vs `pwd`/vep_data_dir ~/.vep


		#Run the merged RAW VCF from the call stats through VEP
		#In either case unpack the data into the home directory where VEP
		#expects to find it
		IS_ZIP=`echo ${VEP_File}|grep -Pic '\.zip$'` ;
		if [ "$IS_ZIP" -eq "1" ] ;
		then
			#it's a zip file
			unzip -d ~ ${VEP_File} ;
		else
			#tar ball
			tar -C ~ -xvzf ${VEP_File}
		fi ;
		/ensembl-tools-release-83/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
		--cache  --vcf  --offline  -i ${mutectMergedRawVCF} --symbol

		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta@sha256:2285809244dabf3aefe11f4e16e91fe3d541b894f8325e4f6c9c98a8dfc81deb"
		memory: "${memoryGB} GB"
		disks: "local-disk  ${diskGB}  HDD"
		}

	output {
		File VEP_Output="variant_effect_output.txt"
		File VEP_Report="variant_effect_output.txt_summary.html"
		}


	}





task mutect_nozzle_task {
        #RUNTIME INPUT PARAMS
        Int diskGB

	#TASK INPUT PARAMS
	File vcfFile
	File vepReportHTML
	File vepVCFFile
	String pairName
	File fcCallStats

	command <<<

		#increase verbosity
		set -x

		#make data local in directory
		mv -vf ${vcfFile} variants.vcf
		mv -vf ${vepVCFFile} variant_effect_output.txt
		mv -vf ${vepReportHTML} vepReport.html
		mv -vf ${fcCallStats} force_called.call_stats.txt

		#count the number of variants
		NUM_VARIANTS=`grep -Pv '^#' variants.vcf |wc -l|grep -Po '\d+'` ;

		#extract variant data for display within HTML
		if [ "$NUM_VARIANTS" -eq "0" ] ; then
			#there are no variants
			touch variant_nozzle_table.dat ; 
		else
			#extract variant data to display
			HEADER_LINE_NUM=`grep -Pn '^#CHROM' variants.vcf  |tr ":" "\t"|cut -f1` ; 
			cat variants.vcf|awk -v HEADER=$HEADER_LINE_NUM '{if(NR>=HEADER) print $0}'|sed -r "s/#CHROM/CHROM/g" | cut -f 1,2,3,4,5|sed -r "s/[ACGT]{10,}/\&lt;TRUNCATED_FOR_LENGTH\&gt;/g"  > variant_nozzle_table.dat
		fi ;

		#make the report
		/R/RunR.sh  /R/makeMutectNozzleReport.R ${pairName} \
		  $NUM_VARIANTS variant_nozzle_table.dat variants.vcf \
		  vepReport.html vepVCF.vcf \
		  force_called.call_stats.txt

		#assure files exist for delocalization
		touch nozzle.html nozzle.RData variants.vcf variant_nozzle_table.dat vepReport.html force_called.call_stats.txt

		>>>
		

	runtime {
		docker: "broadinstitute/mutect_nozzle@sha256:868e36e7058ae26c1d2b4888ccc2fa74be00817ceab1cbb1bfd5c157a8892c39"
		disks: "local-disk  ${diskGB}  HDD"
		memory: "3 GB"
		}


	output {

		#nozzle files
		File nozzleReport="nozzle.html"
		File nozzleRData="nozzle.RData"

		#linked files
		File variantsVCF="variants.vcf"
		File variantTable="variant_nozzle_table.dat"
		File vepReportLocalHTML="vepReport.html"
		File vepVCFOut="variant_effect_output.txt"
		File fcCallStatsLocal="force_called.call_stats.txt"

		}

	}

workflow CallingGroup_Workflow {
	#RUNTIME INPUT PARAMS
	Int preemptible
	Int oncoMemoryGB
	Int VEPMemoryGB

	#WORKFLOW INPUT PARMS
	File tumorBam
	File tumorBamIdx
	File normalBam
	File normalBamIdx
	File refFastaIdx
	File mutectIntervals
	File mutectForcecallIntervals
	File refFasta
	File refFastaDict
	Float fracContam
	File dbSNPVCF
	File cosmicVCF
	Int downsampleToCoverage
	File readGroupBlackList
	File normalPanel
	File oncoDBTarBall
	File VEP_File
	Int nWay
	String pairName
	String caseName
	String ctrlName
  
    # GMQL datasets paths
    String repoPath
    String annotationPath

	# PREPARE FOR SCATTER
	call CallSomaticMutations_131_Prepare_Task {
		input: 
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			mutectIntervals=mutectIntervals,
			nWay=nWay,
			preemptible=preemptible,
			tBamSize=size(tumorBam),
			tBaiSize=size(tumorBamIdx),
			nBamSize=size(normalBam),
			nBaiSize=size(normalBamIdx),
			dbSNPVCFSize=size(dbSNPVCF),
			rgBLSize=size(readGroupBlackList),
			cosmicVCFSize=size(cosmicVCF),
			normalPanelSize=size(normalPanel),
			fcIntervalsSize=size(mutectForcecallIntervals),
			fastaSize=size(refFasta),
			fastaDictSize=size(refFastaDict),
			fastaIdxSize=size(refFastaIdx),
			mutectForcecallIntervals=mutectForcecallIntervals,
			oncoDBSize=size(oncoDBTarBall),
			vepZipSize=size(VEP_File)
		}

	#SCATTER AND ANALYZE
	scatter (idx in CallSomaticMutations_131_Prepare_Task.scatterIndices) {
			
			call Mutect1_Task {
				input:
					contamFloor=0.01,
					tumorBam=tumorBam,
					normalBam=normalBam,
					tumorBamIdx=tumorBamIdx,
					pairName=pairName,
					caseName=caseName,
					ctrlName=ctrlName,
					normalBamIdx=normalBamIdx,
					mutectIntervals=CallSomaticMutations_131_Prepare_Task.interval_files[idx],
					refFasta=refFasta,
					refFastaIdx=refFastaIdx,
					refFastaDict=refFastaDict,
					fracContam=fracContam,
					dbSNPVCF=dbSNPVCF,
					cosmicVCF=cosmicVCF,
					downsampleToCoverage=downsampleToCoverage,
					readGroupBlackList=readGroupBlackList,
					normalPanel=normalPanel,
					preemptible=preemptible,
					diskGB=CallSomaticMutations_131_Prepare_Task.mutectDiskGB
				}

			call Mutect2_Task {
				input:
					contamFloor=0.01,
					tumorBam=tumorBam,
					normalBam=normalBam,
					tumorBamIdx=tumorBamIdx,
					normalBamIdx=normalBamIdx,
					mutectIntervals=CallSomaticMutations_131_Prepare_Task.interval_files[idx],
					refFasta=refFasta,
					pairName=pairName,
					refFastaIdx=refFastaIdx,
					refFastaDict=refFastaDict,
					fracContam=fracContam,
					dbSNPVCF=dbSNPVCF,
					downsampleToCoverage=downsampleToCoverage,
					readGroupBlackList=readGroupBlackList,
					normalPanel=normalPanel,
					preemptible=preemptible,
					pairName=pairName,
					caseName=caseName,
					ctrlName=ctrlName,
					diskGB=CallSomaticMutations_131_Prepare_Task.mutectDiskGB
				}
			}

	call MutectFC_Task {
		input:
			contamFloor=0.01,
			tumorBam=tumorBam,
			normalBam=normalBam,
			tumorBamIdx=tumorBamIdx,
			normalBamIdx=normalBamIdx,
			pairName=pairName,
			caseName=caseName,
			ctrlName=ctrlName,
			mutectIntervals=mutectForcecallIntervals,
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			fracContam=fracContam,
			dbSNPVCF=dbSNPVCF,
			cosmicVCF=cosmicVCF,
			downsampleToCoverage=downsampleToCoverage,
			readGroupBlackList=readGroupBlackList,
			normalPanel=normalPanel,
			preemptible=preemptible,
			diskGB=CallSomaticMutations_131_Prepare_Task.mutectFCDisk
			}


	call GatherAndOncotate_Task {
		input:
			pairName=pairName,
			mutect1_cs=Mutect1_Task.mutect1_cs,
			mutect2_cs=Mutect2_Task.mutect2_cs,
			oncoDBTarBall=oncoDBTarBall,
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			VEP_File=VEP_File,
			preemptible=preemptible,
			memoryGB=oncoMemoryGB,
			diskGB=CallSomaticMutations_131_Prepare_Task.oncotatorDisk
		}
        
  # We call the gmql task after the execution of the mutation calling tasks
    call gmql{
       input:
          inFile=GatherAndOncotate_Task.WXS_Mutation_M1M2_noindels_filetered_oncotated, 		# the mutation data will be the one provided by the GatherAndOncotate_Task
          repoPath=repoPath,
          annotationPath=annotationPath
    }

	call VEP_Task {
		input:
			mutectMergedRawVCF=GatherAndOncotate_Task.VEP_VCF,
			VEP_File=VEP_File,
			preemptible=preemptible,
			memoryGB=VEPMemoryGB,
			diskGB=CallSomaticMutations_131_Prepare_Task.vepDisk
		}

	call mutect_nozzle_task {
		input:
			vcfFile=GatherAndOncotate_Task.WXS_Mutation_M1M2_noindels_filetered_oncotated,
			vepReportHTML=VEP_Task.VEP_Report,
			vepVCFFile=VEP_Task.VEP_Output,
			pairName=pairName,
			fcCallStats=MutectFC_Task.mutectfc_cs,
			diskGB=CallSomaticMutations_131_Prepare_Task.mutectFCDisk
		}

	
	output {
	       MutectFC_Task.mutectfc_cs
	       MutectFC_Task.mutectfc_pw
	       MutectFC_Task.mutectfc_cw
	       GatherAndOncotate_Task.WXS_Mutation_M1M2_noindels_filetered_not_oncotated
	       GatherAndOncotate_Task.WXS_Mutation_M1M2_noindels_filetered_oncotated
	       GatherAndOncotate_Task.oncotator_log
	       GatherAndOncotate_Task.WXS_Mutation_M1
	       GatherAndOncotate_Task.WXS_Mutation_M2
	       VEP_Task.VEP_Output
	       VEP_Task.VEP_Report
	       mutect_nozzle_task.nozzleReport
	       mutect_nozzle_task.nozzleRData
           gmql.out
	       }

	}
