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

		#split intervals
		 /usr/local/jre1.8.0_73/bin/java \
		 	-jar /usr/local/bin/hellbender-all-1.0.jar  CreateEvenIntervals\
			-R ${refFasta} --intervals ${mutectIntervals} --number_of_ways ${nWay} \
			--output_intervals even_intervals.txt

		#turn them into a file of "-L" command-line args to interval_list files
		SPLIT_BASE="split_base";
		NUM_LINES=`wc -l even_intervals.txt|grep -Po '^\s*\d+\s'|grep -Po '\d+'` ;
		for LINE in `seq -w 1 $NUM_LINES`; do 
			HEAD_LINE=`echo $LINE|sed -r "s/^0+//g"` ;
			#USE HEAD_LINE to remove 0-padding/prefix
			echo "Turning line $HEAD_LINE of even_intervals.txt into a file ..." ;
			TGT_FILE="$SPLIT_BASE.$LINE.interval_list" ;
			grep -P '^@' ${mutectIntervals} > $TGT_FILE ; 
			head --lines=$HEAD_LINE even_intervals.txt| tail -1 |tr " " "\n"|grep -Pv '\-L'|grep -P ':'|tr ":" "\t" |tr "-" "\t"|awk '{print $1 "\t" $2 "\t" $3 "\t+\t" NR }' >> $TGT_FILE
			#count number of bases covered and incorporate into size calc for disk; multiply each base mFactor (incorporates call rate and size in bytes used per call)
			BUFFER_FILE=split_base_buffer.$LINE.buffer 
			cat $TGT_FILE|grep -Pv '^@'|awk '{print (($3-$2)+1)*${mFactor} }'|tr "\n" "+"|sed -r "s/\+$//g"|perl -ne 'print eval($_);' > $BUFFER_FILE ; 
			SIZE_FILE=split_base_size.$LINE.size
			cat $BUFFER_FILE | awk '{ print ($1 +${tBamSize}+${tBaiSize}+${nBamSize}+${nBaiSize}+${dbSNPVCFSize}+${rgBLSize}+${cosmicVCFSize}+${normalPanelSize}+${fastaSize}+${fastaDictSize}+${fastaIdxSize}+2000000000)/1000000000 }' | grep -Po '^\d+' > $SIZE_FILE
			done ;

		#generate indices to index sizes and interval files
		#subtract 1 because indexes are zero-based
		find split_base.*|grep -Pn '.'|grep -Po '^\d+' | awk '{print $1-1}' >  indices.dat
		cat split_base_size.* > split_base_sizes_disk.dat 

		#sizes for disks for MutectFC, VEP, oncotator
		NUM_FC_BASES=`cat ${mutectForcecallIntervals}| awk '{print ($3-$2)+1}' |tr "\n" "+"|sed -r "s/\+$//g"|perl -ne 'print eval($_);'` ; 
		echo -ne "$NUM_FC_BASES" | awk '{ print ($1*${mfcFactor}+${tBamSize}+${tBaiSize}+${nBamSize}+${nBaiSize}+${dbSNPVCFSize}+${rgBLSize}+${cosmicVCFSize}+${normalPanelSize}+${fastaSize}+${fastaDictSize}+${fastaIdxSize}+2000000000)/1000000000 }' | grep -Po '^\d+' > fc_size.dat
		NUM_REG_BASES=`cat ${mutectIntervals} | awk '{print (($3-$2)+1)}'|tr "\n" "+"|sed -r "s/\+$//g"|perl -ne 'print eval($_);'` ; 
		echo -ne "$NUM_REG_BASES" |awk '{ print (($1*${mFactor}*10+${oncoDBSize}+${oncoDBSize}*4)+2000000000)/1000000000 }' | grep -Po '^\d+' > onco_disk_size.dat
		echo -ne "$NUM_REG_BASES" |awk '{ print (($1*${mFactor}*10+${vepZipSize}+${vepZipSize}*4)+2000000000)/1000000000 }' | grep -Po '^\d+' > vep_disk_size.dat

		>>>

	output  {
		Array[File] interval_files=glob("split_base.*")
		Array[Int] mutectDiskGBArr=read_lines("split_base_sizes_disk.dat")
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
	File dbSNPVCF
	File cosmicVCF
	Int downsampleToCoverage
	File readGroupBlackList
	File normalPanel

	command <<<
	#increase verbosity
	set -x	

	#variable for normal panel
	NORMAL_PANEL_FLAG_AND_VAL=""
	if [ -s "${normalPanel}" ] ; then
		NORMAL_PANEL_FLAG_AND_VAL="--normal_panel ${normalPanel}" ;
	fi ;

	#mutect 1
	java -jar -Xmx4g /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
	 -L ${mutectIntervals}  --normal_sample_name NORMAL_SAMPLE -I:normal  ${normalBam}  \
	 --tumor_sample_name TUMOR_SAMPLE -I:tumor ${tumorBam}  \
	 --reference_sequence ${refFasta} \
	 --fraction_contamination ${fracContam}  --dbsnp ${dbSNPVCF} \
	 --cosmic ${cosmicVCF}  --read_group_black_list ${readGroupBlackList}  \
	 --out MuTect1.call_stats.txt --coverage_file MuTect1.coverage.wig.txt \
	 --power_file MuTect1.power.wig.txt --downsample_to_coverage ${downsampleToCoverage} \
	 $NORMAL_PANEL_FLAG_AND_VAL

	>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta@sha256:2285809244dabf3aefe11f4e16e91fe3d541b894f8325e4f6c9c98a8dfc81deb"
		memory: "3 GB"
		disks: "local-disk ${diskGB} HDD"
		}


	output {
		File mutect1_cs="MuTect1.call_stats.txt"
		File mutect1_pw="MuTect1.power.wig.txt"
		File mutect1_cw="MuTect1.coverage.wig.txt"
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
	File dbSNPVCF
	File cosmicVCF
	Int downsampleToCoverage
	File readGroupBlackList
	File normalPanel

	command <<<
	#increase verbosity
	set -x

	#variable for normal panel
	NORMAL_PANEL_FLAG_AND_VAL=""
	if [ -s "${normalPanel}" ] ; then
		NORMAL_PANEL_FLAG_AND_VAL="--normal_panel ${normalPanel}" ;
	fi ;	

	#mutect 2
	/usr/local/jre1.8.0_73/bin/java -jar -Xmx4g /usr/local/bin/GenomeAnalysisTK.jar --analysis_type MuTect2 \
	 -L ${mutectIntervals}  -I:normal  ${normalBam}  \
	 -I:tumor ${tumorBam}  \
	 --reference_sequence ${refFasta} \
	 --dbsnp ${dbSNPVCF}  --read_group_black_list ${readGroupBlackList}  \
	 --cosmic ${cosmicVCF} \
	 --out MuTect.call_stats.txt \
	 $NORMAL_PANEL_FLAG_AND_VAL

	>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta@sha256:2285809244dabf3aefe11f4e16e91fe3d541b894f8325e4f6c9c98a8dfc81deb"
		memory: "3 GB"
		disks: "local-disk  ${diskGB}  HDD"
		}

	

	output {
		File mutect2_cs="MuTect.call_stats.txt"
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
	File dbSNPVCF
	File cosmicVCF
	Int downsampleToCoverage
	File readGroupBlackList
	File normalPanel

	command <<<
	#increase verbosity
	set -x

	#variable for normal panel
	NORMAL_PANEL_FLAG_AND_VAL=""
	if [ -s "${normalPanel}" ] ; then
		NORMAL_PANEL_FLAG_AND_VAL="--normal_panel ${normalPanel}" ;
	fi ;		


	java -jar -Xmx4g /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
	 -L ${mutectIntervals}  --normal_sample_name NORMAL_SAMPLE -I:normal  ${normalBam}  \
	 --tumor_sample_name TUMOR_SAMPLE -I:tumor ${tumorBam}  \
	 --reference_sequence ${refFasta} \
	 --fraction_contamination ${fracContam}  --dbsnp ${dbSNPVCF} \
	 --cosmic ${cosmicVCF} \
	 --force_output   --read_group_black_list ${readGroupBlackList}    \
	 --out MuTectFC.call_stats.txt --coverage_file MuTectFC.coverage.wig.txt \
	 --power_file MuTectFC.power.wig.txt --downsample_to_coverage ${downsampleToCoverage} \
	 $NORMAL_PANEL_FLAG_AND_VAL

	>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta@sha256:2285809244dabf3aefe11f4e16e91fe3d541b894f8325e4f6c9c98a8dfc81deb"
		memory: "3 GB"
		disks: "local-disk  ${diskGB}  HDD"
		}

	output {
		File mutectfc_cs="MuTectFC.call_stats.txt"
		File mutectfc_pw="MuTectFC.power.wig.txt"
		File mutectfc_cw="MuTectFC.coverage.wig.txt"
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

	command <<<
		#increase verbosity
		set -x

		#mutect1 call_stats merging
		MUTECT1_CS="MuTect1.call_stats.txt"
		head --lines=2 ${mutect1_cs[0]} > $MUTECT1_CS
		cat ${sep =' ' mutect1_cs} | grep -Pv '#'|grep -Pv '^contig' >> $MUTECT1_CS

		#mutect2 call_stats merging
		MUTECT2_CS="MuTect2.call_stats.txt"
		cat ${mutect2_cs[0]} |grep -P '^#' > $MUTECT2_CS ;
		cat ${sep=' ' mutect2_cs} |grep -Pv '^#' >> $MUTECT2_CS ;

		#convert them to VCFs
		MUTECT1_VCF="MuTect1.call_stats.vcf"
		MUTECT2_VCF="MuTect2.call_stats.vcf"
		python /usr/local/bin/callstats_to_vcf.py "MuTect1.call_stats" $MUTECT1_CS ;
		ln -vs $MUTECT2_CS $MUTECT2_VCF ;
		
		#split MuTect2 outputs into indels and non-indels
		MUTECT2_INDELS="MuTect2.call_stats.indels.vcf"
		MUTECT2_OTHER="MuTect2.call_stats.other.vcf"
		python /usr/local/bin/vcf_partition.py $MUTECT2_VCF $MUTECT2_INDELS $MUTECT2_OTHER

		#merge the M1 and M2 together 
		MUTECT_MERGED_RAW="MuTect.1.2.call_stats.M1All.M2IndelsOnly.raw.vcf"
		MUTECT_MERGED_FILTERED="MuTect.1.2.call_stats.M1All.M2IndelsOnly.filtered.vcf"
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
		File preannotated_vcf="MuTect.1.2.call_stats.M1All.M2IndelsOnly.filtered.vcf"
		File oncotator_log="oncotator.log"
		File oncotator_out="MuTect.1.2.call_stats.M1All.M2IndelsOnly.filtered.vcf.annotated.vcf"
		File mutect1Merged="MuTect1.call_stats.txt"
		File mutect2Merged="MuTect2.call_stats.txt"
		File VEP_VCF="MuTect.1.2.call_stats.M1All.M2IndelsOnly.raw.vcf"
		}


	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta@sha256:2285809244dabf3aefe11f4e16e91fe3d541b894f8325e4f6c9c98a8dfc81deb"
		memory: "${memoryGB} GB"
		disks: "local-disk  ${diskGB}  HDD"
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
	Int mutectDiskGB
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
					tumorBam=tumorBam,
					normalBam=normalBam,
					tumorBamIdx=tumorBamIdx,
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
					diskGB=CallSomaticMutations_131_Prepare_Task.mutectDiskGBArr[idx]
				}

			call Mutect2_Task {
				input:
					tumorBam=tumorBam,
					normalBam=normalBam,
					tumorBamIdx=tumorBamIdx,
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
					diskGB=CallSomaticMutations_131_Prepare_Task.mutectDiskGBArr[idx]
				}
			}

	call MutectFC_Task {
		input:
			tumorBam=tumorBam,
			normalBam=normalBam,
			tumorBamIdx=tumorBamIdx,
			normalBamIdx=normalBamIdx,
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
			vcfFile=GatherAndOncotate_Task.oncotator_out,
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
	       GatherAndOncotate_Task.preannotated_vcf
	       GatherAndOncotate_Task.oncotator_log
	       GatherAndOncotate_Task.mutect1Merged
	       GatherAndOncotate_Task.mutect2Merged
	       VEP_Task.VEP_Output
	       VEP_Task.VEP_Report
	       mutect_nozzle_task.nozzleReport
	       mutect_nozzle_task.nozzleRData
	       }

	}
