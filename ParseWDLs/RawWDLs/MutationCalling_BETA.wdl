task CallSomaticMutations_131_Prepare {
	File tumorBam
	File tumorBamIdx
	File normalBam
	File normalBamIdx
	File refFastaIdx
	File mutectIntervals 

	command <<<
		#increase verbosity
		set -x

		mkdir -pv tumor_split
		mkdir -pv normal_split
		mkdir -pv interval_split
		export SPLIT_INDEX=0 ;
		for SEQUENCE in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT ;
		#for SEQUENCE in `cut -f1 ${refFastaIdx}` ; 
		#for SEQUENCE in `cut -f1 ${refFastaIdx}|grep -Pi '^\s*2[02]\s*$'` ; 
			do 
			CHROM=$SEQUENCE ;
			echo "Now splitting with CHROM=$CHROM ..."

			#split normal bam by chrom
			samtools view ${normalBam} $CHROM -b > normal_split/$CHROM.bam ;
			#split tumor bam by chrom
			samtools view ${tumorBam} $CHROM -b > tumor_split/$CHROM.bam ;
			#split intervals
			INT_SPLIT_EXT=`echo ${mutectIntervals}|grep -Po '\.[^\.]*$'|tr -d "."` ;
			echo "Acquired extension $INT_SPLIT_EXT "
			INT_SPLIT_FILE="interval_split/$CHROM.$INT_SPLIT_EXT" ; 
			echo "To write split intervals to $INT_SPLIT_FILE"
			cat ${mutectIntervals} | grep -P '^@' > $INT_SPLIT_FILE
			cat ${mutectIntervals} | awk -v CHROM_COMP=$CHROM '{if($1==CHROM_COMP) print $0}' >> $INT_SPLIT_FILE

			#see if it's a valid interval file
			NUM_INTS=`grep -Pi "^$SEQUENCE" $INT_SPLIT_FILE |cut -f2,3|grep -Pc '\d+\s+\d+'` ;
			if [ "$NUM_INTS" -gt "0" ] ; then 
				# a valid interval file
				echo "At least one mutect interval detected on sequence $SEQUENCE " ;
				#save integers to file for subsequent iteration for indexing
				echo "$SPLIT_INDEX" >> split_indices.dat ;
				SPLIT_INDEX=`echo $SPLIT_INDEX + 1 |bc` ;
			else 
				# no intervals detected!
				echo "No mutect intervals detected on sequence $SEQUENCE .  So deleting data split on the sequence . " ;
				rm -v normal_split/$CHROM.bam ;
				rm -v tumor_split/$CHROM.bam ; 
				rm -v $INT_SPLIT_FILE ;
			fi ;

		done ;

		touch split_indices.dat ; 

		>>>

	output  {
		Array[File] tumor_bams=glob("tumor_split/*.bam")
		Array[File] normal_bams=glob("normal_split/*.bam")
		Array[File] intervals=glob("interval_split/*")
		Array[Int] split_indices=read_lines("split_indices.dat")
		}

	runtime {
		docker: "broadinstitute/broadmutationcalling_beta"
		memory: "24 GB"
		disks: "local-disk 100 HDD"
		}

	}

task Mutect1Task {

	File tumorBam
	File normalBam
	File mutectIntervals
	File refFasta
	File refFastaIdx
	File refFastaDict
	String fracContam
	File dbSNPVCF
	File cosmicVCF
	String downsampleToCoverage
	File readGroupBlackList
	File normalPanel

	command <<<
	#increase verbosity
	set -x	

	#index the bams first
	samtools index ${tumorBam}
	samtools index ${normalBam}

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
	 --cosmic ${cosmicVCF} \
	 --out MuTect1.call_stats.txt --coverage_file MuTect1.coverage.wig.txt \
	 --power_file MuTect1.power.wig.txt --downsample_to_coverage ${downsampleToCoverage} \
	 $NORMAL_PANEL_FLAG_AND_VAL

	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_beta"
		memory: "24 GB"
		disks: "local-disk 100 HDD"
		}


	output {
		File mutect1_cs="MuTect1.call_stats.txt"
		File mutect1_pw="MuTect1.power.wig.txt"
		File mutect1_cw="MuTect1.coverage.wig.txt"
		}

	}

task Mutect2Task {

	File tumorBam
	File normalBam
	File mutectIntervals
	File refFasta
	File refFastaIdx
	File refFastaDict
	String fracContam
	File dbSNPVCF
	File cosmicVCF
	String downsampleToCoverage
	File readGroupBlackList
	File normalPanel

	command <<<
	#increase verbosity
	set -x

	#index the bams first
	samtools index ${tumorBam}
	samtools index ${normalBam}

	#variable for normal panel
	NORMAL_PANEL_FLAG_AND_VAL=""
	if [ -s "${normalPanel}" ] ; then
		NORMAL_PANEL_FLAG_AND_VAL="--normal_panel ${normalPanel}" ;
	fi ;	

	#mutect 2
	java -jar -Xmx4g /usr/local/bin/GenomeAnalysisTK.jar --analysis_type MuTect2 \
	 -L ${mutectIntervals}  -I:normal  ${normalBam}  \
	 -I:tumor ${tumorBam}  \
	 --reference_sequence ${refFasta} \
	 --dbsnp ${dbSNPVCF} \
	 --cosmic ${cosmicVCF} \
	 --out MuTect.call_stats.txt \
	 $NORMAL_PANEL_FLAG_AND_VAL

	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_beta"
		memory: "24 GB"
		disks: "local-disk 100 HDD"
		}

	

	output {
		File mutect2_cs="MuTect.call_stats.txt"
		}

	}

task MutectFCTask {

	File tumorBam
	File normalBam
	File mutectIntervals
	File refFasta
	File refFastaIdx
	File refFastaDict
	String fracContam
	File dbSNPVCF
	File cosmicVCF
	String downsampleToCoverage
	File readGroupBlackList
	File normalPanel

	command <<<
	#increase verbosity
	set -x

	#index the bams first
	samtools index ${tumorBam}
	samtools index ${normalBam}

	#variable for normal panel
	NORMAL_PANEL_FLAG_AND_VAL=""
	if [ -s "${normalPanel}" ] ; then
		NORMAL_PANEL_FLAG_AND_VAL="--normal_panel ${normalPanel}" ;
	fi ;		

	#mutect force-calling from CallSomaticMutationsForceCalling_45
	java -Xmx4g -jar /usr/local/bin/muTect-qscore.jar --read_group_black_list ${readGroupBlackList} \
	 -rf BadCigar --analysis_type MuTect -L ${mutectIntervals} \
	 --normal_sample_name NORMAL_SAMPLE -I:normal ${normalBam} \
	 --tumor_sample_name TUMOR_SAMPLE -I:tumor ${tumorBam} --reference_sequence ${refFasta}\
	 --dbsnp ${dbSNPVCF}  --cosmic ${cosmicVCF} --out MuTectFC.call_stats.txt \
	 --coverage_file MuTectFC.coverage.wig.txt --power_file MuTectFC.power.wig.txt \
	 --enable_extended_output --enable_qscore_output --downsample_to_coverage ${downsampleToCoverage}\
	 $NORMAL_PANEL_FLAG_AND_VAL

	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_beta"
		memory: "24 GB"
		disks: "local-disk 100 HDD"
		}

	output {
		File mutectfc_cs="MuTectFC.call_stats.txt"
		File mutectfc_pw="MuTectFC.power.wig.txt"
		File mutectfc_cw="MuTectFC.coverage.wig.txt"
		}

	}

task GatherAndOncotateAndVEP {

	Array[File] mutect1_cs
	Array[File] mutect1_pw
	Array[File] mutect1_cw
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

		#obtain the name of the directory for oncodb
		ONCO_DB_DIR_NAME=`gunzip -c ${oncoDBTarBall} |tar -tf /dev/stdin|head -1` ; 

		#unpack the oncodir! (may take a while...)
		tar -xzf ${oncoDBTarBall}

		#Run the merged filtered VCF (from both mutects through Oncotator) 
		/usr/local/lib/python2.7/site-packages/Oncotator-1.8.0.0-py2.7.egg/oncotator/Oncotator.py\
			 -i VCF --db-dir `pwd`/$ONCO_DB_DIR_NAME -o VCF $MUTECT_MERGED_FILTERED $MUTECT_MERGED_FILTERED.annotated.vcf hg19

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
		--cache  --offline  -i $MUTECT_MERGED_RAW --everything

		>>>


	output {
		Array[File] all_outs=glob("*")
		File oncotator_log="oncotator.log"
		File oncotator_out="MuTect.1.2.call_stats.M1All.M2IndelsOnly.filtered.vcf.annotated.vcf"
		File mutect1Merged="MuTect1.call_stats.txt"
		File mutect2Merged="MuTect2.call_stats.txt"
		File VEP_Output="variant_effect_output.txt"
		File VEP_Report="variant_effect_output.txt_summary.html"
		}


	runtime {
		docker: "broadinstitute/broadmutationcalling_beta"
		memory: "104 GB"
		disks: "local-disk 100 HDD"
		
		}

	}

workflow CallingGroupWorkflow {

	File tumorBam
	File tumorBamIdx
	File normalBam
	File normalBamIdx
	File refFastaIdx
	File mutectIntervals
	File mutectForcecallLintervals
	File refFasta
	File refFastaDict
	String fracContam
	File dbSNPVCF
	File cosmicVCF
	String downsampleToCoverage
	File readGroupBlackList
	File? normalPanel
	File oncoDBTarBall
	File VEP_File

	# PREPARE FOR SCATTER
	call CallSomaticMutations_131_Prepare {
		input: 
			tumorBam=tumorBam,
			tumorBamIdx=tumorBamIdx,
			normalBam=normalBam,
			normalBamIdx=normalBamIdx,
			refFastaIdx=refFastaIdx,
			mutectIntervals=mutectIntervals
		}

	#SCATTER AND ANALYZE
	scatter (i in CallSomaticMutations_131_Prepare.split_indices) {
			
			call Mutect1Task {
				input:
					tumorBam=CallSomaticMutations_131_Prepare.tumor_bams[i],
					normalBam=CallSomaticMutations_131_Prepare.normal_bams[i],
					mutectIntervals=CallSomaticMutations_131_Prepare.intervals[i],
					refFasta=refFasta,
					refFastaIdx=refFastaIdx,
					refFastaDict=refFastaDict,
					fracContam=fracContam,
					dbSNPVCF=dbSNPVCF,
					cosmicVCF=cosmicVCF,
					downsampleToCoverage=downsampleToCoverage,
					readGroupBlackList=readGroupBlackList,
					normalPanel=normalPanel
				}

			call Mutect2Task {
				input:
					tumorBam=CallSomaticMutations_131_Prepare.tumor_bams[i],
					normalBam=CallSomaticMutations_131_Prepare.normal_bams[i],
					mutectIntervals=CallSomaticMutations_131_Prepare.intervals[i],
					refFasta=refFasta,
					refFastaIdx=refFastaIdx,
					refFastaDict=refFastaDict,
					fracContam=fracContam,
					dbSNPVCF=dbSNPVCF,
					cosmicVCF=cosmicVCF,
					downsampleToCoverage=downsampleToCoverage,
					readGroupBlackList=readGroupBlackList,
					normalPanel=normalPanel
				}
			}

	call MutectFCTask {
		input:
			tumorBam=tumorBam,
			normalBam=normalBam,
			mutectIntervals=mutectForcecallLintervals,
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			fracContam=fracContam,
			dbSNPVCF=dbSNPVCF,
			cosmicVCF=cosmicVCF,
			downsampleToCoverage=downsampleToCoverage,
			readGroupBlackList=readGroupBlackList,
			normalPanel=normalPanel
			}

	call GatherAndOncotateAndVEP {
		input:
			mutect1_cs=Mutect1Task.mutect1_cs,
			mutect1_pw=Mutect1Task.mutect1_pw,
			mutect1_cw=Mutect1Task.mutect1_cw,
			mutect2_cs=Mutect2Task.mutect2_cs,
			oncoDBTarBall=oncoDBTarBall,
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			VEP_File=VEP_File
		}


	}