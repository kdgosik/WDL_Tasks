task QC_Prepare_Task {
	#RUNTIME INPUT PARAMS
	Int preemptible

	#sizes of files
	Float tBamBytes
	Float tBaiBytes
	Float nBamBytes
	Float nBaiBytes
	Float regionFileBytes
	Float rgBLBytes
	Float capNormDBZipBytes
	Float fastaBytes
	Float fastaDictBytes
	Float fastaIdxBytes
	Float exomeIntervalsBytes
	Float snpSixBytes
	Float hapMapVCFBytes
	Float hapDBForCCBytes
	Float dbSNPVCFBytes
	Float dbSNPVCFIDXBytes
	Float picardHapMapVCFBytes
	Float picardTargetIntervalsBytes
	Float picardBaitIntervalsBytes



	command <<<

	#increase verbosity
	set -x

	#disk sizes calculations
	#qc scatter disk (add two gigabytes for round-up and buffer)
	echo "${tBamBytes}+${tBaiBytes}+${nBamBytes}+${nBaiBytes}+${regionFileBytes}+2000000000" | perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > qc_scatter_disk.dat
	#qc gather disk
	echo "${tBamBytes}+${tBaiBytes}+${nBamBytes}+${nBaiBytes}+${rgBLBytes}+${capNormDBZipBytes}+2000000000"| perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > qc_gather_disk.dat
	#contest disk
	echo "${tBamBytes}+${tBaiBytes}+${nBamBytes}+${nBaiBytes}+${fastaBytes}+${fastaDictBytes}+${fastaIdxBytes}+${exomeIntervalsBytes}+${snpSixBytes}+${hapMapVCFBytes}+2000000000" | perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > contest_disk.dat
	#cross-check lane fingerprints
	echo "${tBamBytes}+${tBaiBytes}+${nBamBytes}+${nBaiBytes}+${hapDBForCCBytes}+2000000000"| perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > cc_lane_fp.dat
	#tumorMM
	echo "${tBamBytes}*2.0+${tBaiBytes}*2.0+${dbSNPVCFBytes}+${dbSNPVCFIDXBytes}+${fastaBytes}+3000000000" | perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > tumor_mm.dat
	#tumor BAM validation
	echo "${tBamBytes}+${tBaiBytes}+10000000000" | perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > tumor_val.dat
	#tumor RGFP
	echo "${tBamBytes}+${tBaiBytes}+${picardHapMapVCFBytes}+2000000000"  | perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > tumor_rgfp.dat
	#tumor HS metrics
	echo "${tBamBytes}+${tBaiBytes}+${picardTargetIntervalsBytes}+${picardBaitIntervalsBytes}+10000000000"| perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > tumor_hsm.dat
	#normalMM
	echo "${nBamBytes}*2.0+${nBaiBytes}*2.0+${dbSNPVCFBytes}+${dbSNPVCFIDXBytes}+${fastaBytes}+3000000000" | perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > nrml_mm.dat
	#normal BAM validation
	echo "${nBamBytes}+${nBaiBytes}+10000000000" | perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > nrml_val.dat
	#normal RGFP
	echo "${nBamBytes}+${nBaiBytes}+${picardHapMapVCFBytes}+2000000000"  | perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > nrml_rgfp.dat
	#normal HS metrics
	echo "${nBamBytes}+${nBaiBytes}+${picardTargetIntervalsBytes}+${picardBaitIntervalsBytes}+10000000000"| perl -ne 'print (eval($_)/1000000000)."\n";'|grep -Po '^\d+' > nrml_hsm.dat


	#split by chromosome
	for CHROM in `seq 1 24` ;
		do
			echo $CHROM >> split_chromes.dat ;

		done ;
	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		memory: "0.1 GB" #specifying memory this low and no CPU will trigger usage of f1-micro
		preemptible: "${preemptible}"
		}

	output {
		Array[Int] split_chroms=read_lines("split_chromes.dat")
		Int qc_scatter_disk=read_int("qc_scatter_disk.dat")
		Int qc_gather_disk=read_int("qc_gather_disk.dat")
		Int contest_disk=read_int("contest_disk.dat")
		Int cc_lane_fp_disk=read_int("cc_lane_fp.dat")
		Int tumor_mm_disk=read_int("tumor_mm.dat")
		Int tumor_validation_disk=read_int("tumor_val.dat")
		Int tumor_rgfp_disk=read_int("tumor_rgfp.dat")
		Int tumor_hsm_disk=read_int("tumor_hsm.dat")
		Int nrml_mm_disk=read_int("nrml_mm.dat")
		Int nrml_validation_disk=read_int("nrml_val.dat")
		Int nrml_rgfp_disk=read_int("nrml_rgfp.dat")
		Int nrml_hsm_disk=read_int("nrml_hsm.dat")
		}
}

task QC_Scatter_Task {
        #RUNTIME INPUT PARAMS
	Int preemptible
	Int diskGB

	#TASK INPUT PARAMS
	File tumorBam
	File tumorBamIdx
	File normalBam
	File normalBamIdx
	File regionFile
	Int CHROM

	command <<<

		#increase verbosity
		set -x

		#scatter by chromsome
		OUTPATH=`echo -ne "tumor.chr${CHROM}.control.rcl"` ;
		java -jar /usr/local/bin/RegionCovPerLane.jar ${tumorBam} ${regionFile} $OUTPATH ${CHROM}
		OUTPATH=`echo -ne "normal.chr${CHROM}.control.rcl"` ;
		java -jar /usr/local/bin/RegionCovPerLane.jar ${normalBam} ${regionFile} $OUTPATH ${CHROM}


		>>>
	
	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "${preemptible}"
		}

	output {
		File tumorRCL="tumor.chr${CHROM}.control.rcl"
		File normalRCL="normal.chr${CHROM}.control.rcl"
		}

	}

task QC_Gather_Task {
        #RUNTIME INPUT PARAMS
	Int preemptible
	Int diskGB

	#TASK INPUT PARAMS
	File tumorBam
	File tumorBamIdx
	File normalBam
	File normalBamIdx
	File regionFile
	File readGroupBlackList
	File captureNormalsDBRCLZip
	Array[File] tumorRCLs
	Array[File] normalRCLs


	command  <<<
	#increase verbosity
	set -x

	###############################################
	# Copy Number QC for Capture
	#Make lane lists for tumor and normal
	#MakeLaneList_12
	java -jar /usr/local/bin/MakeLaneList.jar ${tumorBam}  case.lanelist.txt ;
	java -jar /usr/local/bin/MakeLaneList.jar ${normalBam} control.lanelist.txt ; 

	#Gather region coverages per lane
	cat ${sep=' ' tumorRCLs} > gathered.tumor.rcl
	cat ${sep=' ' normalRCLs} > gathered.normal.rcl

	#run the mat-lab based report
	#CopyNumberQCReport_27
	cp -vfr /CopyNumberQCReport_27/unzip/* . 
	#Make a file of files paths
	#This command aims to list the zip files contents, filtering out all but the file paths and then write the paths to a file
	unzip -l ${captureNormalsDBRCLZip} |sed -r "s/^\s+//g"|sed -r "s/^[0-9]+//g"|sed -r "s/^\s*//g"|sed -r "s/^\S*//g"|sed -r "s/^\s*[0-9:]*\s*//g"|grep -Pv '^Date\s+Time\s+Name'|grep -Pv '^[\-\s]*$'|grep -Pv '\.zip$'|grep -Pv '^files$' > capture_normals_db_wdl
	#this command actually performs the unzipping
	unzip ${captureNormalsDBRCLZip}
	./run_fh_CopyNumberQCReport.sh /opt/MATLAB/MATLAB_Compiler_Runtime/v710/  PairCopyNumQCReport  \
	   gathered.tumor.rcl  gathered.normal.rcl  case.lanelist.txt  control.lanelist.txt \
	     ${readGroupBlackList}  ${regionFile} capture_normals_db_wdl  NA NA
	#zip up the output
	zip CopyNumQCout.zip report.html num.mixups.txt PairCopyNumQCReport* *.png

	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "${preemptible}"
		}

	output {
		#lane lists
		File tumorBamLaneList="case.lanelist.txt"
		File normalBamLaneList="control.lanelist.txt"
		File tumorRCL="gathered.tumor.rcl"
		File normalRCL="gathered.normal.rcl"
		#CopyNumQC
		File CopyNumQCOutZip="CopyNumQCout.zip"
		File CopyNumQCReport="report.html"
		File CopyNumQCReportPNG="PairCopyNumQCReport_CopyNumberQC.png"
		File CopyNumQCMixUps="num.mixups.txt"
		Int CopyNumMixupsInt=read_int("num.mixups.txt")
		}

	}


task PicardMultipleMetrics_Task {
     	#RUNTIME INPUT PARAMS
	Int preemptible
	Int diskGB

	#TASK INPUT PARAMS
	File bam
	File bamIndex
	File refFasta
	File DB_SNP_VCF
	File DB_SNP_VCF_IDX

	command <<<
	#increase verbosity
	set -x


	BASE=`basename ${bam}` ; 
	echo "Processing ${bam} ..." ;

	#run bam through CleanSam to set MAPQ of unmapped reads to zero
	/usr/local/jre1.8.0_73/bin/java -Xmx4g -jar /usr/local/bin/picardtools-2.1.0.jar  CleanSam \
		I=${bam} O=$BASE.unmapped_reads_cleaned.bam

	/usr/local/jre1.8.0_73/bin/java   -Xmx4g  -jar  /usr/local/bin/picardtools-2.1.0.jar   CollectMultipleMetrics \
		CREATE_MD5_FILE=true   I=$BASE.unmapped_reads_cleaned.bam O=$BASE.multiple_metrics R=${refFasta} PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution  PROGRAM=MeanQualityByCycle \
		PROGRAM=CollectBaseDistributionByCycle PROGRAM=CollectSequencingArtifactMetrics  \
		PROGRAM=CollectQualityYieldMetrics PROGRAM=CollectGcBiasMetrics ;

	#oxoG metrics
	/usr/local/jre1.8.0_73/bin/java -jar /usr/local/bin/picard.1.895.jar \
		ConvertSequencingArtifactToOxoG I=$BASE.multiple_metrics    R=${refFasta}

	#zip up reports 
	zip picard_multiple_metrics.zip $BASE.multiple_metrics.*


	#acquire BAM name
	basename ${bam}|sed -r "s/\.[^\.]+$//g" > bam_name_file.txt


	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		memory: "7 GB"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "${preemptible}"
		}


	output {
		String bamName=read_string("bam_name_file.txt")
		File metricsReportsZip="picard_multiple_metrics.zip"
		}

	}


task BAMValidation_Task {
        #RUNTIME INPUT PARAMS
	Int diskGB
	Int preemptible

	#TASK INPUT PARAMS
	File bam
	File bamIndex

	command <<<
		#increase verbosity
		set -x

		#validation verbose
		/usr/local/jre1.8.0_73/bin/java  -Xmx4g  -jar /usr/local/bin/picardtools-2.1.0.jar ValidateSamFile \
			I=${bam} MODE=VERBOSE > validation_verbose.txt 2>&1

		#look for errors found
		NO_ERRORS_FOUND_COUNT=`grep -Pc '^No\s+errors\s+found\s*$'` ;
		echo $NO_ERRORS_FOUND_COUNT > no_errors_found_count.txt ;
	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		memory: "3 GB"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "${preemptible}"
		}

	output {
		File validation_verbose="validation_verbose.txt"
		File no_errors_found_count_file="no_errors_found_count.txt"
		Int no_errors_found_count=read_int("no_errors_found_count.txt")
		}

	}

task crossCheckRGFP_Task {
        #RUNTIME INPUT PARAMS
	Int diskGB
	Int preemptible
	
	#TASK INPUT PARAMS
	File bam
	File bamIndex
	File picardHapMap

	command <<<

	#increase verbosity
	set -x

	#drop from haplotypeDB seq entries which aren't in BAM if there are any found
	PREPPED_HAPLOTYPE_DB=PreppedHaplotypeDB.txt
	/usr/local/bin/filter_not_in_bam_dict.pl ${bam} ${picardHapMap} $PREPPED_HAPLOTYPE_DB	

	#cross check
	/usr/local/jre1.8.0_73/bin/java  -Xmx4g  -jar /usr/local/bin/picardtools-2.1.0.jar \
	CrosscheckReadGroupFingerprints HAPLOTYPE_MAP=$PREPPED_HAPLOTYPE_DB \
	 I=${bam} O="picard_crosscheck_report.txt" EXIT_CODE_WHEN_MISMATCH=0

	#obtain int value for anything not EXPECTED MATCH
	UNEXPECTED_MATCH_COUNT=`cut -f1 picard_crosscheck_report.txt |grep -Pv 'RESULT'|grep -Pv 'EXPECTED.MATCH'|wc -l|grep -Po '\d+'` ;
	echo $UNEXPECTED_MATCH_COUNT > unexpected_match_count.txt

	>>>

	output {
		File picardCCREport="picard_crosscheck_report.txt"
		File unexpectedMatchFile="unexpected_match_count.txt"
		Int unexpectedMatch=read_int("unexpected_match_count.txt")
		}

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		memory: "3 GB"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "${preemptible}"
		}
	}

task HSMetrics_Task {
        #RUNTIME INPUT PARAMS
	Int preemptible
	Int diskGB

	#TASK INPUT PARAMS
	File bam
	File bamIndex
	File picardTargetIntervals
	File picardBaitIntervals

	command <<<

	#increase verbosity
	set -x

	#drop from targets and bait items in their seq list that don't appear in the bam list
	PREPPED_TGTS=TargetIntervalsFiltered.txt
	/usr/local/bin/filter_not_in_bam_dict.pl ${bam} ${picardTargetIntervals} $PREPPED_TGTS
	PREPPED_BAITS=BaitIntervalsFiltered.txt
	/usr/local/bin/filter_not_in_bam_dict.pl ${bam} ${picardBaitIntervals} $PREPPED_BAITS

	#HS metrics
	/usr/local/jre1.8.0_73/bin/java  -Xmx4g  -jar /usr/local/bin/picardtools-2.1.0.jar CollectHsMetrics I=${bam} \
	 BAIT_INTERVALS=$PREPPED_BAITS TARGET_INTERVALS=$PREPPED_TGTS OUTPUT="HSMetrics.txt"

	#get 20x coverage value
	cat HSMetrics.txt | grep -A 10000 -P '^##.METR' | grep -B 1000 -m 1 -P '^\s*$' | grep -Pv '##' | grep -Pv '^\s*$' > temp_table.txt
	TWENTY_X_COL=`head -1 temp_table.txt|tr "\t" "\n"|grep -Pin 'PCT_TARGET_BASES_20X'|grep -Po '^\d+'`;
	cut -f $TWENTY_X_COL temp_table.txt|tail -1 > PCT_TARGET_BASES_20X.txt

	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		memory: "3 GB"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "${preemptible}"
		}


	output {
		File hsMetrics="HSMetrics.txt"
		File pctTgtTwentyXFile="PCT_TARGET_BASES_20X.txt"
		String pctTgtTwentyX=read_string("PCT_TARGET_BASES_20X.txt")
		}

	}





task ContEST_Task {
        #RUNTIME INPUT PARAMS
	Int preemptible
	Int diskGB

	#TASK INPUT PARAMS
	File tumorBam
	File tumorBamIdx
	File normalBam
	File normalBamIdx
	File refFasta
	File refFastaIdx
	File refFastaDict
	File exomeIntervals
	File SNP6Bed
	File HapMapVCF
	String pairName

	command <<<
	#increase verbosity
	set -x

	#run contest 
	java  -Xmx2048m  -Djava.io.tmpdir=/tmp -jar /usr/local/bin/GenomeAnalysisTK.jar \
	-T ContEst  -I:eval ${tumorBam} \
	-I:genotype ${normalBam} -L ${exomeIntervals} -L ${SNP6Bed}   -isr INTERSECTION  \
	-R ${refFasta} -l INFO -pf ${HapMapVCF} -o contamination.af.txt  --trim_fraction 0.03  \
	--beta_threshold 0.05  -br contamination.base_report.txt -mbc 100  --min_genotype_depth 30  --min_genotype_ratio 0.8  

	#Get contamination column and value
	CONTAM_COL=`head -1 contamination.af.txt|tr "\t" "\n"|grep -Pin '^contamination$'|grep -Po '^\d+:'|tr -d ":"` ;
	CONTAM_VAL=`cut -f $CONTAM_COL contamination.af.txt|tail -1` ;
	CONTAM_NM=`cut -f1 contamination.af.txt|tail -1 `;

	#write contam data for validation
	VAL_CONTAM_FLAG=`echo -ne "$CONTAM_VAL"|grep -Pc '\d*\.\d+'`
	if [ "$VAL_CONTAM_FLAG" -eq "1" ] ;
	then
		# multiply by 1%
		CONTAM_VAL=`echo $CONTAM_VAL*0.01|bc -l` ; 
	else
		# use default of 0.02 if no data found
		CONTAM_VAL="0.02" ; 
	fi ;
	echo -ne "${pairName}\t$CONTAM_VAL\n" > contamination_validation.array_free.txt
	echo $CONTAM_VAL > fraction_contamination.txt

	#Contamination validation/consensus
	python /usr/local/populateConsensusContamination_v26/contaminationConsensus.py --pass_snp_qc false \
	--output contest_validation.output.tsv \
	--annotation contamination_percentage_consensus_capture \
	--array contamination_validation.array_free.txt --noarray contamination_validation.array_free.txt


	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		memory: "3 GB"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "${preemptible}"
		}


	output {
		File contamDataFile="contamination_validation.array_free.txt"
		File contestAFFile="contamination.af.txt"
		File contestBaseReport="contamination.base_report.txt"
		File validationOutput="contest_validation.output.tsv"
		String fraction_contamination=read_string("fraction_contamination.txt")
		}

	}


task CrossCheckLaneFingerprints_Task {
     	#RUNTIME INPUT PARAMS
	Int preemptible
	Int diskGB

	#TASK INPUT PARAMS
	File tumorBam
	File normalBam
	File tumorBamIdx
	File normalBamIdx
	String pairName
	File HaplotypeDBForCrossCheck
	String validationStringencyLevel



	command <<<
		#increase verbosity
		set -x


		#ConvertDependentsList[version=8] (custom!)
		#prepare input for the convert command
		echo -ne "${pairName}_Normal\t${normalBam}\n${pairName}_Tumor\t${tumorBam}\n" > conv_file.input.tsv 
		#make the options file
		OPTIONS_FILE=${pairName}_CCLFP.options
		echo -ne "I=${normalBam}\nI=${tumorBam}\n" > $OPTIONS_FILE ; 
		
		#drop from haplotypeDB seq entries which aren't in BAM if there are any found
		PREPPED_HAPLOTYPE_DB=PreppedHaplotypeDB.txt
		/usr/local/bin/filter_not_in_bam_dict.pl ${normalBam} ${HaplotypeDBForCrossCheck} $PREPPED_HAPLOTYPE_DB


		#Cross Check!
		#CrosscheckLaneFingerprints[version=9]
		mkdir -v tmp
		java -Xmx3g -jar /usr/local/bin/CrosscheckReadGroupFingerprints.jar \
		 OPTIONS_FILE=$OPTIONS_FILE  H=$PREPPED_HAPLOTYPE_DB TMP_DIR=`pwd`/tmp QUIET=true \
		   EXIT_CODE_WHEN_MISMATCH=0 OUTPUT=crosscheck.stats.txt  \
		    VALIDATION_STRINGENCY=${validationStringencyLevel}

		#CrosscheckLaneFingerprintsReport[version=10]
		#Run report
		/usr/local/bin/crosscheck_report.pl --stats crosscheck.stats.txt ${normalBam} ${tumorBam} conv_file.input.tsv 

		#obtain Lowest LOD
		grep -Po '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*$' crosscheck.stats.txt |sort -g|head -1 > lowest.lod.txt


	>>>


	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		memory: "3 GB"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "${preemptible}"
		}


	output {
		File crossCheckStats="crosscheck.stats.txt"
		File crossCheckReport="report.html"
		File lowestLODFile="lowest.lod.txt"
		String lowestLODStr=read_string("lowest.lod.txt")
		}
	}


task QC_Nozzle_Report_Task {

	String pairName
	Int preemptible
	String fracContam
	File CCNormalRGFPReport
	File CCTumorRGFPReport
	File tumorBAMValidation
	File normalBAMValidation
	File picardTumorMultipleMetricsZip
	File picardNormalMultipleMetricsZip
	File tumorHSMetrics
	File normalHSMetrics
	File copyNumQCZip
	String tumorBamName
	String normalBamName
	File CCLaneFingerPrintsReportHTML
	File CCLaneFingerPrintsStats

	command <<<
		#increase verbosity
		set -x

		#acquire contEST value
		CONTEST_VALUE="${fracContam}" ;

		#Keep files linked by nozzle report in current directory
		cp -vf ${CCLaneFingerPrintsReportHTML} local_cc_fp_report.html
		cp -vf ${CCLaneFingerPrintsStats} local_cc_fp_stats.txt

		#copy in HS files
		cp -vf ${tumorHSMetrics}  ${tumorBamName}.tumorHSMetrics.txt
		cp -vf ${normalHSMetrics} ${normalBamName}.normlHSMetrics.txt

		#copy in lane RG FP files
		cp -vf ${CCTumorRGFPReport}  ${tumorBamName}.tumor.CCRGFPReport.txt
		cp -vf ${CCNormalRGFPReport} ${normalBamName}.normal.CCRGFPReport.txt

		#build the report
		/R/RunR.sh  /R/makeQCNozzleReport.R "${pairName}" \
			"${tumorBamName}" "${normalBamName}" \
			 $CONTEST_VALUE ${tumorBamName}.tumor.CCRGFPReport.txt ${normalBamName}.normal.CCRGFPReport.txt \
			 ${tumorBAMValidation} ${normalBAMValidation} \
			 ${picardTumorMultipleMetricsZip} ${picardNormalMultipleMetricsZip} \
			 ${tumorBamName}.tumorHSMetrics.txt ${normalBamName}.normlHSMetrics.txt \
			 ${copyNumQCZip} \
			 local_cc_fp_report.html local_cc_fp_stats.txt

		#touch the report
		touch nozzle.html
		touch nozzle.RData

		#get copy num QC PNG and report
		touch PairCopyNumQCReport_CopyNumberQC.png report.html

		#capture PDFs from picard MM
		#issue "touch" command to help capture PDFs!
		for PDF_EXT in base_distribution_by_cycle.pdf gc_bias.pdf insert_size_histogram.pdf quality_by_cycle.pdf quality_distribution.pdf ;
			do
			PDF_N_PATH=${normalBamName}.bam.multiple_metrics.$PDF_EXT ;
			PDF_T_PATH=${tumorBamName}.bam.multiple_metrics.$PDF_EXT ;
			touch $PDF_N_PATH $PDF_T_PATH ;
		done ;

		>>>
	runtime {
		memory: "3 GB" #specifying memory this low and no CPU will trigger usage of f1-micro
		preemptible: "${preemptible}"
		docker : "broadinstitute/qc_nozzle@sha256:75953685680cec2b958bcf3caa333b793fcc17a8976d579bff51ad7600417340"
		disks: "local-disk 2 HDD"
		}

	output {
		#main nozzle data
		File nozzleHTML="nozzle.html"
		File nozzleRData="nozzle.RData"

		#copy num QC
		File copyNumQCPNG="PairCopyNumQCReport_CopyNumberQC.png"
		File copyNumQCHTML="report.html"

		#picard MM pdf reports (tumor)
		File T_B_DIST="${tumorBamName}.bam.multiple_metrics.base_distribution_by_cycle.pdf"
		File T_GC_BIAS="${tumorBamName}.bam.multiple_metrics.gc_bias.pdf"
		File T_INS_SIZE_HIST="${tumorBamName}.bam.multiple_metrics.insert_size_histogram.pdf"
		File T_QUAL_BY_CYCLE="${tumorBamName}.bam.multiple_metrics.quality_by_cycle.pdf"
		File T_QUAL_DIST="${tumorBamName}.bam.multiple_metrics.quality_distribution.pdf"
		#picard MM pdf report (normal)
		File N_B_DIST="${normalBamName}.bam.multiple_metrics.base_distribution_by_cycle.pdf"
		File N_GC_BIAS="${normalBamName}.bam.multiple_metrics.gc_bias.pdf"
		File N_INS_SIZE_HIST="${normalBamName}.bam.multiple_metrics.insert_size_histogram.pdf"
		File N_QUAL_BY_CYCLE="${normalBamName}.bam.multiple_metrics.quality_by_cycle.pdf"
		File N_QUAL_DIST="${normalBamName}.bam.multiple_metrics.quality_distribution.pdf"

		#picard MM png reports(tumor)
		File T_B_ISM_PNG="${tumorBamName}.bam.multiple_metrics.insert_size_metrics.hist.txt.png"
		File T_QUAL_BY_CYCLE_PNG="${tumorBamName}.bam.multiple_metrics.quality_by_cycle_metrics.hist.txt.png"
		File T_QUAL_DIST_PNG="${tumorBamName}.bam.multiple_metrics.quality_distribution_metrics.hist.txt.png"
		#picard MM png reports(normal)
		File N_B_ISM_PNG="${normalBamName}.bam.multiple_metrics.insert_size_metrics.hist.txt.png"
		File N_QUAL_BY_CYCLE_PNG="${normalBamName}.bam.multiple_metrics.quality_by_cycle_metrics.hist.txt.png"
		File N_QUAL_DIST_PNG="${normalBamName}.bam.multiple_metrics.quality_distribution_metrics.hist.txt.png"

	 	#CC lane FP
	 	File ccLaneFPReportOut="local_cc_fp_report.html"
	 	File ccLaneFPStatsOut="local_cc_fp_stats.txt"

	 	#TABULAR DATA (tumor)
		File asmTableTum="${tumorBamName}.bam.multiple_metrics.alignment_summary_metrics.table.txt"
		File bbdTum="${tumorBamName}.bam.multiple_metrics.bait_bias_detail_metrics.table.txt"
		File bbsTum="${tumorBamName}.bam.multiple_metrics.bait_bias_summary_metrics.table.txt"
		File bdbcTum="${tumorBamName}.bam.multiple_metrics.base_distribution_by_cycle_metrics.table.txt"
		File gcdmTum="${tumorBamName}.bam.multiple_metrics.gc_bias.detail_metrics.table.txt"
		File gcSumTum="${tumorBamName}.bam.multiple_metrics.gc_bias.summary_metrics.table.txt"
		File ismTum="${tumorBamName}.bam.multiple_metrics.insert_size_metrics.table.txt"
		File oxoGTum="${tumorBamName}.bam.multiple_metrics.oxog_metrics.table.txt"
		File padmTum="${tumorBamName}.bam.multiple_metrics.pre_adapter_detail_metrics.table.txt"
		File pasmTum="${tumorBamName}.bam.multiple_metrics.pre_adapter_summary_metrics.table.txt"
		File qymTum="${tumorBamName}.bam.multiple_metrics.quality_yield_metrics.table.txt"
		File tumRGFPReport="${tumorBamName}.tumor.CCRGFPReport.txt"
		File tumHSMetricsHist="${tumorBamName}.tumorHSMetrics.txt.hist.txt"
		File tumHSMetricsTbl="${tumorBamName}.tumorHSMetrics.txt.table.txt"

	 	#TABULAR DATA (normal)
		File asmTableNrm="${normalBamName}.bam.multiple_metrics.alignment_summary_metrics.table.txt"
		File bbdNrm="${normalBamName}.bam.multiple_metrics.bait_bias_detail_metrics.table.txt"
		File bbsNrm="${normalBamName}.bam.multiple_metrics.bait_bias_summary_metrics.table.txt"
		File bdbcNrm="${normalBamName}.bam.multiple_metrics.base_distribution_by_cycle_metrics.table.txt"
		File gcdmNrm="${normalBamName}.bam.multiple_metrics.gc_bias.detail_metrics.table.txt"
		File gcSumNrm="${normalBamName}.bam.multiple_metrics.gc_bias.summary_metrics.table.txt"
		File ismNrm="${normalBamName}.bam.multiple_metrics.insert_size_metrics.table.txt"
		File oxoGNrm="${normalBamName}.bam.multiple_metrics.oxog_metrics.table.txt"
		File padmNrm="${normalBamName}.bam.multiple_metrics.pre_adapter_detail_metrics.table.txt"
		File pasmNrm="${normalBamName}.bam.multiple_metrics.pre_adapter_summary_metrics.table.txt"
		File qymNrm="${normalBamName}.bam.multiple_metrics.quality_yield_metrics.table.txt"
		File normalRGFPReport="${normalBamName}.normal.CCRGFPReport.txt"
		File normHSMetricsHist="${normalBamName}.normlHSMetrics.txt.hist.txt"
		File normHSMetricsTbl="${normalBamName}.normlHSMetrics.txt.table.txt"

		}
	}




workflow QC_Workflow {
	#RUNTIME INPUT PARAMS
	Int preemptible


	#WORKFLOW INPUT PARAMS
	File tumorBam
	File normalBam
	File tumorBamIdx
	File normalBamIdx
	File regionFile
	File readGroupBlackList
	File captureNormalsDBRCLZip
	String pairName
	File refFasta
	File refFastaIdx
	File refFastaDict
	File exomeIntervals
	File SNP6Bed
	File HapMapVCF
	File HaplotypeDBForCrossCheck
	String validationStringencyLevel
	File picardHapMap
	File picardBaitIntervals
	File picardTargetIntervals
	File DB_SNP_VCF
	File DB_SNP_VCF_IDX

	#setup for scatter
	call QC_Prepare_Task {
		input:
			preemptible=preemptible,
			tBamBytes=size(tumorBam),
			tBaiBytes=size(tumorBamIdx),
			nBamBytes=size(normalBam),
			nBaiBytes=size(normalBamIdx),
			regionFileBytes=size(regionFile),
			rgBLBytes=size(readGroupBlackList),
			capNormDBZipBytes=size(captureNormalsDBRCLZip),
			fastaBytes=size(refFasta),
			fastaDictBytes=size(refFastaDict),
			fastaIdxBytes=size(refFastaIdx),
			exomeIntervalsBytes=size(exomeIntervals),
			snpSixBytes=size(SNP6Bed),
			hapMapVCFBytes=size(HapMapVCF),
			hapDBForCCBytes=size(HaplotypeDBForCrossCheck),
			dbSNPVCFBytes=size(DB_SNP_VCF),
			dbSNPVCFIDXBytes=size(DB_SNP_VCF_IDX),
			picardHapMapVCFBytes=size(picardHapMap),
			picardTargetIntervalsBytes=size(picardTargetIntervals),
			picardBaitIntervalsBytes=size(picardBaitIntervals)
		}

	#scatter over chromosomes
	scatter(chrom in QC_Prepare_Task.split_chroms)
		{
		call QC_Scatter_Task {
			input:
				tumorBam=tumorBam,
				tumorBamIdx=tumorBamIdx,
				normalBam=normalBam,
				normalBamIdx=normalBamIdx,
				regionFile=regionFile,
				CHROM=chrom,
				preemptible=preemptible,
				diskGB=QC_Prepare_Task.qc_scatter_disk
			}
		}

	#gather over chromosomes
	call QC_Gather_Task {
		input:
			tumorBam=tumorBam,
			tumorBamIdx=tumorBamIdx,
			normalBam=normalBam,
			normalBamIdx=normalBamIdx,
			regionFile=regionFile,
			readGroupBlackList=readGroupBlackList,
			captureNormalsDBRCLZip=captureNormalsDBRCLZip,
			tumorRCLs=QC_Scatter_Task.tumorRCL,
			normalRCLs=QC_Scatter_Task.normalRCL,
			preemptible=preemptible,
			diskGB=QC_Prepare_Task.qc_gather_disk
		}


	call ContEST_Task {
		input : 
			tumorBam=tumorBam,
			tumorBamIdx=tumorBamIdx,
			normalBam=normalBam,
			normalBamIdx=normalBamIdx,
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			exomeIntervals=exomeIntervals,
			SNP6Bed=SNP6Bed,
			HapMapVCF=HapMapVCF,
			pairName=pairName,
			preemptible=preemptible,
			diskGB=QC_Prepare_Task.contest_disk
		}

	call CrossCheckLaneFingerprints_Task {
		input:
			tumorBam=tumorBam,
			normalBam=normalBam,
			tumorBamIdx=tumorBamIdx,
			normalBamIdx=normalBamIdx,
			pairName=pairName,
			HaplotypeDBForCrossCheck=HaplotypeDBForCrossCheck,
			validationStringencyLevel=validationStringencyLevel,
			preemptible=preemptible,
			diskGB=QC_Prepare_Task.cc_lane_fp_disk
		}

	#Picard tasks (tumor and normal)

	###################################
	#tumor
	call PicardMultipleMetrics_Task as tumorMM_Task {
		input:
			bam=tumorBam,
			bamIndex=tumorBamIdx,
			refFasta=refFasta,
			DB_SNP_VCF=DB_SNP_VCF,
			DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
			preemptible=preemptible,
			diskGB=QC_Prepare_Task.tumor_mm_disk
		}
	call BAMValidation_Task as tumorBAMValidation_Task {
		input:
			bam=tumorBam,
			bamIndex=tumorBamIdx,
			preemptible=preemptible,
			diskGB=QC_Prepare_Task.tumor_validation_disk
		}
	call crossCheckRGFP_Task as CCTumorRGFP_Task {
		input:
			bam=tumorBam,
			bamIndex=tumorBamIdx,
			preemptible=preemptible,
			diskGB=QC_Prepare_Task.tumor_rgfp_disk,
			picardHapMap=picardHapMap
		}
	call HSMetrics_Task as TumorHSMetrics_Task {
		input:
			bam=tumorBam,
			bamIndex=tumorBamIdx,
			preemptible=preemptible,
			diskGB=QC_Prepare_Task.tumor_hsm_disk,
			picardTargetIntervals=picardTargetIntervals,
			picardBaitIntervals=picardBaitIntervals
		}


	#####################################
	#normal
	call PicardMultipleMetrics_Task as normalMM_Task {
		input:
			bam=normalBam,
			bamIndex=normalBamIdx,
			refFasta=refFasta,
			DB_SNP_VCF=DB_SNP_VCF,
			DB_SNP_VCF_IDX=DB_SNP_VCF_IDX,
			preemptible=preemptible,
			diskGB=QC_Prepare_Task.nrml_mm_disk
		}
	call BAMValidation_Task as normalBAMValidation_Task {
		input:
			bam=normalBam,
			bamIndex=normalBamIdx,
			preemptible=preemptible,
			diskGB=QC_Prepare_Task.nrml_validation_disk
		}
	call crossCheckRGFP_Task as CCNormalRGFP_Task {
		input:
			bam=normalBam,
			bamIndex=normalBamIdx,
			preemptible=preemptible,
			diskGB=QC_Prepare_Task.nrml_rgfp_disk,
			picardHapMap=picardHapMap
		}
	call HSMetrics_Task as NormalHSMetrics_Task {
		input:
			bam=normalBam,
			bamIndex=normalBamIdx,
			preemptible=preemptible,
			diskGB=QC_Prepare_Task.nrml_hsm_disk,
			picardTargetIntervals=picardTargetIntervals,
			picardBaitIntervals=picardBaitIntervals
		}	

	###############################
	# Nozzle Summary report
	call QC_Nozzle_Report_Task {
		input:
			preemptible=preemptible,
			pairName=pairName,
			fracContam=ContEST_Task.fraction_contamination,
			CCTumorRGFPReport=CCTumorRGFP_Task.picardCCREport,
			CCNormalRGFPReport=CCNormalRGFP_Task.picardCCREport,
			tumorBAMValidation=tumorBAMValidation_Task.validation_verbose,
			normalBAMValidation=normalBAMValidation_Task.validation_verbose,
			picardTumorMultipleMetricsZip=tumorMM_Task.metricsReportsZip,
			picardNormalMultipleMetricsZip=normalMM_Task.metricsReportsZip,
			tumorHSMetrics=TumorHSMetrics_Task.hsMetrics,
			normalHSMetrics=NormalHSMetrics_Task.hsMetrics,
			copyNumQCZip=QC_Gather_Task.CopyNumQCOutZip,
			tumorBamName=tumorMM_Task.bamName,
			normalBamName=normalMM_Task.bamName,
			CCLaneFingerPrintsReportHTML=CrossCheckLaneFingerprints_Task.crossCheckReport,
			CCLaneFingerPrintsStats=CrossCheckLaneFingerprints_Task.crossCheckStats
		}

	output {
	       ContEST_Task.fraction_contamination
	       CrossCheckLaneFingerprints_Task.crossCheckStats
	       CrossCheckLaneFingerprints_Task.crossCheckReport
	       tumorMM_Task.metricsReportsZip
	       tumorBAMValidation_Task.validation_verbose
	       TumorHSMetrics_Task.hsMetrics
	       normalMM_Task.metricsReportsZip
	       normalBAMValidation_Task.validation_verbose
	       NormalHSMetrics_Task.hsMetrics
	       QC_Nozzle_Report_Task.nozzleHTML
	       QC_Nozzle_Report_Task.nozzleRData
	       QC_Nozzle_Report_Task.copyNumQCHTML
	}	

}