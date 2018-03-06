task QC_Task {
	
	File tumorBam
	File tumorBamIdx
	File normalBam
	File normalBamIdx
	File regionFile
	File readGroupBlackList
	File captureNormalsDBRCLZip

	command <<<
	#increase verbosity
	set -x

	###############################################
	# Copy Number QC for Capture
	#Make lane lists for tumor and normal
	#MakeLaneList_12
	java -jar /usr/local/bin/MakeLaneList.jar ${tumorBam}  case.lanelist.txt ;
	java -jar /usr/local/bin/MakeLaneList.jar ${normalBam} control.lanelist.txt ; 

	#make region coverages per lane for tumor
	#RegionCovPerLane_18
	for CHROM in `seq 1 24` ; 
		do 
		echo "Working on $CHROM for ${tumorBam}"
		SCATTER_DIR="scatter.$CHROM" ; 
		mkdir -v $SCATTER_DIR
		OUTPATH=`echo -ne "$SCATTER_DIR/chr$CHROM.case.rcl"` ; 
		echo "OUTPATH is $OUTPATH"
		echo "chrom is $CHROM"
		java -jar /usr/local/bin/RegionCovPerLane.jar ${tumorBam} ${regionFile} $OUTPATH $CHROM   
	done ;
	wait ; 
	/usr/local/bin/rcl.gather.sh  case.rcl ; 
	rm -rf scatter.* ;

	#make region coverages per lane for control
	#RegionCovPerLane_18
	for CHROM in `seq 1 24` ; 
		do 
		echo "Working on $CHROM for ${normalBam}"
		SCATTER_DIR="scatter.$CHROM" ; 
		mkdir -v $SCATTER_DIR
		OUTPATH=`echo -ne "$SCATTER_DIR/chr$CHROM.control.rcl"` ; 
		echo "OUTPATH is $OUTPATH"
		echo "chrom is $CHROM"
		java -jar /usr/local/bin/RegionCovPerLane.jar ${normalBam} ${regionFile} $OUTPATH $CHROM   
	done ;
	wait ; 
	/usr/local/bin/rcl.gather.sh  control.rcl ; 
	rm -rf scatter.* ;

	#run the mat-lab based report
	#CopyNumberQCReport_27
	cp -vfr /CopyNumberQCReport_27/unzip/* . 
	#Make a file of files paths
	#This command aims to list the zip files contents, filtering out all but the file paths and then write the paths to a file
	unzip -l ${captureNormalsDBRCLZip} |sed -r "s/^\s+//g"|sed -r "s/^[0-9]+//g"|sed -r "s/^\s*//g"|sed -r "s/^\S*//g"|sed -r "s/^\s*[0-9:]*\s*//g"|grep -Pv '^Date\s+Time\s+Name'|grep -Pv '^[\-\s]*$'|grep -Pv '\.zip$'|grep -Pv '^files$' > capture_normals_db_wdl
	#this command actually performs the unzipping
	unzip ${captureNormalsDBRCLZip}
	./run_fh_CopyNumberQCReport.sh /opt/MATLAB/MATLAB_Compiler_Runtime/v710/  PairCopyNumQCReport  \
	   case.rcl  control.rcl  case.lanelist.txt  control.lanelist.txt \
	     ${readGroupBlackList}  ${regionFile} capture_normals_db_wdl  NA NA
	#zip up the output
	zip CopyNumQCout.zip report.html num.mixups.txt PairCopyNumQCReport* *.png


	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta"
		memory: "7 GB"
		disks: "local-disk 100 HDD"
		preemptible: 1
		}

	output {
		#lane lists
		File tumorBamLaneList="case.lanelist.txt"
		File normalBamLaneList="control.lanelist.txt"
		File tumorRCL="case.rcl"
		File normalRCL="control.rcl"
		#CopyNumQC
		File CopyNumQCOutZip="CopyNumQCout.zip"
		File CopyNumQCReport="report.html"
		File CopyNumQCReportPNG="PairCopyNumQCReport_CopyNumberQC.png"
		File CopyNumQCMixUps="num.mixups.txt"
		}
	}



task PicardMultipleMetrics {

	File bam
	File bamIndex
	File refFasta
	File HapMapVCF
	File picardTargetIntervals
	File picardHapMap
	File picardBaitIntervals
	File DB_SNP_VCF
	File DB_SNP_VCF_IDX

	command <<<
	#increse verbosity
	set -x

	for BAM in ${bam} ; 
		do
		BASE=`basename $BAM` ; 
		echo "Processing $BAM ..." ; 
		/usr/local/jre1.8.0_73/bin/java   -Xmx4g  -jar  /usr/local/bin/picardtools-2.1.0.jar   CollectMultipleMetrics \
			CREATE_MD5_FILE=true   I=$BAM O=$BASE.multiple_metrics R=${refFasta} PROGRAM=CollectAlignmentSummaryMetrics \
			PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution  PROGRAM=MeanQualityByCycle \
			PROGRAM=CollectBaseDistributionByCycle PROGRAM=CollectSequencingArtifactMetrics  \
			PROGRAM=CollectQualityYieldMetrics PROGRAM=CollectGcBiasMetrics ;
		done ;

	#zip up reports
	BASE=`basename ${bam}` ; 
	zip piccard_multiple_metrics.zip $BASE.multiple_metrics.*

	#collect oxog metrics
	/usr/local/jre1.8.0_73/bin/java  -Xmx4g  -jar /usr/local/bin/picardtools-2.1.0.jar CollectOxoGMetrics \
	I=${bam} \
	O=oxoG_metrics.txt \
	R=${refFasta}

	#validation summary
	/usr/local/jre1.8.0_73/bin/java  -Xmx4g  -jar /usr/local/bin/picardtools-2.1.0.jar ValidateSamFile \
		I=${bam} MODE=SUMMARY > validation_summary.txt  2>&1 
	#validation verbose
	/usr/local/jre1.8.0_73/bin/java  -Xmx4g  -jar /usr/local/bin/picardtools-2.1.0.jar ValidateSamFile \
		I=${bam} MODE=VERBOSE > validation_verbose.txt 2>&1 

	#cross check
	#Exception in thread "main" java.lang.IllegalStateException: Haplotype map file must contain header
	/usr/local/jre1.8.0_73/bin/java  -Xmx4g  -jar /usr/local/bin/picardtools-2.1.0.jar \
	CrosscheckReadGroupFingerprints HAPLOTYPE_MAP=${picardHapMap}  I=${bam} O="picard_crosscheck_report.txt"

	#HS metrics
	/usr/local/jre1.8.0_73/bin/java  -Xmx4g  -jar /usr/local/bin/picardtools-2.1.0.jar CollectHsMetrics I=${bam} \
	 BAIT_INTERVALS=${picardBaitIntervals} TARGET_INTERVALS=${picardTargetIntervals} OUTPUT="HSMetrics.txt"

	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta"
		memory: "7 GB"
		disks: "local-disk 100 HDD"
		preemptible: 1
		}


	output {
		#File fpDetails="fingerprinting_detail_metrics.txt"
		#File fpSummary="fingerprinting_summary_metrics.txt"
		File hsMetrics="HSMetrics.txt"
		File picardCCREport="picard_crosscheck_report.txt"
		File validation_summary="validation_summary.txt"
		File validation_verbose="validation_verbose.txt"
		File metricsReportsZip="piccard_multiple_metrics.zip"
		File oxoG_metrics="oxoG_metrics.txt"

		}

	}





task ContESTTask {
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
	java  -Xmx1024m  -Djava.io.tmpdir=/tmp -cp /usr/local/bin/Queue-1.4-437-g6b8a9e1-svn-35362.jar \
	org.broadinstitute.sting.gatk.CommandLineGATK  -T Contamination  -I:eval ${tumorBam} \
	-I:genotype ${normalBam} -et NO_ET -L ${exomeIntervals} -L ${SNP6Bed}   -isr INTERSECTION  \
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
		docker: "broadinstitute/broadmutationcalling_qc_beta"
		memory: "7 GB"
		disks: "local-disk 100 HDD"
		preemptible: 1
		}


	output {
		File contamDataFile="contamination_validation.array_free.txt"
		File contestAFFile="contamination.af.txt"
		File contestBaseReport="contamination.base_report.txt"
		File validationOutput="contest_validation.output.tsv"
		String fraction_contamination=read_string("fraction_contamination.txt")
		}

	}


task CrossCheckLaneFingerprints {

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
		
		#Cross Check!
		#CrosscheckLaneFingerprints[version=9]
		mkdir -v tmp
		java -Xmx3g -jar /usr/local/bin/CrosscheckReadGroupFingerprints.jar \
		 OPTIONS_FILE=$OPTIONS_FILE  H=${HaplotypeDBForCrossCheck} TMP_DIR=`pwd`/tmp QUIET=true \
		   EXIT_CODE_WHEN_MISMATCH=0 OUTPUT=crosscheck.stats.txt  \
		    VALIDATION_STRINGENCY=${validationStringencyLevel}

		#CrosscheckLaneFingerprintsReport[version=10]
		#Run report
		/usr/local/bin/crosscheck_report.pl --stats crosscheck.stats.txt ${normalBam} ${tumorBam} conv_file.input.tsv 

	>>>


	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta"
		memory: "7 GB"
		disks: "local-disk 100 HDD"
		preemptible: 1
		}


	output {
		File crossCheckStats="crosscheck.stats.txt"
		File crossCheckReport="report.html"
		}
	}



task CollectSequencingArtifactMetrics {

	File bam
	File bamIndex
	File refFasta
	File DB_SNP_VCF
	File DB_SNP_VCF_IDX

	command <<<
		set -x

		java -Xmx3600M -jar  /usr/local/bin/picard.1.895.jar CollectSequencingArtifactMetrics \
		DB_SNP=${DB_SNP_VCF} INPUT=${bam} OUTPUT="sequencing_artifact_metrics.txt"  REFERENCE_SEQUENCE=${refFasta}  MINIMUM_QUALITY_SCORE=20 \
		MINIMUM_MAPPING_QUALITY=30 MINIMUM_INSERT_SIZE=60 MAXIMUM_INSERT_SIZE=600 INCLUDE_UNPAIRED=false \
		TANDEM_READS=false USE_OQ=true CONTEXT_SIZE=1 ASSUME_SORTED=true STOP_AFTER=100000000 VERBOSITY=INFO \
		QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false \
		CREATE_MD5_FILE=false


		>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta"
		memory: "7 GB"
		disks: "local-disk 100 HDD"
		preemptible: 1
		}


	output {
		File SAM_Bait_Detail="sequencing_artifact_metrics.txt.bait_bias_detail_metrics"
		File SAM_Bait_Summary="sequencing_artifact_metrics.txt.bait_bias_summary_metrics"
		File SAM_PreAdapter_Detail="sequencing_artifact_metrics.txt.pre_adapter_detail_metrics"
		File SAM_PreAdapter_Summary="sequencing_artifact_metrics.txt.pre_adapter_summary_metrics"
		}
	}




workflow QC_Workflow {

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

	call QC_Task {
		input : 
			tumorBam=tumorBam,
			normalBam=normalBam,
			tumorBamIdx=tumorBamIdx,
			normalBamIdx=normalBamIdx,
			regionFile=regionFile,
			readGroupBlackList=readGroupBlackList,
			captureNormalsDBRCLZip=captureNormalsDBRCLZip
		}


	call ContESTTask {
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
			pairName=pairName
		}

	call CrossCheckLaneFingerprints {
		input:
			tumorBam=tumorBam,
			normalBam=normalBam,
			tumorBamIdx=tumorBamIdx,
			normalBamIdx=normalBamIdx,
			pairName=pairName,
			HaplotypeDBForCrossCheck=HaplotypeDBForCrossCheck,
			validationStringencyLevel=validationStringencyLevel
		}


	#Piccard
	call PicardMultipleMetrics as tumorMM {
		input:
			bam=tumorBam,
			bamIndex=tumorBamIdx,
			refFasta=refFasta,
			HapMapVCF=HapMapVCF,
			picardHapMap=picardHapMap,
			picardBaitIntervals=picardBaitIntervals,
			picardTargetIntervals=picardTargetIntervals,		
			DB_SNP_VCF=DB_SNP_VCF,
			DB_SNP_VCF_IDX=DB_SNP_VCF_IDX
		}
	call PicardMultipleMetrics as normalMM {
		input:
			bam=normalBam,
			bamIndex=normalBamIdx,
			refFasta=refFasta,
			HapMapVCF=HapMapVCF,
			picardHapMap=picardHapMap,
			picardBaitIntervals=picardBaitIntervals,
			picardTargetIntervals=picardTargetIntervals,		
			DB_SNP_VCF=DB_SNP_VCF,
			DB_SNP_VCF_IDX=DB_SNP_VCF_IDX
		}

	call CollectSequencingArtifactMetrics as tumorArtifacts {
		input:
			bam=tumorBam,
			bamIndex=tumorBamIdx,
			refFasta=refFasta,		
			DB_SNP_VCF=DB_SNP_VCF,
			DB_SNP_VCF_IDX=DB_SNP_VCF_IDX
		}
	call CollectSequencingArtifactMetrics as normalArtifacts {
		input:
			bam=normalBam,
			bamIndex=normalBamIdx,
			refFasta=refFasta,
			DB_SNP_VCF=DB_SNP_VCF,
			DB_SNP_VCF_IDX=DB_SNP_VCF_IDX
		}


	}