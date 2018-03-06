task CallSomaticMutations_131_Prepare_Task {
    #RUNTIME INPUT PARAMS
	Int preemptible
    Int memoryGB
	Int cpu
	Int diskGB	

	#TASK INPUT PARAMS
	File mutectIntervals 
	File refFastaIdx
	File refFasta
	File refFastaDict
	Int nWay
	Float intervalsSize
	Float nBamSize
	Float nBamIdxSize
	Float tBamSize
	Float tBamIdxSize
	Float fastaSize
	Float fastaIdxSize
	Float fastaDictSize
	Float dbSNPSize
	Float cosmicVCFSize
	Float normalPanelSize
	Float oncoDBTarBallSize
	Float oncoDBTarBallFactor
	Float VEP_File_Size
	Float VEP_File_SizeFactor
	Float fcFactor
	Float m1Factor
	Float gatherOncotateFactor
	Float vepFactor

	command <<<
		#increase verbosity
		set -x

		echo "TMPDIR IS $TMPDIR"
		mkdir -pv $TMPDIR

		# run the interval splitter
		/usr/local/bin/python /home/process_monitor.py /usr/local/bin/intervals_prep_script.sh ${nWay} ${mutectIntervals} ${refFasta}

		#disk size calculation for allocation
		#fasta refs size GB
		FASTA_SIZE_BYTES=`echo -ne "${fastaSize}\t${fastaIdxSize}\t${fastaDictSize}\n" | awk '{print $1+$2+$3}' | awk '{printf("%20.20f\n",$1)} ' `
		echo "FASTA_SIZE_BYTES IS $FASTA_SIZE_BYTES"

		#M1 task disk size
		M1_DISK=`echo -ne "$FASTA_SIZE_BYTES\t${intervalsSize}\t${nBamSize}\t${nBamIdxSize}\t${tBamSize}\t${tBamIdxSize}\t${dbSNPSize}\t${cosmicVCFSize}\t${normalPanelSize}\n" | awk '{print ($1+$2+$3+$4+$5+$6+$7+$8+$9)/1000000000}'| awk '{printf("%20.20f\n",$1)} '` ; 
		M1_DISK_WITH_FACTOR=`echo -ne "$M1_DISK\t${m1Factor}\t${tBamSize}\n" | awk '{ print 1+($1+($3/1000000000)*$2) }'| awk '{printf("%20.20f\n",$1)} '` ;
		echo "M1_DISK IS $M1_DISK"
		echo "M1_DISK_WITH_FACTOR IS $M1_DISK_WITH_FACTOR"

		#FC task disk size
		FC_DISK="$M1_DISK" ; 
		FC_DISK_WITH_FACTOR=`echo -ne "$M1_DISK\t${fcFactor}\t${tBamSize}\n" | awk '{ print 1+($1+($3/1000000000)*$2) }'| awk '{printf("%20.20f\n",$1)} '` ;
		echo "FC_DISK IS $FC_DISK"
		echo "FC_DISK_WITH_FACTOR IS $FC_DISK_WITH_FACTOR" ; 

		#oncotation disk size
		ONCO_DISK=`echo -ne "$FASTA_SIZE_BYTES\t${oncoDBTarBallSize}\n" | awk '{print ($1+($2*${oncoDBTarBallFactor}))/1000000000}'| awk '{printf("%20.20f\n",$1)} '` ; 
		ONCO_DISK_WITH_FACTOR=`echo -ne "$ONCO_DISK\t${gatherOncotateFactor}\t${tBamSize}\n" | awk '{ print 1+$1+($3/1000000000)*$2 }'| awk '{printf("%20.20f\n",$1)} '` ;
		echo "ONCO_DISK IS $ONCO_DISK"
		echo "ONCO_DISK_WITH_FACTOR IS $ONCO_DISK_WITH_FACTOR"

		#VEP disk sipze
		VEP_DISK=`echo -ne "$FASTA_SIZE_BYTES\t${VEP_File_Size}\n" | awk '{ print ($1+($2*${VEP_File_SizeFactor}))/1000000000 } ' | awk '{printf("%20.20f\n",$1)} '` ; 
		VEP_DISK_WITH_FACTOR=`echo -ne "$VEP_DISK\t${vepFactor}\t${tBamSize}\n" | awk '{ print 1+$1+($3/1000000000)*$2 }'| awk '{printf("%20.20f\n",$1)} '` ;
		echo "VEP_DISK IS $VEP_DISK" 
		echo "VEP_DISK_WITH_FACTOR IS $VEP_DISK_WITH_FACTOR"

		#write values to files for saving/usage
		echo $M1_DISK_WITH_FACTOR    |grep -Po '^\d+' > m1_disk_giga.dat ;
		echo $FC_DISK_WITH_FACTOR    |grep -Po '^\d+' > fc_disk_giga.dat ;
		echo $ONCO_DISK_WITH_FACTOR  |grep -Po '^\d+' > onco_disk_giga.dat ;
		echo $VEP_DISK_WITH_FACTOR   |grep -Po '^\d+' > vep_disk_giga.dat ;


		>>>

	output  {
		Int m1DiskGiga=read_int("m1_disk_giga.dat") 
		Int fcDiskGiga=read_int("fc_disk_giga.dat")
		Int oncoDiskGiga=read_int("onco_disk_giga.dat")
		Int vepDiskGiga=read_int("vep_disk_giga.dat")
		File dstat_log="dstat.log.txt"
		File df_log="df.log.txt"
		Array[File] interval_files=glob("split_base.*")
		}

	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta:benchmark_2"
		memory: "${memoryGB} GB"
		disks: "local-disk  ${diskGB}  HDD"
		cpu: "${cpu}"
		}
	}




task Mutect1_Task {
	#RUNTIME INPUT PARAMS
	Int preemptible
	Int memoryGB
	Int cpu
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
	/usr/local/bin/python /home/process_monitor.py /usr/local/bin/java -jar -Xmx4g /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
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
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta:benchmark_2"
		memory: "${memoryGB} GB"
		disks: "local-disk ${diskGB} HDD"
		cpu: "${cpu}"
		}


	output {
		File dstat_log="dstat.log.txt"
		File df_log="df.log.txt"
		File mutect1_cs="MuTect1.call_stats.txt"
		File mutect1_pw="MuTect1.power.wig.txt"
		File mutect1_cw="MuTect1.coverage.wig.txt"
		}

	}


task MutectFC_Task {
	#RUNTIME INPUT PARAMS
	Int preemptible
    Int memoryGB
	Int cpu
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


	/usr/local/bin/python /home/process_monitor.py /usr/local/bin/java -jar -Xmx4g /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
	 -L ${mutectIntervals}  --normal_sample_name NORMAL_SAMPLE -I:normal  ${normalBam}  \
	 --tumor_sample_name TUMOR_SAMPLE -I:tumor ${tumorBam}  \
	 --reference_sequence ${refFasta} \
	 --fraction_contamination ${fracContam}  --dbsnp ${dbSNPVCF} \
	 --cosmic ${cosmicVCF} \
	 --force_output \
	 --out MuTectFC.call_stats.txt --coverage_file MuTectFC.coverage.wig.txt \
	 --power_file MuTectFC.power.wig.txt --downsample_to_coverage ${downsampleToCoverage} \
	 $NORMAL_PANEL_FLAG_AND_VAL

	>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta:benchmark_2"
		memory: "${memoryGB} GB"
		disks: "local-disk  ${diskGB}  HDD"
		cpu: "${cpu}"
		}

	output {
		File dstat_log="dstat.log.txt"
		File df_log="df.log.txt"
		File mutectfc_cs="MuTectFC.call_stats.txt"
		File mutectfc_pw="MuTectFC.power.wig.txt"
		File mutectfc_cw="MuTectFC.coverage.wig.txt"
		}

	}

task GatherAndOncotate_Task {
    #RUNTIME INPUT PARAMS
	Int preemptible
	Int memoryGB
	Int cpu
	Int diskGB

	#TASK INPUT PARAMS
	Array[File] mutect1_cs
	File oncoDBTarBall
	File VEP_File


	command <<<
		#increase verbosity
		set -x

		#make M1 FOF
		find ${sep =' ' mutect1_cs} > m1_cs.fof ;

		#run oncotation!
		/usr/local/bin/python /home/process_monitor.py /usr/local/bin/oncotation.sh m1_cs.fof ${oncoDBTarBall} 

		>>>


	output {
		#Array[File] all_outs=glob("*")
		File dstat_log="dstat.log.txt"
		File df_log="df.log.txt"
		File preannotated_vcf="MuTect1.call_stats.vcf.filtered.vcf"
		File oncotator_log="oncotator.log"
		File oncotator_out="MuTect1.call_stats.vcf.filtered.vcf.annotated.vcf"
		File mutect1Merged="MuTect1.call_stats.txt"
		File VEP_VCF="MuTect1.call_stats.vcf"
		}


	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta:benchmark_2"
		memory: "${memoryGB} GB"
		disks: "local-disk  ${diskGB}  HDD"
		cpu: "${cpu}"
		}
	}

task VEP_Task {

    #RUNTIME INPUT PARAMS
	Int preemptible
	Int memoryGB
	Int cpu
	Int diskGB

	#TASK INPUT PARAMS
	File mutectMergedRawVCF
	File VEP_File

	command <<<

		#increase verbosity
		set -x

		#run VEP 
		/usr/local/bin/python /home/process_monitor.py /usr/local/bin/vep_run.sh ${VEP_File} ${mutectMergedRawVCF}

		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "broadinstitute/broadmutationcalling_beta:benchmark_2"
		memory: "${memoryGB} GB"
		disks: "local-disk  ${diskGB}  HDD"
		cpu: "${cpu}"
		}

	output {
		File dstat_log="dstat.log.txt"
		File df_log="df.log.txt"
		File VEP_Output="variant_effect_output.txt"
		File VEP_Report="variant_effect_output.txt_summary.html"
		}


	}





task mutect_nozzle_task {
	#RUNTIME INPUT PARAMS
	Int preemptible
	Int memoryGB
	Int cpu
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

		#make the nozzle!
		/usr/bin/python /home/process_monitor.py /usr/local/bin/nozzle_report_wrapped.sh ${vcfFile} ${vepVCFFile} ${vepReportHTML} ${fcCallStats} ${pairName}

		>>>
		

	runtime {
		docker: "broadinstitute/broadmutationcalling_beta:benchmark_nozzle_1"
		disks: "local-disk  ${diskGB}  HDD"
		memory: "${memoryGB} GB"
		preemptible: "${preemptible}"
		cpu: "${cpu}"
		}


	output {
		#process_monitor.py scripts
		File dstat_log="dstat.log.txt"
		File df_log="df.log.txt"		

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
			intervalsSize=size(mutectIntervals),
			nBamSize=size(normalBam),
			nBamIdxSize=size(normalBamIdx),
			tBamSize=size(tumorBam),
			tBamIdxSize=size(tumorBamIdx),
			fastaSize=size(refFasta),
			fastaIdxSize=size(refFastaIdx),
			fastaDictSize=size(refFastaDict),
			dbSNPSize=size(dbSNPVCF),
			cosmicVCFSize=size(cosmicVCF),
			normalPanelSize=size(normalPanel),
			oncoDBTarBallSize=size(oncoDBTarBall),
			VEP_File_Size=size(VEP_File),
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			mutectIntervals=mutectIntervals,
			nWay=nWay
		}

	#SCATTER AND ANALYZE
	scatter (interval_file in CallSomaticMutations_131_Prepare_Task.interval_files) {
			
			call Mutect1_Task {
				input:
					tumorBam=tumorBam,
					normalBam=normalBam,
					tumorBamIdx=tumorBamIdx,
					normalBamIdx=normalBamIdx,
					mutectIntervals=interval_file,
					refFasta=refFasta,
					refFastaIdx=refFastaIdx,
					refFastaDict=refFastaDict,
					fracContam=fracContam,
					dbSNPVCF=dbSNPVCF,
					cosmicVCF=cosmicVCF,
					downsampleToCoverage=downsampleToCoverage,
					readGroupBlackList=readGroupBlackList,
					normalPanel=normalPanel,
					diskGB=CallSomaticMutations_131_Prepare_Task.m1DiskGiga
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
			diskGB=CallSomaticMutations_131_Prepare_Task.fcDiskGiga
			}

	call GatherAndOncotate_Task {
		input:
			mutect1_cs=Mutect1_Task.mutect1_cs,
			oncoDBTarBall=oncoDBTarBall,
			VEP_File=VEP_File,
			diskGB=CallSomaticMutations_131_Prepare_Task.oncoDiskGiga
		}

	call VEP_Task {
		input:
			mutectMergedRawVCF=GatherAndOncotate_Task.VEP_VCF,
			VEP_File=VEP_File,
			diskGB=CallSomaticMutations_131_Prepare_Task.vepDiskGiga
		}

	call mutect_nozzle_task {
		input:
			vcfFile=GatherAndOncotate_Task.oncotator_out,
			vepReportHTML=VEP_Task.VEP_Report,
			vepVCFFile=VEP_Task.VEP_Output,
			pairName=pairName,
			fcCallStats=MutectFC_Task.mutectfc_cs,
			diskGB=CallSomaticMutations_131_Prepare_Task.m1DiskGiga
		}

	
	output {
	       MutectFC_Task.mutectfc_cs
	       MutectFC_Task.mutectfc_pw
	       MutectFC_Task.mutectfc_cw
	       GatherAndOncotate_Task.preannotated_vcf
	       GatherAndOncotate_Task.oncotator_log
	       GatherAndOncotate_Task.mutect1Merged
	       VEP_Task.VEP_Output
	       VEP_Task.VEP_Report
	       mutect_nozzle_task.nozzleReport
	       mutect_nozzle_task.nozzleRData
	       }

	}