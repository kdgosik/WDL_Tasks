workflow CallSomaticMutations_Workflow {
	String pairID
	String normalName
	File normalBam
	File normalBai
	String tumorName
	File tumorBam
	File tumorBai
	File refFasta
	File refFastaIdx
	File refFastaDict
	File normalPanel
	Int downsampleToCoverage
	File readgroupBlacklist
	File mutectIntervalList
	File contEstIntervalList
	File dbSNPVCF
	File dbSNPVCFIdx
	File cosmicVCF
	File picard_hapMapVCF
	File oncotatorDataSourceTarGz
	String oncotatorMode
	String hdSize

	call PrepareScatter_Task {
		input:
			targetsIntervalList=mutectIntervalList
	}

	call ContEst_Task {
		input:
			pairID=pairID,
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			tumorBam=tumorBam,
			tumorBai=tumorBai,
			normalBam=normalBam,
			normalBai=normalBai,
			contEstIntervalList=contEstIntervalList,
			picard_hapMapVCF=picard_hapMapVCF,
			hdSize=hdSize
	}
    
    call ContEstFraction_Task {
        input:
            contEst=ContEst_Task.contEst
    }

	scatter (idx in PrepareScatter_Task.scatterIndices) {
		call Mutect1_Task {
			input:
				normalName=normalName,
				normalBam=normalBam,
				normalBai=normalBai,
				tumorName=tumorName,
				tumorBam=tumorBam,
				tumorBai=tumorBai,
				refFasta=refFasta,
				refFastaIdx=refFastaIdx,
				refFastaDict=refFastaDict,
				normalPanel=normalPanel,
				downsampleToCoverage=downsampleToCoverage,
				fractionContamination=ContEstFraction_Task.contEstFraction,
				readgroupBlacklist=readgroupBlacklist,
				targetsIntervalList=PrepareScatter_Task.interval_files[idx],
				dbSNPVCF=dbSNPVCF,
				dbSNPVCFIdx=dbSNPVCFIdx,
				cosmicVCF=cosmicVCF,
				hdSize=hdSize
		}		
	}

	call GatherMutect1_Task {
		input:
			mutect1_cw=Mutect1_Task.mutect1_cw,
			mutect1_pw=Mutect1_Task.mutect1_pw,
			mutect1_cs=Mutect1_Task.mutect1_cs
	}

	call summarizeWigFile_Task {
		input:
			pairID=pairID,
			wigFile=GatherMutect1_Task.mutect1_coveragewig
	}

	call CallStatstoMAFLite_Task {
		input:
			callstats=GatherMutect1_Task.mutect1_callstats,
			pairID=pairID
	}

	call Oncotator_Task {
		input:
			MAFLITE=CallStatstoMAFLite_Task.maflite,
			oncotatorDataSourceTarGz=oncotatorDataSourceTarGz,
			oncotatorMode=oncotatorMode,
			pairID=pairID
	}

	output {
		GatherMutect1_Task.mutect1_callstats
		GatherMutect1_Task.mutect1_powerwig
		GatherMutect1_Task.mutect1_coveragewig
		summarizeWigFile_Task.somatic_mutation_covered_bases_file_capture
		summarizeWigFile_Task.somatic_mutation_covered_bases_capture
		CallStatstoMAFLite_Task.maflite
		Oncotator_Task.oncotatedMAF
		ContEst_Task.contEst
		ContEst_Task.contEstTable
     	ContEst_Task.contEstBaseReport
     	ContEstFraction_Task.contEstFraction
	}
}

task ContEst_Task {
	String pairID
	File tumorBam
	File tumorBai
	File normalBam
	File normalBai
	File refFasta
	File refFastaIdx
	File refFastaDict
	File contEstIntervalList
	File picard_hapMapVCF
    String hdSize

    command <<<
    	java -Xmx4096m	-jar /usr/GenomeAnalysisTK.jar \
    	        -T ContEst \
    	        -R ${refFasta} -I:eval ${tumorBam} -I:genotype ${normalBam} \
    	        -l INFO -pf ${picard_hapMapVCF} \
    	        -o ${pairID}.contamination.txt -br ${pairID}.contEst.baseReport.txt \
    	        -L ${contEstIntervalList} -isr INTERSECTION \
    	        --minimum_base_count 100 --trim_fraction 0.03 --beta_threshold 0.05 \
    	        --min_genotype_depth 30 --min_genotype_ratio 0.8

    	awk 'FNR == 2 {print $4}' ${pairID}.contamination.txt > ${pairID}.contEst.txt
    >>>

    runtime {
		docker: "broadinstitute/gatk3:3.8-0"
		memory: "7GB"
		disks: "local-disk " + hdSize + " HDD"
    }

    output {
    	String contEst=read_string("${pairID}.contEst.txt") 
    	File contEstTable="${pairID}.contamination.txt"
    	File contEstBaseReport="${pairID}.contEst.baseReport.txt"
    }
}

task ContEstFraction_Task {
    String contEst

    command <<<
        python /contEstFraction.py --contEst ${contEst}
    >>>

    runtime {
        docker: "vanallenlab/contestfraction:1.0"
        memory: "1GB"
        disks: "local-disk 4 HDD"
    }

    output {
        String contEstFraction=read_string("contEstFraction.txt")
    }
}

task PrepareScatter_Task {
	File targetsIntervalList

	command <<<
		python3 /splitIntervals/splitIntervals.py --intervalHandle ${targetsIntervalList}

		# Create index for split intervals
		numIntervals="$(ls scatter.interval.* | wc -l)"
		seq 0 $((numIntervals - 1)) > indices.dat

	>>>

	output {
		Array[Int] scatterIndices=read_lines("indices.dat")
		Array[File] interval_files=glob("scatter.interval.*")
	}

	runtime {
		docker: "vanallenlab/mutect:1.1.6"
		memory: "1 GB"
	}
}

task Mutect1_Task {
	String normalName
	File normalBam
	File normalBai
	String tumorName
	File tumorBam
	File tumorBai
	File refFasta
	File refFastaIdx
	File refFastaDict
	File normalPanel
	Int downsampleToCoverage
	Float fractionContamination
	File readgroupBlacklist
	File targetsIntervalList
	File dbSNPVCF
	File dbSNPVCFIdx
	File cosmicVCF
	String hdSize

	command <<<
		java -jar -Xmx4g /muTect-1.1.6.jar --analysis_type MuTect -L ${targetsIntervalList} \
		--normal_sample_name ${normalName} -I:normal ${normalBam} \
		--tumor_sample_name ${tumorName} -I:tumor ${tumorBam} \
		--reference_sequence ${refFasta} --normal_panel ${normalPanel} \
		--fraction_contamination ${fractionContamination} --downsample_to_coverage ${downsampleToCoverage} \
		--dbsnp ${dbSNPVCF} --cosmic ${cosmicVCF} --read_group_black_list ${readgroupBlacklist} --enable_extended_output \
		--out Mutect1.call_stats.txt --coverage_file Mutect1.coverage.wig.txt --power_file Mutect1.power.wig.txt
	>>>

	runtime {
		docker: "vanallenlab/mutect:1.1.6"
		memory: "6.5 GB"
		disks: "local-disk " + hdSize + " HDD"
	}

	output {
		File mutect1_cs="Mutect1.call_stats.txt"
		File mutect1_pw="Mutect1.power.wig.txt"
		File mutect1_cw="Mutect1.coverage.wig.txt"
	}
}

task GatherMutect1_Task {
	Array[File] mutect1_cs
	Array[File] mutect1_pw
	Array[File] mutect1_cw

	command <<<
		MUTECT1_CW="Mutect1.coverage.wig.txt"
		cat ${sep = ' ' mutect1_cw} >> $MUTECT1_CW

		MUTECT1_PW="Mutect1.power.wig.txt"
		cat ${sep = ' ' mutect1_pw} >> $MUTECT1_PW

		MUTECT1_CS="Mutect1.call_stats.txt"
		head -2 ${mutect1_cs[0]} > $MUTECT1_CS
		cat ${sep = ' ' mutect1_cs} | grep -Pv '#' | grep -Pv '^contig' >> $MUTECT1_CS
	>>>
	
	runtime {
		docker: "vanallenlab/mutect:1.1.6"
		memory: "2 GB"
	}

	output {
        	File mutect1_coveragewig="Mutect1.coverage.wig.txt"
        	File mutect1_powerwig="Mutect1.power.wig.txt"
			File mutect1_callstats="Mutect1.call_stats.txt"
	}
}

task summarizeWigFile_Task {
	String pairID
	File wigFile

	command <<<
		python /home/summarizeWigFile.py --pairId ${pairID} --wigFile ${wigFile}
	>>>

	output {
		File somatic_mutation_covered_bases_file_capture="${pairID}.somatic_coverage_summary.txt"
		String somatic_mutation_covered_bases_capture=read_string("${pairID}.somatic_coverage_summary.txt")
	}

	runtime {
		docker: "breardon/summarizewigfile:1.0"
		memory: "1 GB"
		disks: "local-disk 2 HDD"
	}
}

task CallStatstoMAFLite_Task {
	File callstats
	String pairID

	command <<<
		build=37
		mode="FSTAR"
		f_threshold=0
		triallelic_mode="REJECT"
		extra_cols="tumor_f,init_t_lod,t_lod_fstar,t_alt_count,t_ref_count,judgement"
		output="${pairID}.maf"

		perl /home/call_stats_to_maflite.pl ${callstats} $build $mode $f_threshold $triallelic_mode $output $extra_cols
	>>>

	output {
		File maflite="${pairID}.maf"
	}

	runtime {
		docker: "breardon/callstatstomaflite:1.0"
		memory: "4 GB"
		disks: "local-disk 2 HDD"
	}
}

task Oncotator_Task {
	File MAFLITE
	File oncotatorDataSourceTarGz
	String oncotatorMode
	String pairID
	File defaultConfig="gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/reference/oncotator/tcgaMAFManualOverrides2.4.config"
	File uniProt="gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/reference/oncotator/tx_exact_uniprot_matches.txt"

	command <<<
		tar zxvf ${oncotatorDataSourceTarGz}
		dbDir=$(basename ${oncotatorDataSourceTarGz} .tar.gz)
		mv $dbDir oncotatorDatasourceDir

		/root/oncotator_venv/bin/oncotator -i MAFLITE -o TCGAMAF \
		--db-dir=oncotatorDatasourceDir/ --tx-mode=${oncotatorMode} \
		--default_config ${defaultConfig} -v ${MAFLITE} ${pairID}.oncotated.maf hg19 \
		--log_name oncotator_firehose.log --prepend --infer-onps -c ${uniProt} --collapse-number-annotations
	>>>

	runtime {
		docker: "broadinstitute/oncotator:1.9.3.0"
		memory: "4 GB"
		disks: "local-disk 50 HDD"
	}

	output {
		File oncotatedMAF="${pairID}.oncotated.maf"
	}
}
