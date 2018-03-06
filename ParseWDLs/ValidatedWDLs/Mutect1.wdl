workflow Mutect1_Workflow {
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

	call PrepareScatter_Task {
		input:
			targetsIntervalList=targetsIntervalList
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
				fractionContamination=fractionContamination,
				readgroupBlacklist=readgroupBlacklist,
				targetsIntervalList=PrepareScatter_Task.interval_files[idx],
				dbSNPVCF=dbSNPVCF,
				dbSNPVCFIdx=dbSNPVCFIdx,
				cosmicVCF=cosmicVCF
		}		
	}

	call GatherMutect1_Task {
		input:
			mutect1_cw=Mutect1_Task.mutect1_cw,
			mutect1_pw=Mutect1_Task.mutect1_pw,
			mutect1_cs=Mutect1_Task.mutect1_cs
	}

	output {
		GatherMutect1_Task.mutect1_callstats
		GatherMutect1_Task.mutect1_powerwig
		GatherMutect1_Task.mutect1_coveragewig
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
		memory: "7 GB"
		disks: "local-disk 100 HDD"
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