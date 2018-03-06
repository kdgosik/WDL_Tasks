task ContEstTask {
  File normalBamHG19
  File normalBamIndexHG19
  File tumorBamHG19
  File tumorBamIndexHG19
  File ReferenceFasta
  File ContESTIntervals
  File HapMapVCF
  File SNP6Bed
  command <<<
    #echo contest input
    echo "Contest setup files are :"
    find `pwd`
    echo -ne "\n\n"
    for BAM in ${normalBamHG19} ${tumorBamHG19} ; do
    	echo "GOT BAM FILE $BAM" ;
    	DN=`dirname $BAM` ;
    	echo "Its directory is $DN" ;
    	echo "The files there are : " ;
    	ls -alh $DN ;
    done ;

    ################################################################
    #ContEST
    /cga/fh/pcawg_pipeline/utils/firehose_module_adaptor/run_module.py --module_libdir /cga/fh/pcawg_pipeline/modules/contest --reference.file ${ReferenceFasta} --intervals ${ContESTIntervals} --sample.id SAMPLE --clean.bam.file ${tumorBamHG19}  --normal.bam  ${normalBamHG19} --genotypes.vcf none --pop.vcf ${HapMapVCF} --force.array.free true --snp.six.bed ${SNP6Bed} --job.spec.memory 2
    #Use proper units by dividing here. Save to contamination.dat
    echo `cat SAMPLE.contamination.txt.firehose`/100.0 | bc -l > contamination.dat
    #echo contest input
    echo "Contest ending files are :"
    find `pwd`
    echo -ne "\n\n"
    #file inspection
    #Note useage of F instead of usual braces with xargs
    find `pwd` -type f | grep -Piv '\.ba[mi]' | xargs -tI F head  --verbose --lines=250 F
  >>>

  runtime {
    docker: "eddiebroad/contest_mutect_oncotator"
    memory: "24 GB"
    disks: "local-disk 100 SSD"
    }

  output {
    File contaminationFile = "contamination.dat"
    File outFile = "stdout.txt"
    File errFile = "stderr.txt"
  }
}

task MutectTask {
  File HapMapVCF
  File MutectIntervals
  File DBSNPVCF
  File normalBamHG19
  File normalBamIndexHG19
  File tumorBamHG19
  File tumorBamIndexHG19
  File ReferenceFasta
  File ReferenceFastaIndex
  File ReferenceFastaDict
  File COSMICVCF
  File ContEstContamination
  command <<<
    ################################################################
    #MuTect
    export FRAC_CONTAM=`cat ${ContEstContamination}`;
    #Read contamination fraction from ContEst and make a validity flag
    export CONTAM_FLAG=`echo "$FRAC_CONTAM>=0"|bc` ;
    if [[ "$CONTAM_FLAG" -ge 1 ]]; then
    	export MUTECT_FRAC_CONTAM="$FRAC_CONTAM" ;
    	echo "Using the ContEst computed fraction contamination value of $MUTECT_FRAC_CONTAM"
    else
    	export MUTECT_FRAC_CONTAM=0.02 ;
    	echo "Using default contamination value of $MUTECT_FRAC_CONTAM due to negative ContEST return value" ;
    fi ;
    java -jar -Xmx4g /cga/fh/pcawg_pipeline/modules/mutect1/muTect-1.1.6.jar --analysis_type MuTect  -L ${MutectIntervals}   --normal_sample_name NORMAL_SAMPLE -I:normal ${normalBamHG19} --tumor_sample_name TUMOR_SAMPLE -I:tumor ${tumorBamHG19}  --reference_sequence  ${ReferenceFasta}   --fraction_contamination $MUTECT_FRAC_CONTAM --dbsnp ${DBSNPVCF} --cosmic ${COSMICVCF} --out MuTect.call_stats.txt --coverage_file MuTect.coverage.wig.txt --power_file MuTect.power.wig.txt --downsample_to_coverage 10000 ;
    ################################################################
    #call stats to MAF lite
    perl  /cga/fh/pcawg_pipeline/modules/callstats_to_maflite/call_stats_to_maflite.pl MuTect.call_stats.txt 37  FSTAR 0 REJECT MuTect.call_stats.maf tumor_f,init_t_lod,t_lod_fstar,t_alt_count,t_ref_count,judgement
  >>>
  output {
    File MAFLiteFile = "MuTect.call_stats.maf"
    File CallStatsFile = "MuTect.call_stats.txt"
  }
  runtime {
    docker: "eddiebroad/contest_mutect_oncotator"
    memory: "24 GB"
    disks: "local-disk 100 SSD"
    }  
}

task OncotatorTask {
  File inputMAFLite
  command <<<
    ################################################################
    #Oncotator
    oncotator --db-dir /cga/fh/pcawg_pipeline/modules/oncotator/empty  -o VCF ${inputMAFLite} SAMPLE.vcf hg19
  >>>

  output {
    File oncotatorVCF = "SAMPLE.vcf"
    File oncotatorLog = "oncotator.log"
  }

  runtime {
    docker: "eddiebroad/contest_mutect_oncotator"
    memory: "24 GB"
    disks: "local-disk 100 SSD"
    }
}

task OncotatorReportTask {
  File oncoVCFFile
  command <<<
    echo "<HTML><HEADER><TITLE>Oncotator VCF</TITLE></HEADER><BODY><H2>Oncotator Report</H2><TABLE BORDER='1'><TR>" > oncotator_out.html
    #make table column headers
    for HEADER in `grep -P '#CHROM' ${oncoVCFFile} |sed -r "s/^#//g"` ; do echo "<TH>$HEADER</TH>"  >> oncotator_out.html      ; done ;
    echo "</TR>" >> oncotator_out.html
    #get VCF chrom ROW
    HR=`grep -Pin '^#CHROM' ${oncoVCFFile} |grep -Po '^\d+'`
    #get num lines
    VCFLINES=`wc -l ${oncoVCFFile} |grep -Po '^\d+'`
    #get num lines for tail
    TAIL_LINES=`echo $VCFLINES-$HR|bc`;
    #make data rows
    tail --lines=$TAIL_LINES  ${oncoVCFFile}  | sed -r "s/\t/<\/TD><TD>/g"|sed -r "s/^/<TR><TD>/g"|sed -r "s/$/<\/TD><\/TR>/g" >> oncotator_out.html
    #close up TABLE
    echo "</TABLE></BODY></HTML>" >> oncotator_out.html
  >>>
  output {
    File oncotatorHTMLReport = "oncotator_out.html"
  }

  runtime {
    docker: "eddiebroad/contest_mutect_oncotator"
    memory: "24 GB"
    disks: "local-disk 100 SSD"
    }
}

workflow ContEstMuTectOncotatorWorkflow {
  File normalBamHG19
  File normalBamIndexHG19
  File tumorBamHG19
  File tumorBamIndexHG19
  File ReferenceFasta
  File ContESTIntervals
  File HapMapVCF
  File SNP6Bed
  File MutectIntervals
  File DBSNPVCF
  File COSMICVCF
  File ReferenceFastaIndex
  File ReferenceFastaDict
  call ContEstTask {
    input: tumorBamHG19=tumorBamHG19, HapMapVCF=HapMapVCF, ReferenceFasta=ReferenceFasta,
					 ContESTIntervals=ContESTIntervals, normalBamIndexHG19=normalBamIndexHG19,
					 normalBamHG19=normalBamHG19, tumorBamIndexHG19=tumorBamIndexHG19, SNP6Bed=SNP6Bed
  }
  call MutectTask {
    input: tumorBamHG19=tumorBamHG19, HapMapVCF=HapMapVCF, ReferenceFastaIndex=ReferenceFastaIndex,
					 MutectIntervals=MutectIntervals, DBSNPVCF=DBSNPVCF, COSMICVCF=COSMICVCF, ReferenceFasta=ReferenceFasta,
					 normalBamIndexHG19=normalBamIndexHG19, normalBamHG19=normalBamHG19, tumorBamIndexHG19=tumorBamIndexHG19,
					 ContEstContamination=ContEstTask.contaminationFile, ReferenceFastaDict=ReferenceFastaDict
  }
  call OncotatorTask {
    input: inputMAFLite=MutectTask.MAFLiteFile
  }
  call OncotatorReportTask {
    input: oncoVCFFile=OncotatorTask.oncotatorVCF
  }

}