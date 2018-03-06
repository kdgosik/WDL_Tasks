workflow contEst_workflow {
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
            picard_hapMapVCF=picard_hapMapVCF
        }

    call ContEstFraction_Task {
        input:
            contEst=ContEst_Task.contEst
    }

     output {
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
		disks: "local-disk 100 HDD"
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
