
workflow gatk_alignbam {

    call gatk_alignbam_task 
}



task gatk_alignbam_task {
    String picard = "/cil/shed/apps/external/picard/current/bin/picard.jar"
    File in_bam
    String sample_name
    String fq1 = "${sample_name}.1.fq"
    String fq2 = "${sample_name}.2.fq"
    String aligned_bam = "${sample_name}.aligned.bam"
    String sorted_bam = "${sample_name}.sorted.bam"
    String marked_bam = "${sample_name}.marked_duplicates.bam"
    String marked_duplicates_metrics = "${sample_name}.marked_duplicates.metrics"
    String out_bam_path = "${sample_name}.bam"
    String out_index_path = "${out_bam_path}.bai"
    String ref_fasta = "ref.fasta"
    File reference_tgz

    String read_group = "\\'@RG\\\\\\\\tID:FLOWCELL_${sample_name}\\\\\\\\tSM:${sample_name}\\\\\\\\tPL:ILLUMINA\\\\\\\\tLB:LIB_${sample_name}\\'"

    String output_disk_gb 
    String boot_disk_gb = "10"
    String ram_gb = "10"
    String cpu_cores = "1"
    String preemptible = "0"
    String debug_dump_flag

    command {
        set -euo pipefail
        ln -sT `pwd` /opt/execution
        ln -sT `pwd`/../inputs /opt/inputs

        /opt/src/algutil/monitor_start.py
python_cmd="
import subprocess
def run(cmd):
    print (cmd)
    subprocess.check_call(cmd,shell=True)

run('echo STARTING tar xvf to unpack reference')
run('date')
run('tar xvf ${reference_tgz}')

run('echo STARTING SamToFastq')
run('date')
run('java -Xmx12G -jar ${picard} SamToFastq INPUT=${in_bam} FASTQ=${fq1} SECOND_END_FASTQ=${fq2} VALIDATION_STRINGENCY=LENIENT')

run('echo STARTING bwa mem')
run('date')
run('bwa mem -t 8 -R ${read_group} ${ref_fasta} ${fq1} ${fq2} | samtools view -bS - > ${aligned_bam}')

run('echo STARTING SortSam')
run('date')
run('java -Xmx8G -jar ${picard} SortSam I=${aligned_bam} O=${sorted_bam} SO=coordinate')

run('echo STARTING MarkDuplicates')
run('date')
run('java -Xmx8G -jar ${picard} MarkDuplicates I=${sorted_bam} O=${marked_bam} M=${marked_duplicates_metrics}')

run('echo STARTING ReorderSam')
run('date')
run('java -Xmx8G -jar ${picard} ReorderSam I=${marked_bam} O=${out_bam_path} R=${ref_fasta}')

run('echo STARTING index')
run('date')
run('samtools index ${out_bam_path}')

"

        echo "$python_cmd"
        set +e
        python -c "$python_cmd"
        export exit_code=$?
        set -e
        echo exit code is $exit_code
        ls

        # create bundle conditional on failure
        if [[ "${debug_dump_flag}" == "always" || ( "${debug_dump_flag}" == "onfail" && $exit_code -ne 0 ) ]]
        then
            echo "Creating debug bundle"
            # tar up the output directory
            touch debug_bundle.tar.gz
            tar cfz debug_bundle.tar.gz --exclude=debug_bundle.tar.gz .
        else
            touch debug_bundle.tar.gz
        fi     
        /opt/src/algutil/monitor_stop.py

        # exit statement must be the last line in the command block 
        exit $exit_code

    }
    output {
        File out_bam = "${out_bam_path}"
        File out_index = "${out_index_path}"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
        File debug_bundle="debug_bundle.tar.gz"
    } runtime {
        docker : "gcr.io/btl-dockers/btl_gatk:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: "${preemptible}"
    }
    parameter_meta {
        picard: "The absolute path to the picard jar to execute."
        in_bam: "The bam file to convert to fastq."
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        out_fq1: "The fastq file containing the first read of each pair."
        out_fq2: "The fastq file containing the second read of each pair"
    }
}

