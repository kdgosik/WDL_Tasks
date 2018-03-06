workflow gatk_tcir {

    call gatk_tcir_task

}



task gatk_tcir_task {
    String gatk_path = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"
    File in_bam
    File in_bam_index
    String sample_name
    File reference_tgz

    String out_bam_fn =  "${sample_name}.tcir.bam"
    String out_bam_index_fn =  "${out_bam_fn}.bai"

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


run('ln -s ${in_bam} in.bam')
run('ln -s ${in_bam_index} in.bam.bai')

run('echo STARTING tar xvf to unpack reference')
run('date')
run('tar xvf ${reference_tgz}')

run('echo STARTING RealignerTargetCreator')
run('date')

run('java -Xmx8G -jar ${gatk_path} -T RealignerTargetCreator -nct 1 -nt 24 -R ref.fasta -I in.bam -o tcir.intervals.list ')


run('echo STARTING IndelRealigner')
run('date')
run('java -Xmx4G -jar ${gatk_path} -T IndelRealigner -nct 1 -nt 1 -R ref.fasta -I in.bam -targetIntervals tcir.intervals.list -o ${out_bam_fn}')



run('echo STARTING index')
run('date')
run('samtools index ${out_bam_fn}')

run('echo DONE')
run('date')
"

        echo "$python_cmd"
        set +e
        python -c "$python_cmd"
        export exit_code=$?
        set -e
        echo exit code is $exit_code
        ls

        # create bundle conditional on failure of the Python section
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
        File out_bam = "${out_bam_fn}"
        File out_bam_index = "${out_bam_index_fn}"
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

    }

}



