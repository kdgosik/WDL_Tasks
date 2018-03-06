workflow gatk_haplotypecaller {

    File ? bqsr_bam
    File ? indelrealigner_bam
    File ? uncleaned_bam
    File ? bqsr_bam_index
    File ? indelrealigner_bam_index
    File ? uncleaned_bam_index

    File in_bam = select_first([bqsr_bam, indelrealigner_bam, uncleaned_bam])
    File in_bam_index = select_first([bqsr_bam_index, indelrealigner_bam_index, uncleaned_bam_index])

    call gatk_haplotypecaller_task {
        input:
        in_bam = in_bam,
        in_bam_index = in_bam_index
    }

}



task gatk_haplotypecaller_task {
    String gatk_path = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"
    File in_bam
    File in_bam_index
    String sample_name

    #Float interval_size
    File ? bqsr_file
    String ? ploidy
    String ? erc
    String ? extra_hc_params



    File reference_tgz
    Array[String] known_sites


    String out_gvcf_fn = "${sample_name}.gvcf"

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

# Add intervals back in when actually scattering
#
#run('''\
#python /opt/src/intervals_creator.py \
#    -r ref.fasta \
#    -i $interval_size \
#    > intervals.list
#''')
#			--intervals intervals.list \
#			--interval_padding 100 \

run('''\

        java -Xmx8G -jar ${gatk_path} \
            -T HaplotypeCaller \
            -nt 1 \
            -R ref.fasta \
            --input_file ${in_bam} \
            ${"-BQSR " + bqsr_file} \
            -ERC ${default="GVCF" erc} \
            -ploidy ${default="2" ploidy} \
            -o ${out_gvcf_fn} \
            -variant_index_type LINEAR \
            -variant_index_parameter 128000 \
            ${default="\n" extra_hc_params}
''')


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
        File out_gvcf = "${out_gvcf_fn}"
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





