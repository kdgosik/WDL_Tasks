task samtools_filter {

    #Inputs and constants defined here
    File input_bam
    File input_bai
    String samtools_view_args

    String output_disk_gb
    String boot_disk_gb = "10"
    String ram_gb = "8"
    String cpu_cores = "2"

    command {
python_cmd="
import subprocess
def run(cmd):
    print(cmd)
    subprocess.check_call(cmd,shell=True)


run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')
run('/opt/src/algutil/monitor_start.py')

# start task-specific calls
##########################
run ('ln -s ${input_bam} input.bam')
run ('ln -s ${input_bai} input.bam.bai')

run ('samtools view -b ${samtools_view_args} input.bam > output_unsorted.bam')
run ('samtools sort output_unsorted.bam output')
run ('samtools index output.bam output.bam.bai')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File output_bam="output.bam"
        File output_bai="output.bam.bai"
    }

    runtime {
        docker : "docker.io/gsaksena/samtools_filter:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: 0
    }


    meta {
        author : "Your Name"
        email : "Your.Email@Address.com"
    }

}

workflow samtools_filter_workflow {
    call samtools_filter
}
