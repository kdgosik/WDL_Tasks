task mutation_validator_postprocess {

    #Inputs and constants defined here
    File input_maf
    String pair_id    

    String output_disk_gb
    String boot_disk_gb = "10"
    String ram_gb = "1"
    String cpu_cores = "1"

    command {
python_cmd="
import subprocess
def run(cmd):
    print(cmd)
    subprocess.check_call(cmd,shell=True)

run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')
run('/opt/src/algutil/monitor_start.py')

run('python /opt/src/mutation_validator_postprocess.py \"${input_maf}\"  \"${pair_id}\".validated_annotated.maf')


#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File output_maf="${pair_id}.validated_annotated.maf"
        File dstat_log="dstat.log"
    }

    runtime {
        docker : "docker.io/gsaksena/mutation_validator_postprocess:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: 3
    }


    meta {
        author : "Gordon Saksena"
        email : "gsaksena@broadinstitute.org"
    }

}

workflow mutation_validator_postprocess_workflow {
    call mutation_validator_postprocess
}
