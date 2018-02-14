task disk_size_calculator{
    String uuid_and_filename
    String uuid=sub(uuid_and_filename, "/.*", "")
    Int preemptible
    String? is_legacy

    command <<<
        python3 /opt/src/disk_size_calculator.py ${is_legacy} ${uuid}
    >>>

    output {
        Int final_disk_size=read_int("file_size.txt")
    }

    runtime {
        docker : "docker.io/broadinstitute/gdc_downloader:1.0"
        preemptible: "${preemptible}"
    }


    meta {
        author : "Ruslana Frazer"
        email : "rfrazer@broadinstitute.org"
    }

}


task gdc_file_downloader {

    #Inputs and constants defined here
    String uuid_and_filename
    Int preemptible
    Int output_disk_gb
    File? gdc_user_token
    
    Int boot_disk_gb = 10
    Int ram_gb = 8
    Int cpu_cores = 8

    command <<<

python_cmd="
import subprocess
def run(cmd):
    subprocess.check_call(cmd,shell=True)

run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')
run('/opt/src/algutil/monitor_start.py')

# start task-specific calls
##########################

run('python3 /opt/src/gdc_downloader.py ${"-t " + gdc_user_token} ${uuid_and_filename}')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd" 
        python -c "$python_cmd"

    >>>

    output {
        File file="${uuid_and_filename}"
	File dstat_log="dstat.log"
    }

    runtime {
        docker : "docker.io/broadinstitute/gdc_downloader:1.0"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: "${preemptible}"
    }


    meta {
        author : "Ruslana Frazer"
        email : "rfrazer@broadinstitute.org"
    }

}

workflow gdc_file_downloader_workflow {
    String uuid_and_filename
    Int preemptible
    Boolean is_legacy

    if (is_legacy){
        String flag_for_legacy_archive="-l"
    }
    
    call disk_size_calculator {input: uuid_and_filename=uuid_and_filename,preemptible=preemptible, is_legacy=flag_for_legacy_archive}

    call gdc_file_downloader {input: uuid_and_filename=uuid_and_filename,preemptible=preemptible, output_disk_gb=disk_size_calculator.final_disk_size}

    output {
    	   gdc_file_downloader.file
    }

}
