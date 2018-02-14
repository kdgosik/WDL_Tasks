task oncotator_minimal {

    #Inputs and constants defined here
    
    File oncoDBTarBall
    File IN
    String OTHER_FLAGS
    String id

    String output_disk_gb
    String boot_disk_gb = "10"
    String ram_gb = "8"
    String cpu_cores = "2"

    command {
python_cmd="
import subprocess, os
def run(cmd):
    subprocess.check_call(cmd,shell=True)

run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')
#run('/opt/src/algutil/monitor_start.py')

# start task-specific calls
##########################
run('touch start.txt')
run('tar xvf \"${oncoDBTarBall}\" ')
pwd=os.getcwd()
from glob import glob
onco_dir = glob('onco*')
onco_db = onco_dir[0] + ' '
other_flags=\"${OTHER_FLAGS}\"
in1='${IN}'
id1='${id}'
run ('ls -latrh')
cmd='/root/oncotator_venv/bin/oncotator -i MAFLITE --db-dir ' + pwd + '/' + onco_db +  '  --tx-mode=EFFECT -o TCGAMAF ' + other_flags + ' ' +in1 + ' ' +  id1 + '.annotated.maf hg19'
print(cmd)
run(cmd)
run ('ls -latrh')
run('cut -f1-43,66-67,104,108-400 \"${id}.annotated.maf\" > \"${id}.maf\"')
run('tar cvfz \"${id}.annotated.maf.gz\" \"${id}.annotated.maf\"')

#########################
# end task-specific calls
#run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File oncotator_full_maf_gz = "${id}.annotated.maf.gz"
        File oncotator_out = "${id}.maf"
    }

    runtime {
        docker: "broadinstitute/oncotator:1.9.0.0"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: 0
    }


    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }

}

workflow oncotator_minimal_workflow {
    call oncotator_minimal
}
