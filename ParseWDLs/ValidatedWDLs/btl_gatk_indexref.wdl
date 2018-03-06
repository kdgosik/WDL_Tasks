
workflow gatk_indexref {

    # Check's if index files exist(using .dict file as marker). Will always output path to localized reference
    # With assumption that it either exists or will be created by IndexReference.
    call IndexReference   
}


task IndexReference {
    File ref_fasta
    #TODO try to autocalculate ref_name from ref_fasta input
    String ref_name

    String output_disk_gb 
    String boot_disk_gb = "10"
    String ram_gb = "3"
    String cpu_cores = "1"
    String preemptible = "0"
    String debug_dump_flag

    String picard = "/cil/shed/apps/external/picard/current/bin/picard.jar"

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

run('ln ${ref_fasta} ref.fasta')
run('bwa index ref.fasta')
run('samtools faidx ref.fasta')
run('java -jar ${picard} CreateSequenceDictionary REFERENCE=ref.fasta O=ref.dict')

run('tar cvf ${ref_name}.tar ref.fasta ref.dict ref.fasta.amb ref.fasta.ann ref.fasta.bwt ref.fasta.fai ref.fasta.pac ref.fasta.sa')
run('gzip -c --best ${ref_name}.tar > ${ref_name}.tgz')
"
    echo "$python_cmd"
    python -c "$python_cmd"
    export exit_code=$?
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
        File reference_tgz = "${ref_name}.tgz"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
        File debug_bundle="debug_bundle.tar.gz"
    } runtime {
        task_name: "IndexReference"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: "${preemptible}"
        
    }
    parameter_meta {
        ref_fasta: "The name of the reference file (without the path)."
        picard: "The path to the picard executable jar file."
        output_dir: "The root directory for where all sample directories and ref index files will be deposited."
        out_file: "The seq dict filename must be generated and is used to construct seq_dict"
        seq_dict: "The sequence dictionary created by Picard CreateSequenceDictionary"
        old_ref: "The full path to the original reference file location."
        ref: "The absolute path of the reference file to be used by the workflow."
    }
}


