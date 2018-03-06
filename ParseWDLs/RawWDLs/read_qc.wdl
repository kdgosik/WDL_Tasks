task read_qc {

    #Inputs and constants defined here
    String parent_ssf
    String labset_ssf
    String sample_id
    String genus
    String species
    String strain
    String? specimen_id
    String? gnumber
    File ref_tarball
    File r1_fq_path
    File r2_fq_path
    String bam_path=""

    Float fastq_size = size(r1_fq_path,"GB")
    String output_disk_gb = ceil(fastq_size * 20)
    String ram_gb = ceil(fastq_size + 1)

    #String output_disk_gb
    String boot_disk_gb = "10"
    #String ram_gb = "6"
    String cpu_cores = "1"


    command {
python_cmd="
import subprocess
def run(cmd):
    print (cmd)
    subprocess.check_call(cmd,shell=True)

run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')
run('/opt/src/algutil/monitor_start.py')

# start task-specific calls
##########################
import csv
import os


#unpack ref genome to known location
cwd = os.getcwd()
ref_genome_dir = os.path.join(cwd,'ref_genome')
os.mkdir(ref_genome_dir)
run('tar xvf \"${ref_tarball}\" -C %s'%(ref_genome_dir) )
#get name of fasta file, checking that there is only one
ref_files = os.listdir(ref_genome_dir)
fasta_file = None
for fn in ref_files:
    if fn.endswith('.fasta'):
        if fasta_file is None:
            fasta_file = fn
        else:
            raise Exception('multiple .fasta files found in reference tarball')
fasta_path = os.path.join(ref_genome_dir, fasta_file)
print (fasta_path)

#make input table based on input arguments
#TODO check that either fasta files or bam file is present.
#sample_id	genus	species	strain	specimen_id	gnumber	ref_path	r1_fq_path	r2_fq_path	bam_path
header = ['sample_id','genus','species','strain','specimen_id','gnumber','ref_path','r1_fq_path','r2_fq_path','bam_path']
line=[]
line.append('${sample_id}')
line.append('${genus}')
line.append('${species}')
line.append('${strain}')
line.append('${specimen_id}')
line.append('${gnumber}')
line.append(fasta_path)
line.append('${r1_fq_path}')
line.append('${r2_fq_path}')
line.append('${bam_path}')

fid = open('input.txt','w')
outwriter = csv.writer(fid, dialect='excel-tab', lineterminator='\n')
outwriter.writerow(header)
outwriter.writerow(line)
fid.close()

#config output directory
sample_id = '${sample_id}'
parent_ssf = '${parent_ssf}'
labset_ssf = '${labset_ssf}'

run('mkdir -p /btl/projects/SSF')
run('ln -s `pwd` /btl/projects/SSF/PCR-Free')
outdir = os.path.join(cwd,parent_ssf,labset_ssf)

# run the actual algorithm
run('/opt/src/ssf_read_qc_setup.pl -s %s -l %s  -t  pcr-free -i  input.txt'%(parent_ssf, labset_ssf))

metric_out_fn = '%s_%s_%s.metrics.txt'%(parent_ssf, labset_ssf, sample_id)
metric_pipeline_fn = os.path.join(outdir, 'qc', '%s.metrics.txt'%sample_id)
report_out_fn = '%s_%s_%s.pdf'%(parent_ssf, labset_ssf, sample_id)
report_pipeline_fn = os.path.join(outdir, 'qc', '%s.pdf'%sample_id)

# hard link, to survive export to host
run('ln  %s %s'%(metric_pipeline_fn, metric_out_fn))
run('ln  %s %s'%(report_pipeline_fn, report_out_fn))

fid = open(metric_pipeline_fn)
indict = csv.DictReader(fid, dialect='excel-tab')
line = indict.next()
fid.close()

estimate_library_size_out_fn = '%s_%s_%s.estimate_library_size.txt'%(parent_ssf, labset_ssf, sample_id)
fid = open(estimate_library_size_out_fn,'w')
fid.write(line['estimate_library_size'])
fid.close()

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"
        export exit_code=$?
        echo exit code is $exit_code

        # create bundle conditional on failure
        if [ $exit_code -eq 0 ]
        then
            echo "Passed"
            touch debug_bundle.tar.gz
        else
            echo "Failed"
            # tar up the output directory
            touch debug_bundle.tar.gz
            tar cfz debug_bundle.tar.gz --exclude=debug_bundle.tar.gz  .
        fi     
        exit $exit_code
    }

    output {
        File read_qc_table="${parent_ssf}_${labset_ssf}_${sample_id}.metrics.txt"
        File read_qc_pdf="${parent_ssf}_${labset_ssf}_${sample_id}.pdf"
        File read_qc_debug_bundle="debug_bundle.tar.gz"
        String read_qc_estimate_library_size = read_string("${parent_ssf}_${labset_ssf}_${sample_id}.estimate_library_size.txt")
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    }

    runtime {
        docker: "gcr.io/btl-dockers/read_qc:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: 1
    }


    meta {
        author : "Gordon Saksena"
        email : "gsaksena@broadinstitute.org"
    }

}

workflow read_qc_workflow {
    call read_qc
}
