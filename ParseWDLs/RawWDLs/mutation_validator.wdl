task mutation_validator {

    #Inputs and constants defined here
    String PAIRID
    String maf_type
    File MAF
    File WEXTUMOR
    File WEXNORMAL
    File WEXTUMORBAI
    File WEXNORMALBAI
    File WGSTUMOR
    File WGSNORMAL
    File WGSTUMORBAI
    File WGSNORMALBAI
    File RNATUMOR
    File RNATUMORBAI
    File TARGTUMOR
    File TARGNORMAL
    File TARGTUMORBAI
    File TARGNORMALBAI
    File LPTUMOR
    File LPNORMAL
    File LPTUMORBAI
    File LPNORMALBAI
    File OTUMOR
    File ONORMAL
    File OTUMORBAI
    File ONORMALBAI

    #String output_disk_gb
    String boot_disk_gb = "10"
    String ram_gb = "6"
    String cpu_cores = "1"
    String preemption

    # compute output_disk_gb


    Float S_WEXTUMOR = size(WEXTUMOR)
    Float S_MAF = size(MAF)
    Float S_WEXNORMAL = size(WEXNORMAL)
    Float S_WEXTUMORBAI = size(WEXTUMORBAI)
    Float S_WEXNORMALBAI = size(WEXNORMALBAI)
    Float S_WGSTUMOR = size(WGSTUMOR)
    Float S_WGSNORMAL = size(WGSNORMAL)
    Float S_WGSTUMORBAI = size(WGSTUMORBAI)
    Float S_WGSNORMALBAI = size(WGSNORMALBAI)
    Float S_RNATUMOR = size(RNATUMOR)
    Float S_RNATUMORBAI = size(RNATUMORBAI)
    Float S_TARGTUMOR = size(TARGTUMOR)
    Float S_TARGNORMAL = size(TARGNORMAL)
    Float S_TARGTUMORBAI = size(TARGTUMORBAI)
    Float S_TARGNORMALBAI = size(TARGNORMALBAI)
    Float S_LPTUMOR = size(LPTUMOR)
    Float S_LPNORMAL = size(LPNORMAL)
    Float S_LPTUMORBAI = size(LPTUMORBAI)
    Float S_LPNORMALBAI = size(LPNORMALBAI)
    Float S_OTUMOR = size(OTUMOR)
    Float S_ONORMAL = size(ONORMAL)
    Float S_OTUMORBAI = size(OTUMORBAI)
    Float S_ONORMALBAI = size(ONORMALBAI)

    Float total_file_size_float = S_MAF +S_WEXTUMOR +S_WEXNORMAL+S_WEXTUMORBAI+S_WEXNORMALBAI+S_WGSTUMOR+S_WGSNORMAL+S_WGSTUMORBAI+S_WGSNORMALBAI+S_RNATUMOR+
        S_RNATUMORBAI+S_TARGTUMOR+S_TARGNORMAL+S_TARGTUMORBAI+S_TARGNORMALBAI+S_LPTUMOR+S_LPNORMAL+S_LPTUMORBAI+S_LPNORMALBAI+S_OTUMOR+S_ONORMAL+S_OTUMORBAI+S_ONORMALBAI
    Float disk_size_float = (total_file_size_float / 1000000000) + 10
    String s = disk_size_float
    String string_before_decimal = sub(s, "\\..*", "")
    String output_disk_gb = string_before_decimal
    #Int disk_size_int = string_before_decimal







    command {
cat <<EOF > pyscript.py

import subprocess,os
def run(cmd):
    print('about to run')
    print(cmd)
    print('')
    subprocess.check_call(cmd,shell=True)

run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')
run('/opt/src/algutil/monitor_start.py')

# start task-specific calls
##########################

run('python /opt/src/convert_maf_dialect.py \"${MAF}\" input_variants.maf Broad')

# split MAF into INDELs and SNPs

run('python /opt/src/filter_tsv.py -i input_variants.maf  -f Variant_Type -v \"SNP|DNP|TNP|MNP\" -e .snp.maf')

run('python /opt/src/filter_tsv.py -i input_variants.maf  -f Variant_Type -v \"INS|DEL\" -e .indel.maf')

#run('mkdir -p softlinked ')

#run('ls -latr ')

CWD = os.getcwd() 
PAIRID = '${PAIRID}'

maf_basename = 'input_variants.maf'
MAFSNP=os.path.join(CWD, maf_basename +'.snp.maf')
MAFINDEL=os.path.join(CWD, maf_basename +'.indel.maf')
PREPROCESSED_FILE = os.path.join(CWD, PAIRID +'.pileup_preprocessing.txt')  

VALMAFSNP = os.path.join(CWD, 'snp_mv', PAIRID +'.snp.validated.maf')
VALMAFINDEL = os.path.join(CWD, 'indel_mv', PAIRID +'.indel.validated.maf')

VALMAF = os.path.join(CWD, PAIRID +'.validated.maf')

input_file_table = dict() # braces not allowed when embedded in wdl
input_file_table['WEXT'] = ['${WEXTUMOR}', '${WEXTUMORBAI}']
input_file_table['WEXN'] = ['${WEXNORMAL}', '${WEXNORMALBAI}']
input_file_table['WGST'] = ['${WGSTUMOR}', '${WGSTUMORBAI}']
input_file_table['WGSN'] = ['${WGSNORMAL}', '${WGSNORMALBAI}']
input_file_table['RNAT'] = ['${RNATUMOR}', '${RNATUMORBAI}']
input_file_table['TARGT'] = ['${TARGTUMOR}', '${TARGTUMORBAI}']
input_file_table['TARGN'] = ['${TARGNORMAL}', '${TARGNORMALBAI}']
input_file_table['LPT'] = ['${LPTUMOR}', '${LPTUMORBAI}']
input_file_table['LPN'] = ['${LPNORMAL}', '${LPNORMALBAI}']
input_file_table['OT'] = ['${OTUMOR}', '${OTUMORBAI}']
input_file_table['ON'] = ['${ONORMAL}', '${ONORMALBAI}']
    

calling_file_table = dict()

cwd = os.getcwd()
for filetype in input_file_table:
    bampath = input_file_table[filetype][0]
    baipath = input_file_table[filetype][1]
    if bampath == '' or bampath.endswith('NULL'):
        calling_file_table[filetype] = 'None'
    else: 
        os.symlink(bampath, filetype + '.bam')
        os.symlink(baipath, filetype + '.bam.bai')
        calling_file_table[filetype] = CWD + '/' + filetype + '.bam'

RNATYPE='hg19'
RNATBAM = calling_file_table['RNAT']
if RNATBAM != 'None':
    run('samtools view -H ' + RNATBAM + ' | grep SN:chr1 > RNACHECK.txt')
    with open('RNACHECK.txt', 'r') as f:
        RNACHECK_first_line = f.readline()
    if len(RNACHECK_first_line)>0:
        print('file is aligned with chr in front')
        RNATYPE='hg19-chr'
    
        



cmd1=' --mafsnp '+ MAFSNP + ' --mafindel ' + MAFINDEL + ' --wextumor ' + calling_file_table['WEXT'] + ' --wexnormal ' + calling_file_table['WEXN'] + ' --wgstumor ' +  calling_file_table['WGST'] + ' --wgsnormal ' + calling_file_table['WGSN']
cmd2=' --rnatumor ' + calling_file_table['RNAT'] + ' --targetedtumor ' + calling_file_table['TARGT'] + ' --targetednormal ' + calling_file_table['TARGN'] + ' --lowpasstumor ' + calling_file_table['LPT'] + ' --lowpassnormal ' + calling_file_table['LPN']  
cmd3=' --othertumor ' + calling_file_table['OT'] + ' --othernormal ' + calling_file_table['ON'] + ' --out ' + '${PAIRID}' + ' --rnatype ' + RNATYPE

cmd='python /opt/src/fh_MutationValidatorPreprocess/validation_wrapper_firehose_library_hack.py ' + cmd1 + cmd2 + cmd3
run(cmd)


cmd = 'mkdir snp_mv && cd snp_mv && python /opt/src/algutil/firehose_module_adaptor/run_module.py --module_libdir /opt/src/fh_MutationValidator \
--mutation.validator.preprocessed.file %s \
--maf_file_to_annotate %s \
--discovery_type.wgs_wex_rna_targeted ${maf_type} \
--pair_id  %s.snp \
--print_discovery_counts  true \
--normal_coverage_threshold  0 \
--job.spec.memory 2'%(PREPROCESSED_FILE,MAFSNP,PAIRID)



run(cmd)

cmd='mkdir indel_mv && cd indel_mv && python /opt/src/algutil/firehose_module_adaptor/run_module.py --module_libdir /opt/src/fh_MutationValidator \
--mutation.validator.preprocessed.file %s \
--maf_file_to_annotate %s \
--discovery_type.wgs_wex_rna_targeted ${maf_type} \
--pair_id  %s.indel \
--print_discovery_counts  true \
--normal_coverage_threshold  0 \
--job.spec.memory 2'%(PREPROCESSED_FILE,MAFINDEL,PAIRID) 

run(cmd)


cmd = 'python /opt/src/fh_tsvCatFiles/tsvcat.py %s %s > %s'%(VALMAFSNP, VALMAFINDEL, VALMAF + '.tmp')

run(cmd)

run('python /opt/src/convert_maf_dialect.py %s %s TCGA'%(VALMAF + '.tmp',VALMAF ))


import time
#time.sleep(999999999)


#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
EOF

        cat pyscript.py 
        python pyscript.py

    }

    output {
        File pileup_preprocessing_txt="${PAIRID}.pileup_preprocessing.txt"
        File validated_maf="${PAIRID}.validated.maf"
        File dstat_log="dstat.log"
    }

    runtime {
        docker : "docker.io/gsaksena/mutation_validator:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: "${preemption}"
    }


    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }

}

workflow mutation_validator_workflow {
    call mutation_validator
}
