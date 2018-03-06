task gatk_cov_pull {
  File recapseg_bed
  File ref_fasta
  File ref_fasta_dict
  File ref_fasta_fai
  File normal_bam
  File normal_bam_idx
  String Disk_GB  
  String entity_id
    command <<<

      #increase verbosity
      set -x

     
        #perform optional padding
        java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar /root/gatk-protected.jar  PadTargets  \
          --targets ${recapseg_bed} --output padded_bed.bed --padding 250   --help false \
          --version false --verbosity INFO --QUIET false 
     

      #Get coverage
      java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar /root/gatk-protected.jar   CalculateTargetCoverage  \
        --output ${entity_id}.coverage.tsv --groupBy SAMPLE  --transform PCOV --targets padded_bed.bed  \
        --targetInformationColumns FULL  --keepduplicatereads true --input ${normal_bam} \
        --reference ${ref_fasta}    --disable_all_read_filters false --interval_set_rule UNION \
        --interval_padding 0 --secondsBetweenProgressUpdates 10.0   \
        --disableSequenceDictionaryValidation false   --createOutputBamIndex true  \
        --help false --version false  --verbosity INFO --QUIET false
        
      >>>
      
        output {
    File gatk_cnv_coverage_file = ("${entity_id}.coverage.tsv")
   
  }
  runtime {
    disks: "local-disk ${Disk_GB} HDD"
    memory: "7 GB"
    docker: "broadinstitute/hellbender_pon"
    preemptible: 3
  }

}
workflow gatk_cov_pullwkflow {
  File recapseg_bed
  File normal_bam
  File normal_bam_idx
  File ref_fasta
  File ref_fasta_dict
  File ref_fasta_fai
  String Disk_GB  
  String entity_id
  call gatk_cov_pull {
    input:
        recapseg_bed=recapseg_bed,
        normal_bam=normal_bam,
        ref_fasta=ref_fasta,
        ref_fasta_dict=ref_fasta_dict,
        ref_fasta_fai=ref_fasta_fai,
        normal_bam_idx=normal_bam_idx,
        Disk_GB=Disk_GB,
        entity_id=entity_id
  }  
}
