task gatk_acnv_disk_calc {

Float bedSize
Float fastaSize
Float fastaDictSize
Float fastaIdxSize
Float common_snp_list_size
Float tumorBamSize
Float tumorBaiSize
Float normalBamSize
Float normalBaiSize
Float ponSize
Float buffer
Int PreemValue

command <<<

#increase verbosity
set -x

echo -ne "(${bedSize}+${fastaSize}+${fastaDictSize}+${fastaIdxSize}+${common_snp_list_size}+${tumorBamSize}+${tumorBaiSize}+${normalBamSize}+${normalBaiSize}+${ponSize}+${buffer})/1000000000" | perl -ne 'print eval($_);'|grep -Po '^\d+' > acnv_disk_db.dat
>>>

  runtime {
    memory: "0.1 GB"
    disks: "local-disk 1 HDD"
    docker: "broadinstitute/hellbender_pon@sha256:7b9264bfd408ecfc34a3bf75f430ebd667bf454c7bd295eba9c095ca60df47d0"
    preemptible: "${PreemValue}"
  }

output {
  Int acnvDiskGB=read_int("acnv_disk_db.dat")
  }
}


task gatk_acnv_only{

  #pipeline files/settings
  File recapseg_bed
  File ref_fasta
  File ref_fasta_dict
  File ref_fasta_fai
  File common_snp_list
  File tumor_bam
  File tumor_bam_idx
  File normal_bam
  File normal_bam_idx
  File PoN
  String entity_id
  Boolean padTargets

  #Runtime attribute settings
  String Disk_GB
  String PreemValue

  parameter_meta {
    recapseg_bed : "a BED of genomic targets indicating loci for analysis"
    ref_fasta: "fasta file for reference genome "
    ref_fasta_fai: "fasta file index for the reference genome (see http://www.htslib.org/doc/faidx.html)"
    ref_fasta_dict: "fasta file dictionary for the reference genome (see https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary)"
    common_snp_list: "a picard interval-list file indicating common heterozygous sites"
    tumor_bam: "tumor BAM file (see https://samtools.github.io/hts-specs/SAMv1.pdf)"
    tumor_bam_idx : "tumor BAM index file (see samtools index command http://www.htslib.org/doc/samtools.html)"
    normal_bam : "normal BAM file (see https://samtools.github.io/hts-specs/SAMv1.pdf)"
    normal_bam_idx: "normal BAM index file (see samtools index command http://www.htslib.org/doc/samtools.html)"
    PoN: "a panel-of-normals file ; generated for example perhaps by this workflow http://gatkforums.broadinstitute.org/gatk/discussion/comment/31332/"
    entity_id: "a string for the name of the pair under analysis used for naming output files"
    padTargets: "a boolean flag indicating whether PadTargets is invoked on the targets before analysis"
    Disk_GB: "size in GB of the working disks"
    PreemValue: "non-negative interger value for preemptible 0 means not preemptible, otherwise 1,2,... is the max number of pre-emptible tries"
  }

  command 
      <<<

      #increase verbosity
      set -x

      if [ "${padTargets}" == "true" ] ; 
      then
        #perform optional padding
        java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar /root/gatk-protected.jar  PadTargets  \
          --targets ${recapseg_bed} --output padded_bed.bed --padding 250   --help false \
          --version false --verbosity INFO --QUIET false ;
          EFFECTIVE_TARGETS="padded_bed.bed" ;
      else
          EFFECTIVE_TARGETS="${recapseg_bed}"
      fi ;

      #Get coverage
      java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar /root/gatk-protected.jar   CalculateTargetCoverage  \
        --output ${entity_id}.coverage.tsv --groupBy SAMPLE  --transform PCOV --targets $EFFECTIVE_TARGETS  \
        --targetInformationColumns FULL  --keepduplicatereads true --input ${tumor_bam} \
        --reference ${ref_fasta}    --disable_all_read_filters false --interval_set_rule UNION \
        --interval_padding 0 --secondsBetweenProgressUpdates 10.0   \
        --disableSequenceDictionaryValidation false   --createOutputBamIndex true  \
        --help false --version false  --verbosity INFO --QUIET false

      #Normalization
      java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar /root/gatk-protected.jar  NormalizeSomaticReadCounts \
        --input ${entity_id}.coverage.tsv   --targets $EFFECTIVE_TARGETS --panelOfNormals ${PoN}  \
        --factorNormalizedOutput ${entity_id}.fnt.tsv   --tangentNormalized ${entity_id}.tn.tsv  \
        --betaHatsOutput  ${entity_id}.betaHats.tsv  --preTangentNormalized  ${entity_id}.preTN.tsv \
        --help  false   --version false --verbosity INFO --QUIET false


      echo "perform seg"
      java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar /root/gatk-protected.jar PerformSegmentation \
        --targets ${entity_id}.tn.tsv --output ${entity_id}.seg --log2Input true \
        --alpha 0.01 --nperm 10000 --pmethod HYBRID --minWidth 2 --kmax 25 --nmin 200 \
        --eta 0.05 --trim 0.025 --undoSplits NONE --undoPrune 0.05 --undoSD 3 \
        --help false --version false --verbosity INFO --QUIET false


      echo "call segs"
      java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar /root/gatk-protected.jar \
       CallSegments  --targets ${entity_id}.tn.tsv --segments ${entity_id}.seg \
       --output ${entity_id}.called --threshold 2.0  --legacy false  \
       --experimental false --help false --version false --verbosity INFO --QUIET false

      echo "get het cov"
      java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar /root/gatk-protected.jar \
        GetHetCoverage  --reference ${ref_fasta} --normal ${normal_bam} --tumor ${tumor_bam} \
        --snpIntervals ${common_snp_list}  --normalHets ${entity_id}.normal.hets.tsv \
        --tumorHets ${entity_id}.tumor.hets.tsv --pvalueThreshold 0.05  --help false  \
        --version false --verbosity INFO --QUIET false --VALIDATION_STRINGENCY LENIENT

      mkdir -v acnv

      echo "allelic cnv"
      java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar /root/gatk-protected.jar \
        AllelicCNV  --tumorHets ${entity_id}.tumor.hets.tsv \
        --tangentNormalized ${entity_id}.tn.tsv --segments ${entity_id}.called \
        --outputPrefix acnv/${entity_id}  --smallSegmentThreshold 3  \
        --numSamplesCopyRatio 100 --numBurnInCopyRatio 50  \
        --numSamplesAlleleFraction 100 --numBurnInAlleleFraction 50  \
        --intervalThresholdCopyRatio 5.0 --intervalThresholdAlleleFraction 2.0 \
        --help false --version false --verbosity INFO --QUIET false
 
      >>>
  runtime {
    disks: "local-disk ${Disk_GB} HDD"
    memory: "7 GB"
    docker: "broadinstitute/hellbender_pon@sha256:7b9264bfd408ecfc34a3bf75f430ebd667bf454c7bd295eba9c095ca60df47d0"
    preemptible: "${PreemValue}"
  }
  output {
    File gatk_cnv_coverage_file = "${entity_id}.coverage.tsv"
    File gatk_cnv_seg_file = "${entity_id}.seg"
    File gatk_cnv_tn_coverage = "${entity_id}.tn.tsv"
    File gatk_cnv_pre_tn_coverage = "${entity_id}.preTN.tsv"
    File gatk_het_ad_normal= "${entity_id}.normal.hets.tsv"
    File gatk_het_ad_tumor= "${entity_id}.tumor.hets.tsv"
    File gatk_acnv_seg_file = "acnv/${entity_id}-sim-final.acs.seg"
  }
}


workflow gatk_acnv_wkfl {
  #pipeline files/settings
  File recapseg_bed
  File ref_fasta
  File ref_fasta_dict
  File ref_fasta_fai
  File common_snp_list
  File tumor_bam
  File tumor_bam_idx
  File normal_bam
  File normal_bam_idx
  File PoN
  String entity_id
  Boolean padTargets

  #Runtime attribute settings
  String PreemValue
  Float buffer

call gatk_acnv_disk_calc {
  input:
      bedSize=size(recapseg_bed),
      fastaSize=size(ref_fasta),
      fastaDictSize=size(ref_fasta_dict),
      fastaIdxSize=size(ref_fasta_fai),
      common_snp_list_size=size(common_snp_list),
      tumorBamSize=size(tumor_bam),
      tumorBaiSize=size(tumor_bam_idx),
      normalBamSize=size(normal_bam),
      normalBaiSize=size(normal_bam_idx),
      ponSize=size(PoN),
      buffer=buffer
  }


  call gatk_acnv_only {
    input:
        recapseg_bed=recapseg_bed,
        ref_fasta=ref_fasta,
        ref_fasta_dict=ref_fasta_dict,
        ref_fasta_fai=ref_fasta_fai,
        common_snp_list=common_snp_list,
        tumor_bam=tumor_bam,
        tumor_bam_idx=tumor_bam_idx,
        normal_bam=normal_bam,
        normal_bam_idx=normal_bam_idx,
        PoN=PoN,
        Disk_GB=gatk_acnv_disk_calc.acnvDiskGB,
        entity_id=entity_id,
        padTargets=padTargets,
        PreemValue=PreemValue
  }

  output {
    gatk_acnv_only.*
  }


}
