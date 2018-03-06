task gatk_acnv_only{

  File jar_fn
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

  command {

    python /opt/get_targets.py ${PoN}

    java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar ${jar_fn} CalculateTargetCoverage --output ${entity_id}.coverage.tsv --groupBy SAMPLE --transform PCOV --targets padded_bed.bed --targetInformationColumns FULL --input ${tumor_bam} --reference ${ref_fasta} --interval_set_rule UNION --interval_padding 0 --secondsBetweenProgressUpdates 10.0 --createOutputBamIndex true --help false --version false --verbosity INFO --QUIET false --disableToolDefaultReadFilters false --disableSequenceDictionaryValidation true --disableReadFilter NotDuplicateReadFilter

    java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar ${jar_fn} NormalizeSomaticReadCounts  --input ${entity_id}.coverage.tsv --targets padded_bed.bed --panelOfNormals ${PoN} --factorNormalizedOutput ${entity_id}.fnt.tsv --tangentNormalized ${entity_id}.tn.tsv --betaHatsOutput  ${entity_id}.betaHats.tsv --preTangentNormalized  ${entity_id}.preTN.tsv  --help false --version false --verbosity INFO --QUIET false

    java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar ${jar_fn} PerformSegmentation  --tangentNormalized ${entity_id}.tn.tsv --output ${entity_id}.seg --log2Input true  --alpha 0.01 --nperm 10000 --pmethod HYBRID --minWidth 2 --kmax 25 --nmin 200 --eta 0.05 --trim 0.025 --undoSplits NONE --undoPrune 0.05 --undoSD 3 --help false --version false --verbosity INFO --QUIET false

    mkdir plotting

    java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar ${jar_fn} CallSegments  --tangentNormalized ${entity_id}.tn.tsv --segments ${entity_id}.seg --output ${entity_id}.called --legacy false --help false --version false --verbosity INFO --QUIET false

    java -Djava.library.path=/usr/lib/jni/ -jar ${jar_fn} PlotSegmentedCopyRatio --tangentNormalized ${entity_id}.tn.tsv --preTangentNormalized  ${entity_id}.preTN.tsv --sequenceDictionaryFile ${ref_fasta_dict} --outputPrefix ${entity_id} --segments  ${entity_id}.seg --output plotting/ --log2Input true --help false --version false --verbosity INFO --QUIET false

    java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar ${jar_fn} GetHetCoverage  --reference ${ref_fasta} --normal ${normal_bam} --tumor ${tumor_bam} --snpIntervals ${common_snp_list}  --normalHets ${entity_id}.normal.hets.tsv --tumorHets ${entity_id}.tumor.hets.tsv --pvalueThreshold 0.05  --help false --version false --verbosity INFO --QUIET false --VALIDATION_STRINGENCY LENIENT

    mkdir acnv

    java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar ${jar_fn} AllelicCNV --tumorHets ${entity_id}.tumor.hets.tsv --tangentNormalized ${entity_id}.tn.tsv --segments ${entity_id}.called --outputPrefix acnv/${entity_id}  --smallSegmentThreshold 3 --numSamplesCopyRatio 100 --numBurnInCopyRatio 50 --numSamplesAlleleFraction 100 --numBurnInAlleleFraction 50 --intervalThresholdCopyRatio 5.0 --intervalThresholdAlleleFraction 2.0  --help false --version false --verbosity INFO --QUIET false

    mkdir plotting_acnv

    java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar ${jar_fn} ConvertACNVResults --tumorHets  ${entity_id}.tumor.hets.tsv --segments acnv/${entity_id}-sim-final.seg --outputDir "acnv" --tangentNormalized ${entity_id}.tn.tsv || true #hide failure on titan making 

    java -Xmx7g -Djava.library.path=/usr/lib/jni/ -jar ${jar_fn} PlotACNVResults --hets ${entity_id}.tumor.hets.tsv --tangentNormalized ${entity_id}.tn.tsv --segments acnv/${entity_id}-sim-final.seg --sequenceDictionaryFile ${ref_fasta_dict} --outputPrefix ${entity_id} --output plotting_acnv

  }

  runtime {
    disks: "local-disk 100 HDD"
    memory: "7 GB"
    docker: "dlivitzbroad/cnvvdev"
    preemptible: 3
  }

  output {
    File gatk_cnv_coverage_file = "${entity_id}.coverage.tsv"
    File gatk_cnv_seg_file = "${entity_id}.seg"
    File gatk_cnv_tn_coverage = "${entity_id}.tn.tsv"
    File gatk_cnv_pre_tn_coverage = "${entity_id}.preTN.tsv"
    File gatk_het_ad_normal = "${entity_id}.normal.hets.tsv"
    File gatk_het_ad_tumor = "${entity_id}.tumor.hets.tsv"
    Array[File] gatk_cnv_all_plots = glob("plotting/*.png")
    File gatk_acnv_plot = "plotting_acnv/${entity_id}_ACNV.png"
    File gatk_acnv_seg_file = "acnv/${entity_id}-sim-final.acs.seg"
    File gatk_acnv_raw_seg_file="acnv/${entity_id}-sim-final.seg"
    File gatk_af_param_fie="acnv/${entity_id}-sim-final.af.param"
  }
}

workflow gatk_acnv_wkfl {
  call gatk_acnv_only {} 
}