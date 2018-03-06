workflow TumorOnly  {

  File intervals 
  File ref_fasta
  File ref_fasta_idx
  File ref_dict
  File tumor_bam
  File tumor_bam_idx
  String tumor_sample_name
  File GATK_protected_jar
  File GATK_protected_jar_with_npe_fix

  ### Mutect 2 section ###
  File unmatched_normal_bam
  File unmatched_normal_bam_idx
  File cosmic
  File cosmicIndex
  File dbsnp
  File dbsnpIndex
  File mutect_pon
  File exac
  File exac_idx
  File onekg
  Int medium_disk
  Int large_disk
  String GATK_jar
  String picard_jar
  File onco_ds_tar_gz
  Array[String] artifact_modes

  call SplitIntervals {
    input:
      intervals = intervals,
      disk_size = medium_disk,
      picard_jar = picard_jar
  }

  scatter (subIntervals in SplitIntervals.scattered_intervals) {
    call M2 {
      input: 
        ref_fasta = ref_fasta,
        ref_fasta_idx = ref_fasta_idx,
        ref_dict = ref_dict,
        tumor = tumor_bam,        
        tumorIndex = tumor_bam_idx,
        normal = unmatched_normal_bam,
        normalIndex = unmatched_normal_bam_idx,
        vcfBasename = tumor_sample_name,
        intervals = subIntervals,
        disk_size = large_disk,
        dbsnp = dbsnp,
        mutect_pon = mutect_pon,
        cosmicIndex = cosmicIndex,
        dbsnpIndex = dbsnpIndex,
        cosmic = cosmic,
        GATK_jar = GATK_jar
    } 

    call AnnotateVariants {
      input:
        ref_fasta = ref_fasta,
        ref_fasta_idx = ref_fasta_idx,
        ref_dict = ref_dict,
        m2_vcf = M2.output_vcf,
        m2_vcf_idx = M2.output_vcf_index,
        output_basename = tumor_sample_name,
        intervals = subIntervals,
        disk_size = medium_disk,
        onekg = onekg,
        exac = exac,
        exac_idx = exac_idx,
        GATK_jar = GATK_jar
    }

    call VariantFiltration as ExacFiltration {
      input:
        ref_fasta = ref_fasta,
        ref_fasta_idx = ref_fasta_idx,
        ref_dict = ref_dict,
        input_vcf = AnnotateVariants.output_vcf,
        input_vcf_index = AnnotateVariants.output_vcf_index,
        output_basename = tumor_sample_name,
        jexl_expression = "exac.AC > 0",
        filter_name = "in_exac",
        disk_size = medium_disk,
        GATK_jar = GATK_jar
    }

    call VariantFiltration as OneKGFilteration {
      input:
        ref_fasta = ref_fasta,
        ref_fasta_idx = ref_fasta_idx,
        ref_dict = ref_dict,
        input_vcf = ExacFiltration.output_vcf,
        input_vcf_index = ExacFiltration.output_vcf_index,
        output_basename = tumor_sample_name,
        jexl_expression = "onekg.AC > 0",
        filter_name = "in_1kg",
        disk_size = medium_disk,
        GATK_jar = GATK_jar
    }
  }

  call GatherVCFs {
    input:
      input_vcfs = OneKGFilteration.output_vcf,
      output_vcf_name = tumor_sample_name,
      disk_size = medium_disk,
      picard_jar = picard_jar
  }

  call CollectSequencingArtifactMetricsPicard {
    input:
      picard_jar = picard_jar,
      bam_file = tumor_bam,
      ref_fasta = ref_fasta,
      ref_fasta_idx = ref_fasta_idx,
      output_prepend = tumor_sample_name,
  }

  call ConvertPicardSequencingArtifactMetricsToGATK {
    input:
      entity_id=tumor_sample_name,
      pre_adapter_detail_metrics_from_picard=CollectSequencingArtifactMetricsPicard.pre_adapter_detail_metrics
  }

  call FilterByOrientationBias {
    input:
      gatk4_jar = GATK_protected_jar_with_npe_fix,
      output_prepend = tumor_sample_name,
      m2_vcf = GatherVCFs.output_vcf,
      pre_adapter_detail_metrics = ConvertPicardSequencingArtifactMetricsToGATK.gatk_pre_adapter_detail_metrics,
      artifact_modes = artifact_modes
  }

  call CollectMetrics {
    input:
      input_vcf = FilterByOrientationBias.orientation_bias_vcf,
      input_vcf_index = FilterByOrientationBias.orientation_bias_vcf_index,
      cosmic = cosmic,
      cosmicIndex = cosmicIndex,
      intervals = intervals,
      ref_dict = ref_dict,
      output_basename = tumor_sample_name,
      disk_size = medium_disk,
      picard_jar = picard_jar
  }

  ### Oncotator on Mutect ouput ###
  call oncotate_m2_to_maf {
    input:
      m2_vcf = FilterByOrientationBias.orientation_bias_vcf,
      entity_id = tumor_sample_name,
      onco_ds_tar_gz = onco_ds_tar_gz
  }

  ### GATK-CNV Section ###
  # Workflow input files
    File cnv_pon
    Int padding
    File targets

    # Input parameters of the PerformSegmentation tool
    Float seg_param_alpha
    Int seg_param_nperm
    String seg_param_pmethod
    Int seg_param_minWidth
    Int seg_param_kmax
    Int seg_param_nmin
    Float seg_param_eta
    Float seg_param_trim
    String seg_param_undoSplits
    Float seg_param_undoPrune
    Int seg_param_undoSD

    # CalculateTargetCoverage options
    Boolean disable_all_read_filters
    Boolean keep_duplicate_reads
    Boolean disable_sequence_dictionary_validation
    String transform
    String grouping

    # Workflow output directories and other options
    String plots_dir
    String conversion_dir
    Boolean enable_gc_correction
    Boolean isWGS
    Int wgsBinSize

    # Java maximum memory options
    Int calculate_target_coverage_memory
    Int normalize_somatic_read_count_memory
    Int whole_genome_coverage_memory

  call PadTargets {
    input:
        target_file=targets,
        gatk_jar=GATK_protected_jar,
        isWGS=isWGS,
        mem=1,
        padding=padding
  }

  call CalculateTargetCoverage as TumorCalculateTargetCoverage {
    input:
        entity_id=tumor_sample_name,
        padded_target_file=PadTargets.padded_target_file,
        input_bam=tumor_bam,
        bam_idx=tumor_bam_idx,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_idx,
        ref_fasta_dict=ref_dict,
        gatk_jar=GATK_protected_jar,
        disable_sequence_dictionary_validation=disable_sequence_dictionary_validation,
        disable_all_read_filters=disable_all_read_filters,
        keep_duplicate_reads=keep_duplicate_reads,
        transform=transform,
        grouping=grouping,
        isWGS=isWGS,
        mem=calculate_target_coverage_memory
  }  

  call WholeGenomeCoverage as TumorWholeGenomeCoverage {
    input:
        entity_id=tumor_sample_name,
        target_file=PadTargets.padded_target_file,
        input_bam=tumor_bam,
        bam_idx=tumor_bam_idx,
        coverage_file=TumorCalculateTargetCoverage.gatk_coverage_file,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_idx,
        ref_fasta_dict=ref_dict,
        gatk_jar=GATK_protected_jar,
        isWGS=isWGS,
        wgsBinSize=wgsBinSize,
        mem=whole_genome_coverage_memory
  }

  call AnnotateTargets as TumorAnnotateTargets {
    input:
        entity_id=tumor_sample_name,
        gatk_jar=GATK_protected_jar,
        target_file=TumorWholeGenomeCoverage.gatk_target_file,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_idx,
        ref_fasta_dict=ref_dict,
        enable_gc_correction=enable_gc_correction,
        mem=4
  }

  call CorrectGCBias as TumorCorrectGCBias {
    input:
        entity_id=tumor_sample_name,
        gatk_jar=GATK_protected_jar,
        coverage_file=TumorWholeGenomeCoverage.gatk_coverage_file,
        annotated_targets=TumorAnnotateTargets.annotated_targets,
        enable_gc_correction=enable_gc_correction,
        mem=4
  }

  call NormalizeSomaticReadCounts as TumorNormalizeSomaticReadCounts {
    input:
        entity_id=tumor_sample_name,
        coverage_file=TumorCorrectGCBias.gatk_cnv_coverage_file_gcbias,
        padded_target_file=TumorWholeGenomeCoverage.gatk_target_file,
        pon=cnv_pon,
        gatk_jar=GATK_protected_jar,
        mem=normalize_somatic_read_count_memory
  }

  call PerformSegmentation as TumorPerformSeg {
    input:
        entity_id=tumor_sample_name,
        gatk_jar=GATK_protected_jar,
        tn_file=TumorNormalizeSomaticReadCounts.tn_file,
        seg_param_alpha=seg_param_alpha,
        seg_param_nperm=seg_param_nperm,
        seg_param_pmethod=seg_param_pmethod,
        seg_param_minWidth=seg_param_minWidth,
        seg_param_kmax=seg_param_kmax,
        seg_param_nmin=seg_param_nmin,
        seg_param_eta=seg_param_eta,
        seg_param_trim=seg_param_trim,
        seg_param_undoSplits=seg_param_undoSplits,
        seg_param_undoPrune=seg_param_undoPrune,
        seg_param_undoSD=seg_param_undoSD,
        mem=2
  }

  call Caller as TumorCaller {
    input:
        entity_id=tumor_sample_name,
        gatk_jar=GATK_protected_jar,
        tn_file=TumorNormalizeSomaticReadCounts.tn_file,
        seg_file=TumorPerformSeg.seg_file,
        mem=2
  }

#Removed trailing / in output dir to fix call caching
  call PlotSegmentedCopyRatio {
    input:
        entity_id=tumor_sample_name,
        gatk_jar=GATK_protected_jar,
        tn_file=TumorNormalizeSomaticReadCounts.tn_file,
        pre_tn_file=TumorNormalizeSomaticReadCounts.pre_tn_file,
        called_file=TumorCaller.called_file,
        ref_fasta_dict=ref_dict,
        output_dir="${plots_dir}CopyRatio_Plots/${tumor_sample_name}",
        mem=4
  }

  ### Oncotator on CNV output ###

  call oncotate_cnv {
    input:
      entity_id=tumor_sample_name,
      target_seg_gt_file=TumorCaller.called_file
  }

}

#Must be hardcoded outputs due to difficulties with glob.
task SplitIntervals {
  Int disk_size
  String picard_jar
  File intervals

  command {
    set -e
    mkdir scattered

    java -jar ${picard_jar} IntervalListTools \
      I=${intervals} \
      O=scattered \
      SCATTER_COUNT=30
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
    memory: "4 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Array[File] scattered_intervals = ["scattered/temp_0001_of_30/scattered.interval_list",
    "scattered/temp_0002_of_30/scattered.interval_list",
    "scattered/temp_0003_of_30/scattered.interval_list",
    "scattered/temp_0004_of_30/scattered.interval_list",
    "scattered/temp_0005_of_30/scattered.interval_list",
    "scattered/temp_0006_of_30/scattered.interval_list",
    "scattered/temp_0007_of_30/scattered.interval_list",
    "scattered/temp_0008_of_30/scattered.interval_list",
    "scattered/temp_0009_of_30/scattered.interval_list",
    "scattered/temp_0010_of_30/scattered.interval_list",
    "scattered/temp_0011_of_30/scattered.interval_list",
    "scattered/temp_0012_of_30/scattered.interval_list",
    "scattered/temp_0013_of_30/scattered.interval_list",
    "scattered/temp_0014_of_30/scattered.interval_list",
    "scattered/temp_0015_of_30/scattered.interval_list",
    "scattered/temp_0016_of_30/scattered.interval_list",
    "scattered/temp_0017_of_30/scattered.interval_list",
    "scattered/temp_0018_of_30/scattered.interval_list",
    "scattered/temp_0019_of_30/scattered.interval_list",
    "scattered/temp_0020_of_30/scattered.interval_list",
    "scattered/temp_0021_of_30/scattered.interval_list",
    "scattered/temp_0022_of_30/scattered.interval_list",
    "scattered/temp_0023_of_30/scattered.interval_list",
    "scattered/temp_0024_of_30/scattered.interval_list",
    "scattered/temp_0025_of_30/scattered.interval_list",
    "scattered/temp_0026_of_30/scattered.interval_list",
    "scattered/temp_0027_of_30/scattered.interval_list",
    "scattered/temp_0028_of_30/scattered.interval_list",
    "scattered/temp_0029_of_30/scattered.interval_list",
    "scattered/temp_0030_of_30/scattered.interval_list"
    ]
  }
}

task M2 {
  File tumor
  File tumorIndex
  File normal
  File normalIndex
  String vcfBasename
  File intervals
  File ref_fasta 
  File ref_fasta_idx 
  File ref_dict 
  File mutect_pon
  File dbsnp 
  File dbsnpIndex
  File cosmic 
  File cosmicIndex
  Int disk_size
  String GATK_jar

  command {
  	java -Xmx4g -jar ${GATK_jar} -T MuTect2 \
      -R ${ref_fasta} \
      -I:tumor ${tumor} \
      -I:normal ${normal} \
      --dbsnp ${dbsnp} \
      --cosmic ${cosmic} \
      --normal_panel ${mutect_pon} \
      -L ${intervals} \
      -o ${vcfBasename}.vcf.gz
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
    memory: "5 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = "${vcfBasename}.vcf.gz"
    File output_vcf_index = "${vcfBasename}.vcf.gz.tbi"
  }
}

# Adopted from https://github.com/broadinstitute/dsde-pipelines/blob/master/genomes_in_the_cloud/single_sample/PairedSingleSampleWf.wdl
# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs

task GatherVCFs {
  Array[File] input_vcfs
  String output_vcf_name
  Int disk_size
  String picard_jar

  command {
    java -Xmx2g -jar ${picard_jar} \
      MergeVcfs \
      INPUT=${sep=' INPUT=' input_vcfs} \
      OUTPUT=${output_vcf_name}.vcf
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = "${output_vcf_name}.vcf"
    File output_vcf_index = "${output_vcf_name}.vcf.idx"
  }
}

# Annotate the Mutect2 vcf with allele counts in exac and thousand genomes vcf
task AnnotateVariants {
  File ref_fasta
  File ref_fasta_idx
  File ref_dict
  File m2_vcf
  File m2_vcf_idx
  String output_basename
  File intervals
  Int disk_size
  String GATK_jar

  File exac 
  File exac_idx 
  File onekg 

  command {
    java -jar ${GATK_jar} -T VariantAnnotator \
      -R ${ref_fasta} \
      -V ${m2_vcf} \
      --resource:exac ${exac} \
      --resource:onekg ${onekg} \
      --expression exac.AC \
      --expression onekg.AC \
      -o ${output_basename}.annotated.vcf.gz \
      --resourceAlleleConcordance \
      -L ${intervals}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
    memory: "5 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = "${output_basename}.annotated.vcf.gz"
    File output_vcf_index = "${output_basename}.annotated.vcf.gz.tbi"
  }
}

task VariantFiltration {
    File ref_fasta
    File ref_fasta_idx
    File ref_dict
    File input_vcf
    File input_vcf_index
    String output_basename
    String jexl_expression
    String filter_name
    Int disk_size
    String GATK_jar

    command {
      java -jar ${GATK_jar} \
        -T VariantFiltration \
        -R ${ref_fasta} \
        -V ${input_vcf} \
        -o ${output_basename}.filtered.vcf.gz \
        --filterExpression "${jexl_expression}" \
        --filterName "${filter_name}"
    }
    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      memory: "5 GB"
      disks: "local-disk " + disk_size + " HDD"
    }
    output {
        File output_vcf = "${output_basename}.filtered.vcf.gz"
        File output_vcf_index = "${output_basename}.filtered.vcf.gz.tbi"
    }
}

task CollectSequencingArtifactMetricsPicard {
    String picard_jar
    File bam_file
    String output_prepend
    File ref_fasta
    File ref_fasta_idx

    command {
      java -jar ${picard_jar} CollectSequencingArtifactMetrics \
        I=${bam_file} O=${output_prepend}  R=${ref_fasta} VALIDATION_STRINGENCY=SILENT
    }
    
    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      memory: "16 GB"
      disks: "local-disk " + 500 + " HDD"
    }
    
    output {
        File pre_adapter_detail_metrics = "${output_prepend}.pre_adapter_detail_metrics"
        File pre_adapter_summary_metrics = "${output_prepend}.pre_adapter_summary_metrics"
        File bait_bias_detail_metrics = "${output_prepend}.bait_bias_detail_metrics"
        File bait_bias_summary_metrics = "${output_prepend}.bait_bias_summary_metrics"
    }
}

task ConvertPicardSequencingArtifactMetricsToGATK {
    # picard.analysis.artifacts.SequencingArtifactMetrics$PreAdapterDetailMetrics

    String entity_id
    File pre_adapter_detail_metrics_from_picard

    command {
        set -x
        sed -r "s/picard\.analysis\.artifacts\.SequencingArtifactMetrics\\\$PreAdapterDetailMetrics/org\.broadinstitute\.hellbender\.tools\.picard\.analysis\.artifacts\.SequencingArtifactMetrics\$PreAdapterDetailMetrics/g" ${pre_adapter_detail_metrics_from_picard} > ${entity_id}.gatk.pre_adapter_detail_metrics
    }

    output {
        File gatk_pre_adapter_detail_metrics = "${entity_id}.gatk.pre_adapter_detail_metrics"
    }
    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      memory: "4 GB"
      disks: "local-disk " + 500 + " HDD"
    }
}

task FilterByOrientationBias {
    File gatk4_jar
    String output_prepend
    File m2_vcf
    File pre_adapter_detail_metrics
    Array[String] artifact_modes

    command {
      java -jar ${gatk4_jar} FilterByOrientationBias -A ${sep=" -A " artifact_modes} \
        -V ${m2_vcf} -P ${pre_adapter_detail_metrics} --output ${output_prepend}.ob_filtered.vcf
    }

    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      memory: "5 GB"
      disks: "local-disk " + 500 + " HDD"
    }

    output {
        File orientation_bias_vcf = "${output_prepend}.ob_filtered.vcf"
        File orientation_bias_vcf_index = "${output_prepend}.ob_filtered.vcf.idx"
        File orientation_bias_vcf_summary = "${output_prepend}.ob_filtered.vcf.summary"
    }
}

task CollectMetrics {
  File input_vcf
  File input_vcf_index
  File cosmic
  File cosmicIndex
  File intervals
  File ref_dict
  String output_basename
  Int disk_size
  String picard_jar

  command {
    java -jar ${picard_jar} \
      CollectVariantCallingMetrics \
      INPUT=${input_vcf} \
      OUTPUT=${output_basename} \
      DBSNP=${cosmic} \
      TARGET_INTERVALS=${intervals} \
      SEQUENCE_DICTIONARY=${ref_dict}
  }
  runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      memory: "5 GB"
      disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File detail_metrics="${output_basename}.variant_calling_detail_metrics"
    File summary_metrics="${output_basename}.variant_calling_summary_metrics"
  }
}

# Pad the target file. This was found to help sensitivity and specificity. This step should only be altered
# by advanced users. Note that by changing this, you need to have a PoN that also reflects the change.
#
# Note that when isWGS is true, this task is still called by the workflow.
# In that case, an empty target file is created and passed to the CalculateTargetCoverage 
# task to satisfy input and output requirements.
# Regardless of the value of isWGS, the output of PadTargets is passed to WholeGenomeCoverage, 
# which then outputs the target file accordingly (see below)
task PadTargets {
    File target_file
    Int padding
    File gatk_jar
    Boolean isWGS
    Int mem

    command {
        if [ ${isWGS} = false ]; \
          then java -Xmx${mem}g -jar ${gatk_jar} PadTargets --targets ${target_file} --output targets.padded.tsv \
                --padding ${padding} --help false --version false --verbosity INFO --QUIET false; \
          else touch targets.padded.tsv; \
        fi
    }

    output {
        File padded_target_file = "targets.padded.tsv"
    }

    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      zones: "us-central1-a us-central1-b"
      disks: "local-disk 200 SSD"
      memory: "6G"
    }
}

# Calculate the target proportional coverage
task CalculateTargetCoverage {
    String entity_id
    File padded_target_file
    String transform
    String grouping
    Boolean keep_duplicate_reads
    Boolean disable_all_read_filters
    File input_bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File gatk_jar
    Boolean disable_sequence_dictionary_validation
    Boolean isWGS
    Int mem

    # Note that when isWGS is true, this task is still called by the workflow.
    # In that case, an empty coverage file is created and passed to the WholeGenomeCoverage
    # task to satisfy input and output requirements.
    command <<<
        if [ ${isWGS} = false ]
          then
              java -Xmx${mem}g -jar ${gatk_jar} CalculateTargetCoverage --output ${entity_id}.coverage.tsv \
                --groupBy ${grouping} --transform ${transform} --targets ${padded_target_file} --targetInformationColumns FULL \
                --input ${input_bam} --reference ${ref_fasta} --disableAllReadFilters ${disable_all_read_filters} \
                $(if [ ${keep_duplicate_reads} = true ]; then echo " --disableReadFilter NotDuplicateReadFilter "; else echo ""; fi) \
                --interval_set_rule UNION --interval_padding 0 \
                --secondsBetweenProgressUpdates 10.0 --disableSequenceDictionaryValidation ${disable_sequence_dictionary_validation} \
                --createOutputBamIndex true --help false --version false --verbosity INFO --QUIET false
          else
              touch ${entity_id}.coverage.tsv
        fi
    >>>

    output {
        File gatk_coverage_file = "${entity_id}.coverage.tsv"
    }

    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      zones: "us-central1-a us-central1-b"
      disks: "local-disk 200 SSD"
      memory: "6G"
    }
}

# Calculate coverage on Whole Genome Sequence using Spark.
# This task automatically creates a target output file.
task WholeGenomeCoverage {
    String entity_id
    File coverage_file
    File target_file
    File input_bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File gatk_jar
    Boolean isWGS
    Int wgsBinSize
    Int mem

    # If isWGS is set to true, the task produces WGS coverage and targets that are passed to downstream tasks
    # If not, coverage and target files (received from upstream) for WES are passed downstream
    # Note from Megan: I had to drop the -s from the ln to get this to work when isWGS is false 
    command {
        if [ ${isWGS} = true ]; \
          then java -Xmx${mem}g -jar ${gatk_jar} SparkGenomeReadCounts --outputFile ${entity_id}.coverage.tsv \
                --reference ${ref_fasta} --input ${input_bam} --sparkMaster local[1] --binsize ${wgsBinSize}; \
          else ln ${coverage_file} ${entity_id}.coverage.tsv; ln ${target_file} ${entity_id}.coverage.tsv.targets.tsv; \
        fi
    }

    output {
        File gatk_coverage_file = "${entity_id}.coverage.tsv"
        File gatk_target_file = "${entity_id}.coverage.tsv.targets.tsv"
    }

    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      zones: "us-central1-a us-central1-b"
      disks: "local-disk 200 SSD"
      memory: "6G"
    }
}

# Add new columns to an existing target table with various targets
# Note that this task is optional 
task AnnotateTargets {
    String entity_id
    File target_file
    File gatk_jar
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    Boolean enable_gc_correction
    Int mem

    # If GC correction is disabled, then an empty file gets passed downstream
    command {
        if [ ${enable_gc_correction} = true ]; \
          then java -Xmx${mem}g -jar ${gatk_jar} AnnotateTargets --targets ${target_file} --reference ${ref_fasta} --output ${entity_id}.annotated.tsv; \
          else touch ${entity_id}.annotated.tsv; \
        fi
    }

    output {
        File annotated_targets = "${entity_id}.annotated.tsv"
    }

    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      zones: "us-central1-a us-central1-b"
      disks: "local-disk 200 SSD"
      memory: "6G"
    }
}

# Correct coverage for sample-specific GC bias effects
# Note that this task is optional 
task CorrectGCBias {
    String entity_id
    File coverage_file
    File annotated_targets
    File gatk_jar
    Boolean enable_gc_correction
    Int mem

    # If GC correction is disabled, then the coverage file gets passed downstream unchanged
    command {
        if [ ${enable_gc_correction} = true ]; \
          then java -Xmx${mem}g -jar ${gatk_jar} CorrectGCBias --input ${coverage_file} \
           --output ${entity_id}.gc_corrected_coverage.tsv --targets ${annotated_targets}; \
          else ln -s ${coverage_file} ${entity_id}.gc_corrected_coverage.tsv; \
        fi
    }

    output {
        File gatk_cnv_coverage_file_gcbias = "${entity_id}.gc_corrected_coverage.tsv"
    }

    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      zones: "us-central1-a us-central1-b"
      disks: "local-disk 200 SSD"
      memory: "6G"
    }
}

# Perform tangent normalization (noise reduction) on the proportional coverage file.
#add tn file
task NormalizeSomaticReadCounts {
    String entity_id
    File coverage_file
    File padded_target_file
    File pon
    File gatk_jar
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} NormalizeSomaticReadCounts --input ${coverage_file} \
         --targets ${padded_target_file} --panelOfNormals ${pon} --factorNormalizedOutput ${entity_id}.fnt.tsv --tangentNormalized ${entity_id}.tn.tsv \
         --betaHatsOutput ${entity_id}.betaHats.tsv --preTangentNormalized ${entity_id}.preTN.tsv  --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File tn_file = "${entity_id}.tn.tsv"
        File pre_tn_file = "${entity_id}.preTN.tsv"
        File betahats_file = "${entity_id}.betaHats.tsv"
    }
    
    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      zones: "us-central1-a us-central1-b"
      disks: "local-disk 200 SSD"
      memory: "6G"
    }
}

# Segment the tangent normalized coverage profile.
task PerformSegmentation {
    String entity_id
    Float seg_param_alpha
    Int seg_param_nperm
    String seg_param_pmethod
    Int seg_param_minWidth
    Int seg_param_kmax
    Int seg_param_nmin
    Float seg_param_eta
    Float seg_param_trim
    String seg_param_undoSplits
    Float seg_param_undoPrune
    Int seg_param_undoSD
    File gatk_jar
    File tn_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} PerformSegmentation --tangentNormalized ${tn_file} \
         --output ${entity_id}.seg --log2Input true  --alpha ${seg_param_alpha} --nperm ${seg_param_nperm} \
         --pmethod ${seg_param_pmethod} --minWidth ${seg_param_minWidth} --kmax ${seg_param_kmax} \
         --nmin ${seg_param_nmin} --eta ${seg_param_eta} --trim ${seg_param_trim} --undoSplits ${seg_param_undoSplits} \
         --undoPrune ${seg_param_undoPrune} --undoSD ${seg_param_undoSD} --help false --version false \
         --verbosity INFO --QUIET false
    }

    output {
        File seg_file = "${entity_id}.seg"
    }

    runtime {
      docker: "mshand/genomesinthecloud:CNV-RPackages-1.2.4"
      zones: "us-central1-a us-central1-b"
      disks: "local-disk 200 SSD"
      memory: "6G"
    }
}

# Make calls (amp, neutral, or deleted) on each segment.
task Caller {
    String entity_id
    File gatk_jar
    File tn_file
    File seg_file
    Int mem

    command {
        java -Xmx${mem}g -jar ${gatk_jar} CallSegments --tangentNormalized ${tn_file} \
         --segments ${seg_file} --output ${entity_id}.called  --legacy false \
         --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File called_file="${entity_id}.called"
    }

    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
      zones: "us-central1-a us-central1-b"
      disks: "local-disk 200 SSD"
      memory: "6G"
    }
}

# Create plots of coverage data and copy-ratio estimates
task PlotSegmentedCopyRatio {
    String entity_id
    File gatk_jar
    File tn_file
    File pre_tn_file
    File called_file
    File ref_fasta_dict
    String output_dir
    Int mem

    command {
        mkdir -p ${output_dir} && \
        java -Xmx${mem}g -jar ${gatk_jar} PlotSegmentedCopyRatio --tangentNormalized ${tn_file} \
         --preTangentNormalized ${pre_tn_file} --segments ${called_file} \
         -SD ${ref_fasta_dict} \
         --output ${output_dir} --outputPrefix ${entity_id}
    }

    output {
        File segments_plot="${output_dir}/${entity_id}_FullGenome.png"
        File before_after_normalization_plot="${output_dir}/${entity_id}_Before_After.png"
        File before_after_cr_lim_4="${output_dir}/${entity_id}_Before_After_CR_Lim_4.png"
    }

    runtime {
      docker: "mshand/genomesinthecloud:CNV-RPackages-1.2.4"
      zones: "us-central1-a us-central1-b"
      disks: "local-disk 200 SSD"
      memory: "6G"
    }
}

task oncotate_m2_to_maf {
    File m2_vcf
    String entity_id
    File onco_ds_tar_gz

    command {
        # Downlaod and untar the db-dir
        set -e
        echo "Using given tar file: ${onco_ds_tar_gz}"
        tar zxvf ${onco_ds_tar_gz}
        ln -s oncotator_v1_ds_April052016_trimmed onco_dbdir

        /root/oncotator_venv/bin/oncotator --db-dir onco_dbdir/ -c $HOME/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt  -v \
        ${m2_vcf} ${entity_id}.oncotated.maf hg19 \
        -i VCF -o TCGAMAF --collapse-number-annotations --log_name oncotator.log \
            --skip-no-alt \
            -a Center:broadinstitute.org \
            -a tumor_barcode:${entity_id} \
            -a normal_barcode:NA \
            -a status:Somatic \
            -a NCBI_Build:37 \
            -a Strand:+ \
            -a source:WXS \
            -a phase:Phase_I \
            -a sequencer: \
            -a Tumor_Validation_Allele1: \
            -a Tumor_Validation_Allele2: \
            -a Match_Norm_Validation_Allele1: \
            -a Match_Norm_Validation_Allele2: \
            -a Verification_Status: \
            -a Validation_Status:Untested \
            -a Validation_Method:none \
            -a Score: \
            -a BAM_file: \
            -a Match_Norm_Seq_Allele1: \
            -a Match_Norm_Seq_Allele2:

            # $ cat /xchip/cga/reference/annotation/db/tcgaMAFManualOverrides2.4.config
            # [manual_annotations]
            # override:NCBI_Build=37,Strand=+,Center=broad.mit.edu,source=WXS,status=Somatic,phase=Phase_I,sequencer=Illumina GAIIx,Tumor_Validation_Allele1=,Tumor_Validation_Allele2=,Match_Norm_Validation_Allele1=,Match_Norm_Validation_Allele2=,Verification_Status=,Validation_Status=Untested,Validation_Method=none,Score=,BAM_file=,Match_Norm_Seq_Allele1=,Match_Norm_Seq_Allele2=
    }

    runtime {
        docker: "broadinstitute/oncotator:1.9.2.0"
        disks: "local-disk 400 SSD"
        bootDiskSizeGb: "20"
    }

    output {
        File oncotated_m2_maf="${entity_id}.oncotated.maf"
    }
}

task oncotate_cnv {
    String entity_id
    File target_seg_gt_file

    command {
        /root/oncotator_venv/bin/oncotator --db-dir /root/onco_dbdir/ -c /root/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt -u file:///root/onco_cache/ -r -v ${target_seg_gt_file} ${entity_id}.per_target.oncotated.txt hg19 -i SEG_FILE -o SIMPLE_TSV
    }

    output {
        File oncotated_target_seg_gt_file = "${entity_id}.per_target.oncotated.txt"
    }

    runtime {
        docker: "broadinstitute/oncotator:1.9.2.0-eval-gatk-protected"
        memory: "2GB"
        disks: "local-disk 400 SSD"
        bootDiskSizeGb: "20"
    }
}
