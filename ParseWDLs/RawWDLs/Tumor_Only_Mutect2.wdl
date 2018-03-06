#Must be hardcoded outputs due to difficulties with glob.
task SplitIntervals {
  Int disk_size
  String picard_jar
  File intervals

  command {
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
  File pon
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
      --normal_panel ${pon} \
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

workflow TumorOnly  {

  File intervals 
  File ref_fasta
  File ref_fasta_idx
  File ref_dict
  File tumor_bam
  File tumor_bam_idx
  File unmatched_normal_bam
  File unmatched_normal_bam_idx
  String tumor_sample_name
  File cosmic
  File cosmicIndex
  File dbsnp
  File dbsnpIndex
  File pon
  File exac
  File exac_idx
  File onekg
  Int medium_disk
  Int large_disk
  String GATK_jar
  String picard_jar

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
        pon = pon,
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

  call CollectMetrics {
    input:
      input_vcf = GatherVCFs.output_vcf,
      input_vcf_index = GatherVCFs.output_vcf_index,
      cosmic = cosmic,
      cosmicIndex = cosmicIndex,
      intervals = intervals,
      ref_dict = ref_dict,
      output_basename = tumor_sample_name,
      disk_size = medium_disk,
      picard_jar = picard_jar
  }
}
