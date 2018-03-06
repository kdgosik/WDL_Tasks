task IntervalListToArrayInterval {
    File interval_list

    command <<<
        grep -v @ ${interval_list} | awk -F'\t' '{print $1":"$2"-"$3}' > interval_per_line
      >>>

    runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }

    output {
        Array[String] array = read_lines("interval_per_line")
    }
}

task IntToIntArray {
    Int number

    command <<<
        python <<CODE
        for i in range(${number}):
          print(i)
        CODE
      >>>

    runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }

    output {
        Array[Int] array = read_lines(stdout())
    }
}

task FileToStringArray {
    File input_file

    command {

    }

    runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }

    output {
        Array[String] array = read_lines(input_file)
    }
}

task SplitByLine {
  File input_file

  command {
    split -l1 -a10 ${input_file} 'split-'
  }
  output {
    Array[File] files = glob("split-*")
  }
    runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }
}

## TODO - add useInbreedingCoefficientFilter as workflow input and use true and false wdl function to input correct command
# filterExpression and filtername should only be added to the commandline if useInbreedingCoefficientFilter is true
task GenotypeFilterAndMakeSitesOnlyVcf {
  File gvcf

  File interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File dbsnp_vcf
  File dbsnp_vcf_index

  String variant_filtered_vcf_filename
  String variant_filtered_vcf_index_filename

  String split_interval

  String sites_only_vcf_filename
  String sites_only_vcf_index_filename

  Int disk_size
  
  command <<<
    /usr/gitc/tabix ${gvcf}

    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx11500m \
      -jar /usr/gitc/GATK36.jar \
      -T GenotypeGVCFs \
      --disable_auto_index_creation_and_locking_when_reading_rods \
      -R ${ref_fasta} \
      -o tmp.unfiltered.vcf.gz \
      -V ${gvcf} \
      -L ${split_interval} \
      -D ${dbsnp_vcf}

    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx11500m \
      -jar /usr/gitc/GATK35.jar \
      -T VariantFiltration \
      --disable_auto_index_creation_and_locking_when_reading_rods \
      --filterExpression "InbreedingCoeff < -0.3" \
      --filterName InbreedingCoeff \
      -R ${ref_fasta} \
      -o ${variant_filtered_vcf_filename} \
      -V tmp.unfiltered.vcf.gz \
      -L ${interval_list}

    java -Xmx11500m -jar /usr/gitc/picard.jar \
      MakeSitesOnlyVcf \
      INPUT=${variant_filtered_vcf_filename} \
      OUTPUT=${sites_only_vcf_filename}

    >>>
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.2-1466113830"
    memory: "13 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File variant_filtered_vcf = "${variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "${variant_filtered_vcf_index_filename}"
    File sites_only_vcf = "${sites_only_vcf_filename}"
    File sites_only_vcf_index = "${sites_only_vcf_index_filename}"
  }
}

## TODO - Combine following two Variant Recalibrators into one once map iteration similar to array 'sep' or looping is implemented
# the -resource file inputs need to be localised and pulled into the workflow and currently this is the only viable way of doing it
# until one of the above mentioned syntax is implemented

task SNPSVariantRecalibrator {
  String recalibration_filename
  String recalibration_index_filename
  String tranches_filename
  String recalibration_plots_rscript_filename
  String mode

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  Array[File] sites_only_variant_filtered_vcf_list
  Array[File] sites_only_variant_filtered_vcf_index_list

  File interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File hapmap_resource_vcf
  File omni_resource_vcf
  File one_thousand_genomes_resource_vcf
  File dbsnp_resource_vcf
 
  File hapmap_resource_vcf_index
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf_index
  File dbsnp_resource_vcf_index
 

  Int disk_size

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx32000m \
      -jar /usr/gitc/GATK35.jar \
      -T VariantRecalibrator \
      -R ${ref_fasta} \
      -input ${sep=' -input 'sites_only_variant_filtered_vcf_list} \
      -L ${interval_list} \
      --disable_auto_index_creation_and_locking_when_reading_rods \
      -recalFile ${recalibration_filename} \
      -tranchesFile ${tranches_filename} \
      -allPoly \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode ${mode} \
      --maxGaussians 6 \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_resource_vcf} \
      -resource:omni,known=false,training=true,truth=true,prior=12 ${omni_resource_vcf} \
      -resource:1000G,known=false,training=true,truth=false,prior=10 ${one_thousand_genomes_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbsnp_resource_vcf} \
      -rscriptFile ${recalibration_plots_rscript_filename}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.2-1466113830"
    memory: "208 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_index_filename}"
    File tranches = "${tranches_filename}"
    File recalibration_plots_rscript = "${recalibration_plots_rscript_filename}"
    File recalibration_plots_pdf = "${recalibration_plots_rscript_filename}.pdf"
  }
}

task IndelsVariantRecalibrator {
  String recalibration_filename
  String recalibration_index_filename
  String tranches_filename
  String recalibration_plots_rscript_filename
  String mode

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  Array[File] sites_only_variant_filtered_vcf_list
  Array[File] sites_only_variant_filtered_vcf_index_list

  File interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File mills_resource_vcf
  File axiomPoly_resource_vcf
  File dbsnp_resource_vcf
  File mills_resource_vcf_index
  File axiomPoly_resource_vcf_index
  File dbsnp_resource_vcf_index

  Int disk_size

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx32000m \
      -jar /usr/gitc/GATK35.jar \
      -T VariantRecalibrator \
      -R ${ref_fasta} \
      -input ${sep=' -input 'sites_only_variant_filtered_vcf_list} \
      -L ${interval_list} \
      --disable_auto_index_creation_and_locking_when_reading_rods \
      -recalFile ${recalibration_filename} \
      -tranchesFile ${tranches_filename} \
      -allPoly \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode ${mode} \
      --maxGaussians 4 \
      -resource:mills,known=false,training=true,truth=true,prior=12 ${mills_resource_vcf} \
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomPoly_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp_resource_vcf} \
      -rscriptFile ${recalibration_plots_rscript_filename}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.2-1466113830"
    memory: "208 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_index_filename}"
    File tranches = "${tranches_filename}"
    File recalibration_plots_rscript = "${recalibration_plots_rscript_filename}"
  }
}

task ApplyRecalibration {
  String recalibrated_vcf_filename
  File input_vcf
  File input_vcf_index
  File interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File indels_recalibration
  File indels_recalibration_index
  File indels_tranches
  File snps_recalibration
  File snps_recalibration_index
  File snps_tranches

  Float indel_filter_level
  Float snp_filter_level

  Int disk_size

  command <<<
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m \
      -jar /usr/gitc/GATK35.jar \
      -T ApplyRecalibration \
      --disable_auto_index_creation_and_locking_when_reading_rods \
      -R ${ref_fasta} \
      -o tmp.indel.recalibrated.vcf.gz \
      -input ${input_vcf} \
      -L ${interval_list} \
      -recalFile ${indels_recalibration} \
      -tranchesFile ${indels_tranches} \
      -ts_filter_level ${indel_filter_level} \
      -mode INDEL &&

    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m \
      -jar /usr/gitc/GATK35.jar \
      -T ApplyRecalibration \
      --disable_auto_index_creation_and_locking_when_reading_rods \
      -R ${ref_fasta} \
      -o ${recalibrated_vcf_filename} \
      -input tmp.indel.recalibrated.vcf.gz \
      -L ${interval_list} \
      -recalFile ${snps_recalibration} \
      -tranchesFile ${snps_tranches} \
      -ts_filter_level ${snp_filter_level} \
      -mode SNP
  >>>
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.2-1466113830"
    memory: "7 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibrated_vcf = "${recalibrated_vcf_filename}"
    File recalibrated_vcf_index = "${recalibrated_vcf_filename}.tbi"
  }
}

task CalculateGenotypePosteriors {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File input_vcf
  File input_vcf_index
  File variant_priors_vcf
  File variant_priors_vcf_index
  String output_vcf_filename
  Int disk_size 

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m \
    -jar /usr/gitc/GATK35.jar \
    -R ${ref_fasta} \
    -T CalculateGenotypePosteriors \
    -V ${input_vcf} \
    --supporting ${variant_priors_vcf} \
    --numRefSamplesIfNoCall 2500 \
    -o ${output_vcf_filename}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.2-1466113830"
    memory: "7 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}


# we would need to be able to variably set the runtime cpu value to mimic this.
# GatherVcfs CREATE_INDEX must be false - index creation not supported for gzipped files
task GatherVcfs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_name
  
  Int disk_size

  command <<<
    java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx1000m \
      -jar /usr/gitc/picard.jar \
      GatherVcfs \
      INPUT=${sep=' INPUT=' input_vcfs} \
      OUTPUT=${output_vcf_name} \
      CREATE_INDEX=false

    /usr/gitc/tabix ${output_vcf_name}
  >>>
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.2-1466113830"
    memory: "3 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " HDD"
  }
}

# on-prem the CollectVariantCallingMetrics THREAD_COUNT is a range of 8-16 threads depending on the number of slots SGE gives the job
# we would need to be able to variably set the runtime cpu value to mimic this.
task CollectVariantCallingMetrics {
  File input_vcf
  File input_vcf_index
  
  String metrics_filename_prefix
  File dbsnp_vcf
  File dbsnp_vcf_index
  File interval_list
  File ref_dict
  Int disk_size

  command <<<
    java -Xmx6000m -jar /usr/gitc/picard.jar \
      CollectVariantCallingMetrics \
      INPUT=${input_vcf} \
      DBSNP=${dbsnp_vcf} \
      SEQUENCE_DICTIONARY=${ref_dict} \
      OUTPUT=${metrics_filename_prefix} \
      THREAD_COUNT=32 \
      TARGET_INTERVALS=${interval_list}
  >>>
  output {
    File detail_metrics_file = "${metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${metrics_filename_prefix}.variant_calling_summary_metrics"
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.2-1466113830"
    memory: "7 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " HDD"
  }
}

workflow JointGenotyping {
  String callset_name

  Array[String] snp_recalibration_tranche_values
  Array[String] snp_recalibration_annotation_values
  Array[String] indel_recalibration_tranche_values
  Array[String] indel_recalibration_annotation_values

  File call_interval_list
  File eval_interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File dbsnp_vcf
  File dbsnp_vcf_index
  File hapmap_resource_vcf
  File omni_resource_vcf
  File one_thousand_genomes_resource_vcf
  File dbsnp_resource_vcf = dbsnp_vcf
  File mills_resource_vcf
  File axiomPoly_resource_vcf
  File panel_for_posteriors_vcf

  File hapmap_resource_vcf_index
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf_index
  File dbsnp_resource_vcf_index = dbsnp_vcf_index
  File mills_resource_vcf_index
  File axiomPoly_resource_vcf_index
  File panel_for_posteriors_vcf_index


  Float snp_filter_level
  Float indel_filter_level

  Int small_disk
  Int medium_disk
  Int large_disk

  File combine_gvcf_output_file
  File split_interval_list
  Int num_partitions


 call SplitByLine as gvcf_array_split { input: input_file=combine_gvcf_output_file } 


  call IntToIntArray {
    input:
      number = num_partitions
  }

  call IntervalListToArrayInterval {
    input:
      interval_list = split_interval_list
  }

  scatter (idx in IntToIntArray.array) {

    call GenotypeFilterAndMakeSitesOnlyVcf {
      input:
        gvcf = read_lines(gvcf_array_split.files[idx])[0],
        split_interval = IntervalListToArrayInterval.array[idx],
        interval_list = call_interval_list,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        variant_filtered_vcf_filename = callset_name + ".variant_filtered.vcf.gz",
        variant_filtered_vcf_index_filename = callset_name + ".variant_filtered.vcf.gz.tbi",
        sites_only_vcf_filename = callset_name + ".sites_only.variant_filtered.vcf.gz",
        sites_only_vcf_index_filename = callset_name + ".sites_only.variant_filtered.vcf.gz.tbi",
        disk_size = medium_disk
    }
  }

  call SNPSVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf_list = GenotypeFilterAndMakeSitesOnlyVcf.sites_only_vcf,
      sites_only_variant_filtered_vcf_index_list = GenotypeFilterAndMakeSitesOnlyVcf.sites_only_vcf_index,
      interval_list = call_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      recalibration_filename = callset_name + ".snps.recal",
      recalibration_index_filename = callset_name + ".snps.recal.idx",
      tranches_filename = callset_name + ".snps.tranches",
      recalibration_plots_rscript_filename = callset_name + ".snps.recalibration_plots.rscript",
      recalibration_tranche_values = snp_recalibration_tranche_values,
      recalibration_annotation_values = snp_recalibration_annotation_values,
      mode = "SNP",
      hapmap_resource_vcf = hapmap_resource_vcf,
      hapmap_resource_vcf_index = hapmap_resource_vcf_index,
      omni_resource_vcf = omni_resource_vcf,
      omni_resource_vcf_index = omni_resource_vcf_index,
      one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
      one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      disk_size = small_disk
  }

  call IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf_list = GenotypeFilterAndMakeSitesOnlyVcf.sites_only_vcf,
      sites_only_variant_filtered_vcf_index_list = GenotypeFilterAndMakeSitesOnlyVcf.sites_only_vcf_index,
      interval_list = call_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      recalibration_filename = callset_name + ".indels.recal",
      recalibration_index_filename = callset_name + ".indels.recal.idx",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_plots_rscript_filename = callset_name + ".indels.recalibration_plots.rscript",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mode = "INDEL",
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      disk_size = small_disk
  }

  scatter (idx in IntToIntArray.array) {

    call ApplyRecalibration {
      input:
        recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
        input_vcf = GenotypeFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
        input_vcf_index = GenotypeFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
        interval_list = call_interval_list,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = SNPSVariantRecalibrator.recalibration,
        snps_recalibration_index = SNPSVariantRecalibrator.recalibration_index,
        snps_tranches = SNPSVariantRecalibrator.tranches,
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level,
        disk_size = small_disk
    }

    call CalculateGenotypePosteriors {
      input:
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        input_vcf = ApplyRecalibration.recalibrated_vcf,
        input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
        variant_priors_vcf = panel_for_posteriors_vcf,
        variant_priors_vcf_index = panel_for_posteriors_vcf_index,
        output_vcf_filename = callset_name + ".filtered.posteriors." + idx + ".vcf.gz",
        disk_size = small_disk
   }  

  }
 
  call GatherVcfs {
    input:
      input_vcfs = ApplyRecalibration.recalibrated_vcf,
      input_vcfs_indexes = ApplyRecalibration.recalibrated_vcf_index,
      output_vcf_name = callset_name + ".vcf.gz",
      disk_size = large_disk
  }

  call CollectVariantCallingMetrics {
      input:
        input_vcf = GatherVcfs.output_vcf,
        input_vcf_index = GatherVcfs.output_vcf_index,
        metrics_filename_prefix = callset_name,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict,
        disk_size = small_disk
  }

  output {
    GatherVcfs.*
    CollectVariantCallingMetrics.*
  }
}
