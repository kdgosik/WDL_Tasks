task recapseg_coverage {
  File bed_file
  File bam_file
  String sample_id


  command {
    samtools index ${bam_file} && recapseg coverage ${bed_file} ${bam_file} ${sample_id}.coverage
  }
  runtime {
    docker: "adunford/recapseg"
    defaultDisks: "local-disk 100 SSD"

  }
  output {
    File recapseg_coverage_file = "${sample_id}.coverage"
  }
}

task recapseg_pc {
  File cov_file 
  File bed_file
  File bam_file
  String sample_id
  Float percent_zero


  command {
    python /src/recapseg/capseg/tools/generate_sample_read_group_file/generate_sample_read_group_file.py ${bam_file} ${sample_id}.sample_rg.tsv && recapseg pc ${cov_file}  ${bed_file} ${sample_id}.sample_rg.tsv ${sample_id} -pzf ${percent_zero} 
  }
  runtime {
    docker: "adunford/recapseg"
    defaultDisks: "local-disk 100 SSD"

  }
  output {
    File recapseg_pcov = "${sample_id}.pcov"
    File recapseg_cr_stat = "${sample_id}.pcov.cr.stat"
  }
}


task recapseg_tumor_pcov {
  File pon_file
  File recapseg_pcov
  File recapseg_cr_stat
  String sample_id
  File bed_file
  String is_plotting
  String no_pon_reduction

  command {
    recapseg tumor_pcov ${pon_file} ${recapseg_pcov} ${recapseg_cr_stat} ${sample_id} ${bed_file} recapseg ${is_plotting} ${no_pon_reduction}
  }
  runtime {
    docker: "adunford/recapseg"
    memory: "12 GB"
    defaultDisks: "local-disk 100 SSD"

  }

  output {
    File recapseg_probe_file = "recapseg.tn.${sample_id}.tsv"
    File recapseg_seg_file = "recapseg.${sample_id}.seg"
  }
}



workflow recapseg_workflow {
  File BAM
  File ReCapSeg_Bed
  File Coverage_Bed
  String id
  
  call recapseg_coverage {input: sample_id = id, bam_file=BAM, bed_file=Coverage_Bed}
  call recapseg_pc {input: cov_file = recapseg_coverage.recapseg_coverage_file, bed_file=ReCapSeg_Bed, bam_file=BAM, sample_id=id }
  call recapseg_tumor_pcov {input: bed_file=ReCapSeg_Bed, sample_id=id, recapseg_pcov = recapseg_pc.recapseg_pcov, recapseg_cr_stat = recapseg_pc.recapseg_cr_stat}
}