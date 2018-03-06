task plink {
  File input_vcf
  String basename

  command <<<
    plink --vcf ${input_vcf} --out ${basename}
  >>>

  runtime {
    disks: "local-disk 50 HDD"
    memory: "3500 MB"
    docker: "jrose77/plinkdocker"
  }

  output {
    File bed_file = "${basename}.bed"
    File bim_file = "${basename}.bim"
    File fam_file = "${basename}.fam"
  }
}

workflow plink_workflow {
  File input_vcf
  String basename

  call plink {
    input:
      input_vcf=input_vcf,
      basename=basename
  }

  output {
    plink.*
  }
}