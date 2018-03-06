task custom {
  File script
  String input1
  File? input2
  File? input3
  Array[File]? input4
  String docker
  String memory
  String disk
  String preempt

  command {
     chmod +x "${script}"
     ${script} "${input1}" ${input2} ${input3} ${sep=" " input4}
  }

  runtime {
    memory: "${memory} GB"
    docker: "${docker}"
    disks: "local-disk  ${disk}  HDD"
    preemptible: "${preempt}"
  }

  output {
    Array[File] all_outs=glob("*")
  }
}

workflow custm_wkfl {
  call custom
}