task addIndexPath {
  String bamPath

  command {
    echo ${bamPath} | sed 's/bam$/bai/'
  }
  runtime {
    docker: "ldgauthier/delly_v0.7.7"
    memory: "2 GB"
    cpu: "1"
    disks: "local-disk 10 SSD"
  }
  output {
    String crai_or_bai_path = read_string(stdout())
  }
}

workflow VIRGOaddBamIndexes {
  String input_bam

  call addIndexPath as add {
    input:
      bamPath = input_bam
  }

  output {
    add.*
  }
    
}