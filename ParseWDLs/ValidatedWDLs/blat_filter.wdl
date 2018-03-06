task realign {
  File bam
  File bai
  File maf
  String id

  command {
    python /opt/realign.py ${bam} ${maf} ${id}
}

  runtime {
    memory: "12 GB"
    disks: "local-disk 100 HDD" 
    bootDiskSizeGb: 30
    #docker: "dlivitzbroad/blat_filter_v2"
    docker: "gcr.io/broad-firecloud-itools/blat_filter_v2:mate_fix_2"
    preemptible: 1
  }

  output {
    File blat_results = "${id}.blat.maf"
    File debug_results = "${id}.blat.rejected.maf"
    File blat_all_maf = "${id}.blat.all.maf"
  }
}

workflow realign_wkfl {
  call realign
}