import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:01_pesr_clustering_single_algorithm/versions/2/plain-WDL/descriptor" as single

workflow cluster_pesr {
  Array[File] manta_vcfs
  Array[File] delly_vcfs
  Array[File] melt_vcfs
  File contigs
  String batch
  
  Int dist
  Float frac
  File blacklist
  Int svsize
  String flags

  call single.cluster_pesr_algorithm as cluster_manta {
    input:
      vcfs=manta_vcfs,
      batch=batch,
      algorithm="manta",
      contigs=contigs,
      dist=dist,
      frac=frac,
      blacklist=blacklist,
      svsize=svsize,
      flags=flags
  }

  call single.cluster_pesr_algorithm as cluster_delly {
    input:
      vcfs=delly_vcfs,
      batch=batch,
      algorithm="delly",
      contigs=contigs,
      dist=dist,
      frac=frac,
      blacklist=blacklist,
      svsize=svsize,
      flags=flags
  }
  
  call single.cluster_pesr_algorithm as cluster_melt {
    input:
      vcfs=melt_vcfs,
      batch=batch,
      algorithm="melt",
      contigs=contigs,
      dist=dist,
      frac=frac,
      blacklist=blacklist,
      svsize=svsize,
      flags=flags
  }

  output {
    File manta_vcf = cluster_manta.clustered_vcf
    File delly_vcf = cluster_delly.clustered_vcf
    File melt_vcf = cluster_melt.clustered_vcf
  }
}