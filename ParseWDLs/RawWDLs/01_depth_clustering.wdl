import "https://api.firecloud.org/ga4gh/v1/tools/Talkowski-SV:01_depth_clustering_by_chrom/versions/1/plain-WDL/descriptor" as dibc

workflow cluster_depth {
  File del_bed
  File dup_bed
  File contigs
  Float frac
  String flags
  String batch

  call dibc.bedcluster_by_chrom as cluster_DELs {
    input:
      batch=batch,
      svtype="DEL",
      bed=del_bed,
      contigs=contigs,
      frac=frac,
      flags=flags
  }

  call dibc.bedcluster_by_chrom as cluster_DUPs {
    input:
      batch=batch,
      svtype="DUP",
      bed=dup_bed,
      contigs=contigs,
      frac=frac,
      flags=flags
  }

  call make_rdtest_bed {
    input:
      dels=cluster_DELs.clustered_bed,
      dups=cluster_DUPs.clustered_bed,
      batch=batch,
  }

  call make_depth_vcf {
    input:
      bed=make_rdtest_bed.bed,
      batch=batch,
      contigs=contigs
  }

  output {
    File clustered_vcf = make_depth_vcf.vcf
  }
}

task make_rdtest_bed { 
  File dels
  File dups
  File script
  String batch

  command <<<
    cat \
        <(python3 ${script} ${dels} | sed '1d') \
        <(python3 ${script} ${dups} | sed '1d') \
      | sort -k1,1V -k2,2n \
      | cat <(echo -e "#chrom start end name samples svtype" | sed -e 's/ /\t/g') - \
      > ${batch}.depth.bed;
  >>>
  
  output {
    File bed = "${batch}.depth.bed"
  }
  
  runtime {
    docker: "msto/sv-pipeline"
  }
}

task make_depth_vcf {
  File bed
  File contigs
  String batch
  
  command <<<
    cut -f5 ${bed} | sed -e '1d' -e 's/,/\n/g' | sort -u > samples.list;
    svtk rdtest2vcf --contigs ${contigs} ${bed} samples.list ${batch}.depth.vcf.gz;
  >>>

  output {
    File vcf = "${batch}.depth.vcf.gz"
  }
  
  runtime {
    docker: "msto/sv-pipeline"
  }
}