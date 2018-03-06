workflow cluster_pesr_algorithm {
  Array[File] vcfs
  File contigs
  String batch
  String algorithm

  # VCFcluster parameters
  Int dist
  Float frac
  File blacklist
  Int svsize
  String flags

  Array[Array[String]] contiglist = read_tsv(contigs)

  scatter (contig in contiglist) {
    call vcfcluster {
      input:
        vcfs=vcfs,
        batch=batch,
        algorithm=algorithm,
        chrom=contig[0],
        dist=dist,
        frac=frac,
        blacklist=blacklist,
        svsize=svsize,
        flags=flags
    }
  }

  call concat_vcfs {
    input:
      vcfs=vcfcluster.clustered_vcf,
      batch=batch,
      algorithm=algorithm
  }
  
  output {
    File clustered_vcf = concat_vcfs.vcf
  }
}

task vcfcluster {
  Array[File] vcfs
  String batch
  String algorithm
  String chrom

  # VCFcluster parameters
  Int dist
  Float frac
  File blacklist
  Int svsize
  String flags

  command <<<
    for f in ${sep=' ' vcfs}; do tabix -p vcf -f $f; done;
    tabix -p bed ${blacklist};

    svtk vcfcluster ${write_tsv(vcfs)} stdout \
        -r ${chrom} \
        -p ${batch}_${algorithm}_${chrom} \
        -d ${dist} \
        -f ${frac} \
        -x ${blacklist} \
        -z ${svsize} \
        ${flags} \
      | vcf-sort -c \
      | bgzip -c > ${batch}.${algorithm}.${chrom}.vcf.gz
  >>>

  output {
    File clustered_vcf="${batch}.${algorithm}.${chrom}.vcf.gz"
  }
  
  runtime {
    docker: "msto/sv-pipeline"
  }
}

task concat_vcfs {
  Array[File] vcfs
  String batch
  String algorithm

  command {
    vcf-concat ${sep=' ' vcfs} | vcf-sort -c | bgzip -c > ${batch}.${algorithm}.vcf.gz;
    tabix -p vcf ${batch}.${algorithm}.vcf.gz;
  }

  output {
    File vcf="${batch}.${algorithm}.vcf.gz"
    File idx="${batch}.${algorithm}.vcf.gz.tbi"
  }
  
  runtime {
    docker: "msto/sv-pipeline"
  }
}