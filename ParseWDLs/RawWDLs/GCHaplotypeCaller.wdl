workflow HaplotypeCallerWf {

	File wgs_calling_interval_list = "gs://broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list"
	Int haplotype_scatter_count = 50
	Int break_bands_at_multiples_of = 1000000

	File input_bam
	File input_bam_index
	String sample_name
	File ref_dict
	File ref_fasta
	File ref_fasta_index
	Int disk_size


	# Break the calling interval_list into sub-intervals
  	# Perform variant calling on the sub-intervals, and then gather the results
	call ScatterIntervalList {
	 	input:
	 		interval_list = wgs_calling_interval_list,
	    	scatter_count = haplotype_scatter_count,
	      	break_bands_at_multiples_of = break_bands_at_multiples_of,
	      	disk_size = disk_size
	}

  	# Call variants in parallel over WGS calling intervals
	scatter (index in range(ScatterIntervalList.interval_count)) {

		call HaplotypeCaller {
			input:
				input_bam = input_bam,
				input_bam_index = input_bam_index,
				interval_list = ScatterIntervalList.out[index],
				gvcf_basename = sample_name,
				ref_dict = ref_dict,
				ref_fasta = ref_fasta,
				ref_fasta_index = ref_fasta_index,
				disk_size = disk_size
		}
	}

	call MergeVCFs {
        input:
          input_vcfs = HaplotypeCaller.output_gvcf,
          input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
          output_vcf_name = sample_name + ".g.vcf.gz",
          disk_size = disk_size
    }


    output {
        File gvcf = MergeVCFs.output_vcf
        File gvcf_index = MergeVCFs.output_vcf_index
    }

	
}


# should really do contamination but

## HaplotypeCaller from SingleSample Pipeline
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size

  # We use interval_padding 500 below to make sure that the HaplotypeCaller has context on both sides around
  # the interval because the assembly uses them.
  #
  # Using PrintReads is a temporary solution until we update HaploypeCaller to use GATK4. Once that is done,
  # HaplotypeCaller can stream the required intervals directly from the cloud.
  command {

    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m \
      -jar /usr/gitc/GATK35.jar \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${gvcf_basename}.vcf.gz \
      -I ${input_bam} \
      -L ${interval_list} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -variant_index_parameter 128000 \
      -variant_index_type LINEAR \
      -contamination 0 \
      --read_filter OverclippedRead
  }
  runtime {
    preemptible: 0
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:7bc64948a0a9f50ea55edb8b30c710943e44bd861c46a229feaf121d345e68ed"

  }

  output {
    File output_gvcf = "${gvcf_basename}.vcf.gz"
    File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
  }
}

# This task calls picard's IntervalListTools to scatter the input interval list into scatter_count sub interval lists
# Note that the number of sub interval lists may not be exactly equal to scatter_count.  There may be slightly more or less.
# Thus we have the block of python to count the number of generated sub interval lists.
task ScatterIntervalList {
  File interval_list
  Int scatter_count
  Int break_bands_at_multiples_of
  Int disk_size

  command <<<
    set -e
    mkdir out
    java -Xms1g -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT=${scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=${break_bands_at_multiples_of} \
      INPUT=${interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
    Int interval_count = read_int(stdout())
  }
  runtime {
    memory: "2 GB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:7bc64948a0a9f50ea55edb8b30c710943e44bd861c46a229feaf121d345e68ed"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_name
  Int disk_size

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
      MergeVcfs \
      INPUT=${sep=' INPUT=' input_vcfs} \
      OUTPUT=${output_vcf_name}
  }

  runtime {
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:7bc64948a0a9f50ea55edb8b30c710943e44bd861c46a229feaf121d345e68ed"
    cpu: "1"
  }

  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}