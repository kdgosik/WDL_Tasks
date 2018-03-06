## Copyright Broad Institute, 2017
## 
## This WDL pipeline implements SV calling with DELLY2 by Tobias Rausch (https://github.com/dellytools/delly)
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in mapped BAM format
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# PUBLIC #

task GatherBCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_name
  
  command <<<
    bcftools concat ${sep=' ' input_vcfs} -n -o ${output_vcf_name} 
    bcftools index ${output_vcf_name}
  >>>
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.csi"
  }
  runtime {
    docker: "biocontainers/bcftools"
    memory: "2 GB"
    cpu: 1
    disks: "local-disk 50 HDD"
    preemptible: 3
  }
}

task DellyDiscoverEventtype {
  File input_bam
  File input_bam_index
  String output_vcf_basename
  File ref_fasta
  File ref_fasta_index
  File intervals_tsv
  String eventtype

  command {
    /usr/delly_v0.7.7/delly_v0.7.7_CentOS5.4_x86_64bit \
	call -t ${eventtype} \
	-g ${ref_fasta} \
	-x ${intervals_tsv} \
	-o ${output_vcf_basename}.${eventtype}.bcf \
	-n \
	${input_bam}
  }
  runtime {
    docker: "ldgauthier/delly_v0.7.7"
    memory: "2 GB"
    cpu: "1"
    disks: "local-disk 300 HDD"
    preemptible: 3
  }
  output {
    File output_vcf = "${output_vcf_basename}.${eventtype}.bcf"
    File output_vcf_index = "${output_vcf_basename}.${eventtype}.bcf.csi"
  }
}

# PUBLIC #
# WORKFLOW DEFINITION
workflow DellyGermlineGenomeWorkflow {

  # PUBLIC #
  String sample_name
  String input_bam
  String input_bam_index
  File ref_fasta
  File ref_fasta_index
  Array[File] scatterIntervals

  scatter (shardIntervals in scatterIntervals) {

  call DellyDiscoverEventtype as DellyDEL {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
	 ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
	 intervals_tsv = shardIntervals,
	 output_vcf_basename = sample_name,
	 eventtype = "DEL"
  }

  call DellyDiscoverEventtype as DellyDUP {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
	ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
	intervals_tsv = shardIntervals,
	output_vcf_basename = sample_name,
	eventtype = "DUP"
  }

  call DellyDiscoverEventtype as DellyINV {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
	ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
	intervals_tsv = shardIntervals,
	output_vcf_basename = sample_name,
	eventtype = "INV"
  }
}

  call GatherBCFs as gatherDEL {
    input:    
    input_vcfs = DellyDEL.output_vcf,
    input_vcfs_indexes = DellyDEL.output_vcf_index,
    output_vcf_name = "${sample_name}.DEL.bcf"
}

  call GatherBCFs as gatherDUP {
    input:
    input_vcfs = DellyDUP.output_vcf,
    input_vcfs_indexes = DellyDUP.output_vcf_index,
    output_vcf_name = "${sample_name}.DUP.bcf"
}

  call GatherBCFs as gatherINV {
    input:
    input_vcfs = DellyINV.output_vcf,
    input_vcfs_indexes = DellyINV.output_vcf_index,
    output_vcf_name = "${sample_name}.INV.bcf"
}


  output {
    gatherDEL.*
    gatherDUP.*
    gatherINV.*
  }
}
