## Copyright Broad Institute, 2017
## 
## This WDL pipeline implements SV calling with Illumina's Manta software
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
task Manta {
  File input_bam
  File input_bam_index
  String output_vcf_basename
  File ref_fasta
  File ref_fasta_index

  command {
    /usr/bin/manta/bin/configManta.py \
    --bam ${input_bam} \
    --referenceFasta ${ref_fasta} \
    --runDir . &&

    ./runWorkflow.py \
    --mode local \
    --jobs 32 \
    --memGb 96 &&

    mv results/variants/diploidSV.vcf.gz ${output_vcf_basename}.diploidSV.vcf.gz &&
    mv results/variants/diploidSV.vcf.gz.tbi ${output_vcf_basename}.diploidSV.vcf.gz.tbi
  }
  runtime {
    docker: "eitanbanks/manta-1.0.3"
    memory: "100 GB"
    cpu: "32"
    disks: "local-disk 300 HDD"
    preemptible: 3
  }
  output {
    File output_vcf = "${output_vcf_basename}.diploidSV.vcf.gz"
    File output_vcf_index = "${output_vcf_basename}.diploidSV.vcf.gz.tbi"
  }
}

# PUBLIC #
# WORKFLOW DEFINITION
workflow MantaGermlineGenomeWorkflow {

  # PUBLIC #
  String sample_name
  String input_bam
  String input_bam_index
  File ref_fasta
  File ref_fasta_index

  call Manta as MantaSV {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
	ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
	output_vcf_basename = sample_name
  }

  output {
    MantaSV.*
  }
}
