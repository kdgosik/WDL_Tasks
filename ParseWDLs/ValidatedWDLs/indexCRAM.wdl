## Copyright Broad Institute, 2017
## 
## Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
## 
## WDL file to index a single CRAM file
##
## Requirements/expectations:
## - Whole-genome sequencing data in mapped CRAM format
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

#Task to index CRAM file
task indexCRAM {
  File input_cram
  String cram_index_prefix

  command {
    samtools index ${input_cram} ${cram_index_prefix}.cram.crai
  }
  runtime {
    docker: "biocontainers/samtools"
    memory: "4 GB"
    cpu: "1"
    disks: "local-disk 300 SSD"
  }
  output {
    File cram_index = "${cram_index_prefix}.cram.crai"
  }
}

#Workflow to index CRAM file
workflow indexCRAMWorkflow {
  File input_cram
  String cram_index_prefix

  call indexCRAM as indexCRAM {
    input:
      input_cram = input_cram,
      cram_index_prefix = cram_index_prefix
  }
  output {
    indexCRAM.*
  } 
}