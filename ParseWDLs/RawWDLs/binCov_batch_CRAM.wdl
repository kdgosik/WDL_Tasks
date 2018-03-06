## Copyright Broad Institute, 2017
## 
## Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
## 
## This WDL pipeline implements batched sequencing coverage metadata collection with binCov_batch.sh
##
## Requirements/expectations :
## - Whole-genome sequencing data in mapped CRAM format
## - CRAM index (.crai)
## - Exclusion list (BED-formatted) for regions to ignore
## - List of contigs to be evaluated
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# Task to run binCov_batch on a single sample on a list of chromosomes
task binCov {
  # binCov_batch arguments
  File input_cram
  File input_cram_index
  String sampleID
  File contigList
  File exclusion_bed

  command {
    binCov_batch.sh -z -C \
    -I ${input_cram_index} \
    -b 100 \
    -m nucleotide \
    -L ${contigList} \
    -x ${exclusion_bed} \
    -v 0.05 \
    ${input_cram} \
    ${sampleID} \
    ${sampleID}_binCov_batch_out
  }
  runtime {
    docker: "rlcollins/bincov"
    memory: "3 GB"
    cpu: 1
    disks: "local-disk 300 HDD"
    preemptible: 3
  }
  output {
    File bincov_out = "${sampleID}_binCov_batch_out.tar.gz"
  }
}

# Workflow to run binCov on a single sample
workflow binCovBatchWorkflow {

  # Workflow arguments
  File input_cram
  File input_cram_index
  String sampleID
  File exclusion_bed
  File contigList

  call binCov as binCov {
      input:
        input_cram = input_cram,
        input_cram_index = input_cram_index,
      	sampleID = sampleID,
        contigList = contigList,
        exclusion_bed = exclusion_bed
  }
  output {
    binCov.*
  }
}
