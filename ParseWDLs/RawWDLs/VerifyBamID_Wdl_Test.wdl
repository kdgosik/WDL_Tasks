# PRIVATE #
task CheckContamination {
  File input_bam
  File input_bam_index
  File contamination_sites_ud
  File contamination_sites_bed
  File contamination_sites_mu
  File ref_fasta
  File ref_fasta_index
  String output_prefix
  Int disk_size
  Int preemptible_tries

  # Having to do this as a 2-step command in heredoc syntax, adding a python step to read the metrics
  # This is a hack until read_object() is supported by Cromwell.
  # It relies on knowing that there is only one data row in the 2-row selfSM TSV file
  # Piping output of verifyBamId to /dev/null so only stdout is from the python command
  command <<<
    set -e

    /usr/gitc/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output ${output_prefix} \
    --BamFile ${input_bam} \
    --Reference ${ref_fasta} \
    --UDPath ${contamination_sites_ud} \
    --MeanPath ${contamination_sites_mu} \
    --BedPath ${contamination_sites_bed} \
    1>/dev/null

    python3 <<CODE
    import csv
    import sys
    with open('${output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
    # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event. 
    # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between 
    # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"])/0.75)
        i = i + 1
    # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
    # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "2 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "broadinstitute/verify-bam-id:c8a66425c312e5f8be46ab0c41f8d7a1942b6e16"
  }
  output {
    File selfSM = "${output_prefix}.selfSM"

    # I would like to do the following, however:
    # The object is read as a string
    # explicit string->float coercion via float(), as shown below, is supported by Cromwell
    # the interim value cannot be stored as a string and then assigned to a float. Variables intialized in output cannot be dereferenced in output.
    # Float contamination = float(read_object(${output_prefix} + ".selfSM").FREEMIX) / 0.75

    # In the interim, get the value from the python hack above:
    Float contamination = read_float(stdout())
  }
}

# PUBLIC #
workflow CheckContaminationWorkflow {
  File input_bam
  File input_bam_index
  File contamination_sites_ud
  File contamination_sites_bed
  File contamination_sites_mu
  File ref_fasta
  File ref_fasta_index
  String output_prefix
  Int disk_size
  Int preemptible_tries
  
  call CheckContamination {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_prefix = output_prefix,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }
}
