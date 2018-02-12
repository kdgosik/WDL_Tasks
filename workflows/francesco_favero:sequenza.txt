task SequenzaTask {
  String sampleName
  File normalBam
  File tumorBam
  File ReferenceFastaGz
  File ReferenceGcWig
  command <<<
  set -x
  sequenza-pipeline \
      --sample-id ${sampleName} \
      --normal-bam ${normalBam} \
      --tumor-bam ${tumorBam} \
      --reference-gz ${ReferenceFastaGz} \
      --gc_wig ${ReferenceGcWig}
  >>>

  runtime {
    docker: "sequenza/sequenza"
    memory: "24 GB"
    disks: "local-disk 100 SSD"
    }

  output {
    File logs = "${sampleName}_logs.tar.gz"
    File binSeqz = "${sampleName}_seqz_bin.tar.gz"
    File fullSeqz = "${sampleName}_parts_seqz.tar.gz"
    File scnaRes = "${sampleName}_sequenza.tar.gz"
  }
}

workflow SequenzaWorkflow {
  String sampleName
  File normalBam
  File tumorBam
  File ReferenceFastaGz
  File ReferenceGcWig

  call SequenzaTask {
    input: sampleName=sampleName,
           normalBam=normalBam, tumorBam=tumorBam,
           ReferenceFastaGz=ReferenceFastaGz,
           ReferenceGcWig=ReferenceGcWig
  }
}