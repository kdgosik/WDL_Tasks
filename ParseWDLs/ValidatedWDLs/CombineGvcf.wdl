
task TsvToArrayOfArrays {
    File tsv

    command {

    }

    runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }

    output {
        Array[Array[String]] output_matrix = read_tsv("${tsv}")
    }
}

task IntToIntArray {
    Int number

    command <<<
        python <<CODE
        for i in range(${number}):
          print(i)
        CODE
      >>>

    runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }

    output {
        Array[Int] array = read_lines(stdout())
    }
}

task FileToStringArray {
    File input_file

    command {

    }

    runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }

    output {
        Array[String] array = read_lines(input_file)
    }
}


task TileDBCombineFoFn {
  String output_gvcf_filename
  String partition

  Array[File] files

  File ref_fasta
  File ref_fasta_index
  File vid_mapping_file
  File vcf_header
  File sample_name_map

  Int size_per_column_partition
  Int disk_size

  command <<<
  set -e
  # this is here to deal with the JES bug where commands may be run twice
  rm -rf gvcfs
  mkdir gvcfs

  mv ${sep=" " files} gvcfs

  python <<CODE
from collections import OrderedDict
import json
import os
import re
import subprocess
import sys

w = open("/usr/gitc/callset.json", "w")
array_of_gvcfs = subprocess.Popen(['ls', 'gvcfs'], stdout=subprocess.PIPE).communicate()[0].strip().split("\n")
array_length = len(array_of_gvcfs)

with open("${sample_name_map}") as sample_map:
  file_to_sample_map = dict((k, v) for k, v in (map(str, line.strip().split("\t")) for line in sample_map))

for gvcf in array_of_gvcfs:
  filename = gvcf
  subprocess.call(['/usr/gitc/tabix', "gvcfs/" + filename])

vcf_counter = 0
callsets = OrderedDict()
for gvcf in array_of_gvcfs:
  filename = gvcf
  mapping_filename = re.sub("\.\d+\.g\.vcf", ".g.vcf", filename)
  vcf_data = OrderedDict()
  vcf_data["row_idx"] = vcf_counter
  # this is the samples index in the gvcf, because we dont generate multi-sample gvcfs this is always 0
  vcf_data["idx_in_file"] = 0
  vcf_data["filename"] = "gvcfs/" + filename
  sample_name = file_to_sample_map[mapping_filename]
  callsets[sample_name] = vcf_data
  vcf_counter += 1

final_json = OrderedDict()
final_json["callsets"] = callsets

w.write(json.JSONEncoder().encode(final_json))
w.close()

w = open("/usr/gitc/loader.json", "w")

column_partitions = [ ${partition} ]

final_json = OrderedDict()
final_json["column_partitions"] = column_partitions
final_json["callset_mapping_file"] = "/usr/gitc/callset.json"
# this value is currently 30000 * number of gvcfs in the callset
final_json["size_per_column_partition"] = ${size_per_column_partition}
final_json["treat_deletions_as_intervals"] = True
final_json["vcf_header_filename"] = "${vcf_header}"
final_json["vid_mapping_file"] = "${vid_mapping_file}"
final_json["reference_genome"] = "${ref_fasta}"
final_json["num_parallel_vcf_files"] = 1
final_json["do_ping_pong_buffering"] = True
# if this is set to True, then the number of cpus for this task should be increased to 3
final_json["offload_vcf_output_processing"] = True
final_json["discard_vcf_index"] = True
final_json["produce_combined_vcf"] = True
# zipped output
final_json["vcf_output_format"] = "z"
final_json["vcf_output_filename"] = "${output_gvcf_filename}"

w.write(json.JSONEncoder().encode(final_json))
w.close()
CODE

    (/usr/bin/time -v /usr/gitc/vcf2tiledb /usr/gitc/loader.json)

    >>>
  runtime {
    docker: "broadinstitute/tile-db:1.1.3"
    memory: "15 GB"
    cpu: "3"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_gvcf = "${output_gvcf_filename}"
  }
}

task ArrayToFile {
    Array[String] array

    command {
    }

    runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }

    output {
        File write_output = write_lines(array)
    }
}

workflow JointGenotyping {
  File ref_fasta
  File ref_fasta_index

  # this maps the genomic position based on contigs to the tiledb interpretation of the same position
  # chr 3 pos 300 for tiledb would be ((size of chr1 + size of chr2 + 300) - 1)  (0-based)
  File vid_mapping_file
  File vcf_header

  Array[String] gvcf_string_array

  Int size_per_column_partition

  Int medium_disk

  File sample_map

  File gvcf_array

  Int num_partitions
  File tiledb_partition_list_file

  call IntToIntArray {
      input:
        number = num_partitions
  }

  scatter (idx in IntToIntArray.array) {

    call TileDBCombineFoFn {
      input:
        output_gvcf_filename = "output.vcf.gz",
        files = read_tsv(gvcf_array)[idx],
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        partition = read_lines(tiledb_partition_list_file)[idx],
        size_per_column_partition = size_per_column_partition,
        vid_mapping_file = vid_mapping_file,
        vcf_header = vcf_header,
        sample_name_map = sample_map,
        disk_size = medium_disk
    }
  }

  call ArrayToFile {
    input:
      array = TileDBCombineFoFn.output_gvcf
  }

  output {
    ArrayToFile.*
  }
}
