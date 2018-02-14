#This WDL was adapted from https://github.com/NCIP/Trinity_CTAT/tree/master/docker/

#This task requires unpacking a directory that may have subdirectories. Future tasks will need the 
#file structure to remain intact. Therefore it moves all of the unpacked files into a single directory
#and renames them so that their order in glob will be deterministic. It also outputs an Array of strings
#which are the file paths for each of these moved and renamed files. This allows future tasks to move
#them back to where they were originally.
task UnpackRefDir {
    File genome_lib_dir

    command {
        mkdir refDir
        mkdir newRefDir
        tar xvf ${genome_lib_dir} -C refDir 

        i=1
        for f in `find refDir -type f`; do
            echo $f >> list_of_files.txt
            mv $f newRefDir/$(printf "%0.3d" $i)

            ((i++))
        done
    }
    
    output {
        Array[File] reference_dir_files = glob("newRefDir/*")
        Array[String] reference_dir_paths = read_lines("list_of_files.txt")
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
        disks: "local-disk 100 HDD"
        memory: "1200 MB"
    }
}

#Future tasks will take the unpacked files and move them back to where they were originally. It goes through each of 
#the paths and has to iterate through all of the files (because bash's ${} is not possible in a WDL command, so 
#the arrays can't be accessed by index explicitly). When a path and file match it moves the file to the destination 
#the path specifies.
task StarFusion {
    File left_fq
    File right_fq
    Array[File] reference_dir_files
    Array[String] reference_dir_paths

    command {
        i=0
        for path in "${sep='" "' reference_dir_paths}"; do
            j=0
            for refPath in "${sep='" "' reference_dir_files}"; do
                refFile=$refPath
                if [ $j -eq $i ]; then 
                    break 
                fi
                ((j++))
            done
            mkdir -p `dirname $path`
            mv $refFile $path
            ((i++))
        done

        genomeDir=$(ls -d refDir/*/)

        /usr/local/src/STAR-Fusion-v1.0.0/STAR-Fusion \
            --left_fq ${left_fq} \
            --right_fq ${right_fq} \
            --output_dir outDir \
            --genome_lib_dir $genomeDir
    }

    output {
        File out_file = "outDir/star-fusion.fusion_candidates.final.abridged.FFPM"
        File chimeric_file = "outDir/Chimeric.out.junction"
        Array[File] all_out_files = glob("outDir/*")
    }

    runtime {
        docker: "trinityctat/ctatfusion:1.0.0"
        disks: "local-disk 100 SSD"
        memory: "45G"
        cpu: "4"
    }
}

task FusionInspector {
    File star_fusion_out
    Array[File] reference_dir_files
    Array[String] reference_dir_paths
    File left_fq
    File right_fq

    command {
        i=0
        for path in "${sep='" "' reference_dir_paths}"; do
            j=0
            for refPath in "${sep='" "' reference_dir_files}"; do
                refFile=$refPath
                if [ $j -eq $i ]; then 
                    break 
                fi
                ((j++))
            done
            mkdir -p `dirname $path`
            mv $refFile $path
            ((i++))
        done

        genomeDir=$(ls -d refDir/*/)

        /usr/local/src/FusionInspector-v1.0.1/FusionInspector \
            --fusions ${star_fusion_out} \
            --genome_lib $genomeDir \
            --left_fq ${left_fq} \
            --right_fq ${right_fq} \
            --out_dir outDir \
            --out_prefix finspector \
            --align_utils STAR --no_cleanup
    }

    output {
        Array[File] out_file = glob("outDir/*")
        File out_fusions = "outDir/finspector.fusion_predictions.final.abridged.FFPM"
    }

    runtime {
        docker: "trinityctat/ctatfusion:1.0.0"
        disks: "local-disk 100 SSD"
        memory: "45G"
        cpu: "4"
    }
}

task RevertSam {
  File input_bam
  String output_name

  command {
    java -Xmx1000m -jar /usr/gitc/picard.jar \
    RevertSam \
    INPUT=${input_bam} \
    OUTPUT_BY_READGROUP=false \
    VALIDATION_STRINGENCY=LENIENT \
    ATTRIBUTE_TO_CLEAR=FT \
    SORT_ORDER=queryname \
    OUTPUT=${output_name}.bam 
  }
  output {
    File unmapped_bam = "${output_name}.bam"
  }

  runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
        disks: "local-disk 100 HDD"
        memory: "1200 MB"
    }
}

task SamToFastq {
  File input_bam
  String sample

  command {
    java -jar /usr/gitc/picard.jar \
    SamToFastq I=${input_bam} \
    F=${sample}.1.fastq F2=${sample}.2.fastq \
    INTERLEAVE=false NON_PF=true \
    CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2
  }
  output {
    Array[File] fastqs = ["${sample}.1.fastq", "${sample}.2.fastq"]
  }

  runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
        disks: "local-disk 100 HDD"
        memory: "2G"
    }
}

#FusionFilter is what was used to generate the reference used in this pipeline.
workflow trinity {
    File genome_lib
    File input_bam
    String sample

    call UnpackRefDir {
        input: 
            genome_lib_dir = genome_lib
    }

    call StarFusion {
        input: 
            left_fq = SamToFastq.fastqs[0],
            right_fq = SamToFastq.fastqs[1],
            reference_dir_files = UnpackRefDir.reference_dir_files,
            reference_dir_paths = UnpackRefDir.reference_dir_paths
    }

    call RevertSam {
        input: input_bam = input_bam,
            output_name = sample+".unmapped"
    }

    call SamToFastq {
        input: input_bam = RevertSam.unmapped_bam,
            sample = sample
    }

    call FusionInspector {
        input:
            star_fusion_out = StarFusion.out_file,
            left_fq = SamToFastq.fastqs[0],
            right_fq = SamToFastq.fastqs[1],
            reference_dir_files = UnpackRefDir.reference_dir_files,
            reference_dir_paths = UnpackRefDir.reference_dir_paths
    }
}
