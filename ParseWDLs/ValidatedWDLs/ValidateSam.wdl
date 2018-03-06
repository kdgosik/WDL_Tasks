task errors {
	File input_bam
	String output_name
	String mem_size
    Int disk_size

	command <<<
	/usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java -Xmx2g -jar /opt/picard-tools/picard.jar ValidateSamFile \
		I=${input_bam} \
		OUTPUT=${output_name}.errors.list \
		MODE=VERBOSE \
		IGNORE_WARNINGS=true
	>>>

	runtime {
   	docker: "elcinchu27/docker-picard"
    memory: mem_size
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    }
    output {
      File errors_list="${output_name}.errors.list"
    }
}

task warnings {
	File input_bam
	String output_name
	String mem_size
    Int disk_size

	command <<<
	/usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java -Xmx2g -jar /opt/picard-tools/picard.jar ValidateSamFile \
		I=${input_bam} \
		OUTPUT=${output_name}.warnings_errors.list \
		MODE=VERBOSE
	>>>

	runtime {
   	docker: "elcinchu27/docker-picard"
    memory: mem_size
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    }
    output {
      File warnings_errors="${output_name}.warnings_errors.list"
    }
}

workflow Validation{
	String mem_size
    Int disk_size
    File input_bam
	String output_name

	call errors {
		input:
			input_bam=input_bam,
			output_name=output_name,
			mem_size=mem_size,
      		disk_size=disk_size
	}

	call warnings {
		input:
			input_bam=input_bam,
			output_name=output_name,
			mem_size=mem_size,
      		disk_size=disk_size
	}

	output {
		File errors_list=errors.errors_list
		File warnings_errors=warnings.warnings_errors
	}

}
