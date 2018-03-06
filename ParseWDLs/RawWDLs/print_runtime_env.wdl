task print_runtime_env_task {
	Int preemptible_tries
	Int disk_size
	Int memory

	command <<<
		echo -n "Number of Processors: " > runtime_info.txt; 
		nproc >> runtime_info.txt
		cat /proc/meminfo | grep MemTotal >> runtime_info.txt
		df >> runtime_info.txt
	>>>

	output {
		File outFile="runtime_info.txt"
		}

	runtime {
		docker: "ubuntu:14.04.4"
		memory: memory + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
		}
}

workflow print_runtime_env_workflow {
	Int preemptible_tries
	Int disk_size
	Int memory

	call print_runtime_env_task {
		input:
			memory = memory,
			disk_size = disk_size,
			preemptible_tries = preemptible_tries
	}

}