task fof_usage_task {
    File fof
    Array[File] my_files=read_lines(fof)

    command {
    #do stuff with arrays below
    #....
    }
    runtime {
        docker : "ubuntu:16.04"
        disks: "local-disk 50 HDD"
        memory: "2 GB"
    }
    output {
	Array[File] array_of_files = my_files
     }  
}

workflow fof_usage_wf {
    File file_of_files
    call fof_usage_task {
	input:
    fof = file_of_files
    }
    
    output {
	Array[File] array_output = fof_usage_task.array_of_files
    }
}