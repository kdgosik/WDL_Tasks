task HistXtractor {

	File inputFile
    Int preemptible
    Int memorySize
	String tempFileName = sub(inputFile, "\\.svs$", ".seg.txt")
	String outputFileName = sub(tempFileName, ".*/", "")
	command <<<
		/code/run_HistXtractor.sh /usr/local/MATLAB/MATLAB_Runtime/v901/ ${inputFile} 20 4096 ./ true true 
        
		>>>

	output {
		#This is the cell data output file
		File segFile = "${outputFileName}"
	}

	runtime {
    	preemptible: "${preemptible}"
		docker: "broadinstitute/histxtractor:1.0"
		memory: "${memorySize} GB"
		disks: "local-disk 20 HDD"
	}

}



workflow processorWorkflow {
	File inputFile
    Int preemptible
    Int memorySize

	call HistXtractor {
		input:
			inputFile=inputFile,
            preemptible=preemptible,
            memorySize=memorySize
	}
	
	output {
		HistXtractor.segFile
	}
}