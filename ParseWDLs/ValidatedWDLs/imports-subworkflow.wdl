task TheTask {
  command { foo }

  runtime {
    docker: "python:2.7"
  }

  output {
	String output1 = "this is version 1."
  }
}

workflow TheWorkflow {
  call TheTask
}