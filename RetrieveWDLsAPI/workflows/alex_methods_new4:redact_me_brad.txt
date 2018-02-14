task RevCompTask 
	{
	File inDNA
	command
		{
		rev_comper.py `cat ${inDNA}`
		}
	output
		{
		String outRC=read_string(stdout())
		}
 	meta {author : "Eddie Salinas"}
	runtime { docker: "eddiebroad/biopython" }
	}
	
workflow RevCompWorkflow
	{
	File inDNA
	call RevCompTask {
		input:inDNA=inDNA
		}
	output {
		RevCompTask.outRC
		}
	}