
task token_make
	{

	File refTar
	File bam
	File bamIdx
	#use 2.1 multiplier because we need space to hold the tar and also everything once it's unpacked
	#use constant to hold output
	#add 10 at end as buffer (logs, wiggle room, etc)
	Float sizeOfTokenFile=3.3
	Float totSizeFloat=size(refTar,"G")*2.1+size(bam,"G")+size(bamIdx,"G")+sizeOfTokenFile+10
	String stringDiskSizeGB=sub(totSizeFloat,"\\..*","")
    Int diskGB=stringDiskSizeGB
    Int preemptible
	

	command <<<

		#increase verbosity and adjust error tolerance
		set -eux -o pipefail

		#unzip the ref tarbal
		DB_DIR=`tar -tvf ${refTar} |grep -P '/$'|grep -Po '\s\S+$'|tr -d " "` ;
		DB_DIR="$PWD/$DB_DIR" ; 
		tar -xvf ${refTar} ; 

		#OUTPUT token file
		TOKEN_OUT=`basename ${bam}`.token.dat ; 

		#scan for tokens
		java -Xmx4g -cp /xchip/cga_home/jhess/pileups_for_Keren/gpat.jar:/xchip/cga_home/jhess/pileups_for_Keren/htsjdk.jar org.broadinstitute.cga.tools.seq.GetAllPileupsAsTokens ${bam} $TOKEN_OUT $DB_DIR ;

		>>>

	runtime {
		preemptible: "${preemptible}"
		docker : "broadinstitute/blat_filtering:token_pon1"
        disks: "local-disk ${diskGB} HDD"
		}

	output
		{

		File token_file=basename(bam)+".token.dat"
		}


	}





task aggregateIntermediateTokensAndMakePon {

	File fof
	Array[File] tokenFiles=read_lines(fof)
	Int numTokenFiles=length(tokenFiles)
	Float sizeOfTokenFile=size(tokenFiles[0],"G")
	Float sizeOfPon=50.0
	String setName
	#have 10 as buffer here
	Float totSizeFloat=sizeOfTokenFile*numTokenFiles+sizeOfPon+10
	String stringDiskSizeGB=sub(totSizeFloat,"\\..*","")
	Int diskGB=stringDiskSizeGB
	Int preemptible

	command <<<

		#increase verbosity and adjust error tolerance
		set -eux -o pipefail

		#source matlab 
		source /matlab_source_file_2013a.sh 

		#make the FOF (file of files) for the token-maker
		cp -vf ${write_lines(tokenFiles)} "${setName}token_files_to_aggregate.fof"

		#simply display the contents of the FOF
		echo "FOF for Review"
		cat "${setName}token_files_to_aggregate.fof"

		#make the intermediate pon
		/xchip/cga_home/jhess/pileups_for_Keren/agg_tok_mcc/aggregate_tokens_files "${setName}token_files_to_aggregate.fof" "${setName}.final_summed_tokens.hist.bin"

		>>>

	runtime {
		docker : "broadinstitute/blat_filtering:token_pon2"
	    disks: "local-disk ${diskGB} HDD"
	    preemptible: "${preemptible}"
		}

	output {
		File fofFile="${setName}token_files_to_aggregate.fof" 
		File ponFile="${setName}.final_summed_tokens.hist.bin"
		}

	}



task intermediate_gather {
	
	Array[String] inGSURLs
	File inGSURLSFOF=write_lines(inGSURLs)
	Int intermediateScatterWidth
	Int preemptible

	command <<<

	#increase verbosity and adjust error tolerance
	set -eux -o pipefail

	#show input
	echo "Input:"
	grep -Pin '.' ${inGSURLSFOF} 

	#split into batches
	split --lines=${intermediateScatterWidth} ${inGSURLSFOF}  split_

	#show output
	head -v --lines=${intermediateScatterWidth} split_*

	>>>

	runtime {
		docker : "ubuntu:16.04"
	    disks: "local-disk 10 HDD"
	    preemptible: "${preemptible}"
		}

	output {
		Array[File] out_batches=glob("split_*")
		}

	}

task merge_pons {
	Array[File] pon_files
	Int numPonFiles=length(pon_files)
	Float sizeOfPonFile=size(numPonFiles[0],"G")
	String setName
	#have 10 as buffer here
	Float totSizeFloat=sizeOfPonFile*numPonFiles+sizeOfPonFile+10
	String stringDiskSizeGB=sub(totSizeFloat,"\\..*","")
	Int diskGB=stringDiskSizeGB
	Int preemptible

	command <<<
	#increase verbosity and adjust error tolerance
	set -eux -o pipefail

	#source matlab 
	source /matlab_source_file_2013a.sh 

	#make the FOF (file of files) for the pon-merger
	cp -vf ${write_lines(pon_files)} "${setName}_files_to_aggregate.fof"

	#show input
	cat "${setName}_files_to_aggregate.fof" | grep -Pin '.'

	#merge!!!
	/xchip/cga_home/jhess/pileups_for_Keren/agg_tok_mcc/sum_token_hists "${setName}_files_to_aggregate.fof" "${setName}.final_summed_tokens.hist.bin"


	>>>

	runtime {
		docker : "broadinstitute/blat_filtering:token_pon2"
	    disks: "local-disk ${diskGB} HDD"
	    preemptible: "${preemptible}"
		}

	output {
		File fofFile="${setName}_files_to_aggregate.fof" 
		File ponFile="${setName}.final_summed_tokens.hist.bin"
		}


	}


workflow token_wf
	{

	#bams, indices, and reference data
	Array[File] bamIndices
	Array[File] bams
	File refTar
	Array[Int] indices=range(length(bamIndices))
	Int preemptible_scatter
	Int preemptible_gather
	Int intermediateScatterWidth
	String setName

	#scatter
    scatter (index in indices)
    	{
        call token_make
        	{
            input:
            	#if lengths of the arrays don't match, there will be problems
				bam=bams[index],
				bamIdx=bamIndices[index],
                refTar=refTar,
                preemptible=preemptible_scatter
            }
        }

    #intermediate gather
	call intermediate_gather {
		input:
			inGSURLs=token_make.token_file,
			intermediateScatterWidth=intermediateScatterWidth,
			preemptible=preemptible_gather
		}
		

	#make a bunch of PONs (intermediate for batch)
	scatter(ifof in intermediate_gather.out_batches)
		{
	    #build the intermediate PONs
	    call aggregateIntermediateTokensAndMakePon {
	    	input:
	    		fof=ifof,
	    		preemptible=preemptible_gather,
	    		setName=setName
	    	}
		}

	#merge the intermediate PONs
	call merge_pons {
		input:
			pon_files=aggregateIntermediateTokensAndMakePon.ponFile,
			setName=setName,
			preemptible=preemptible_gather
		}


}
