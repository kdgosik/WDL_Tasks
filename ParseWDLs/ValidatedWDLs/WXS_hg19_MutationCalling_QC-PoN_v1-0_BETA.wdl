
task qc_prepare_task_oneBam_scatters {

	Array[String] bam_list
	File bam_list_fof=write_lines(bam_list)
	Array[String] bai_list
	File bai_list_fof=write_lines(bai_list)
	Int? scatterDiskGB
	Int defaultScatterDiskGB

	command <<<
		set -x
		#generate indices to index sizes and interval files
		#subtract 1 because indexes are zero-based
		cat ${bai_list_fof} |grep -Pn '.'|grep -Po '^\d+' | awk '{print $1-1}' >  indices.dat
		mv -vi ${bam_list_fof} bam_list.txt
		mv -vi ${bai_list_fof} bai_list.txt
		#1000 is the default for the scatter disk sizes ; cat over indices.dat so that the array is of the same length
		cat indices.dat | awk '{ print "${if defined(scatterDiskGB) then scatterDiskGB else defaultScatterDiskGB}"   } ' > scatter_sizes.dat ;
		>>>

	runtime {
		docker : "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
			disks: "local-disk 1 HDD"
			memory: "0.1 GB"
		}

	output {
		Array[Int] scatterDiskSizes=read_lines("scatter_sizes.dat")
		Array[Int] scatterIndices=read_lines("indices.dat")
		Array[String]  bams_list=read_lines("bam_list.txt")
		Array[String] bais_list=read_lines("bai_list.txt")
	}
}


task QC_Scatter_Task_OneBam {

	File bam
	File bamIdx
	String outname=basename(bam)+".rcl"
	File regionFile
	Int preemptible
	Int diskGB

	command <<<
		set -x

		#split by chromosome
		for CHROM in `seq 1 24` ; do 
			OUTPATH="normal.chr$CHROM.control.rcl" ;
			echo "Using OUTPATH=$OUTPATH FOR CHROM=$CHROM" ; 
			java -jar /usr/local/bin/RegionCovPerLane.jar ${bam} ${regionFile} $OUTPATH $CHROM
			cat $OUTPATH >> ${outname} ;
		done ;

		>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		disks: "local-disk ${diskGB} HDD"
		cpu: "1"
		memory: "4G"
		preemptible: "${preemptible}"
		}		

	output {
		File RCL="${outname}"
		}

	}


task make_qc_pon_from_bam_rcls {

	Array[File] make_qc_pon_from_bam_rcls
	String zipName

	command <<<
		set -x
		mkdir -v normal_db
		find . -iname "*.rcl" -exec mv -v {} normal_db   \;
		zip -r ${zipName}.zip normal_db
		>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		disks: "local-disk 20 HDD"
		}

	output {
		File zipDB="${zipName}.zip"
		}




	}




task merge_qc_rcl_pons {
	
	File thisPon
	File? otherPon
	String mergePonName

	command <<<

	set -x
	OTHER_PON_COUNT=`echo -ne "${otherPon}" | wc -c` ; 

	if [ "$OTHER_PON_COUNT" -eq "0" ] ;
	then 
		#not merging with existing PON
		echo "NOT performing PON merging...."
		mv -vi ${thisPon} ${mergePonName}.zip
	else
		echo "Performing PON merging with ${thisPon} and ${otherPon} ..." ;
		unzip ${thisPon}
		unzip ${otherPon}
		MERGE_DIR=${mergePonName}_mergedPON
		mkdir -v $MERGE_DIR
		find . -iname "*.rcl" -exec mv -v {} $MERGE_DIR \;
		zip -r ${mergePonName}.zip $MERGE_DIR
	fi ;

	>>>

	runtime {
		docker: "broadinstitute/broadmutationcalling_qc_beta@sha256:b2caf8864681b54b9c9825822be47f221ca577f84faa0b3240f7010712b9dfd3"
		disks: "local-disk 20 HDD"
		}

	output {
		File  mergeZipDB="${mergePonName}.zip"
		}

}




workflow qc_pon_make
	{
   	File regionFile
   	String set_id


	call qc_prepare_task_oneBam_scatters 

	scatter (idx in qc_prepare_task_oneBam_scatters.scatterIndices) {

			call QC_Scatter_Task_OneBam {
				input:
					bam=qc_prepare_task_oneBam_scatters.bams_list[idx],
					bamIdx=qc_prepare_task_oneBam_scatters.bais_list[idx],
					regionFile=regionFile,
					diskGB=qc_prepare_task_oneBam_scatters.scatterDiskSizes[idx]
				}
			}

	call make_qc_pon_from_bam_rcls  {
		input:
			make_qc_pon_from_bam_rcls=QC_Scatter_Task_OneBam.RCL,
			zipName=set_id
		}

	call merge_qc_rcl_pons {
		input:
			thisPon=make_qc_pon_from_bam_rcls.zipDB,
			mergePonName=set_id
		}


	output {
		make_qc_pon_from_bam_rcls.zipDB
		merge_qc_rcl_pons.mergeZipDB
	}	

    
    }

