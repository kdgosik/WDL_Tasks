task top_genes_for_cluster {
	File expfile
	String selectedgenes ##Optional in FH, but not really optional
	String outputprefix

	command <<<
		/R/RunR.sh -f main /src/Topgenes_v1.R -s./ -m${expfile} -u${selectedgenes} -o${outputprefix}
	>>>

	runtime {
		docker: "broadgdac/top_genes_for_cluster:46"
    	memory: "12 GB" 
    	disks: "local-disk 100 SSD"
		}


	output {
		File expression_file="${outputprefix}.expclu.gct"
		String outoutputprefix="${outputprefix}"
	}
}



task nmf_consensus_clustering {
	File expfile
	Int k_int
	Int k_final
	String outputprefix

	runtime { 
		docker: "broadgdac/nmf_consensus_clustering:51"
    		memory: "12 GB" 
    		disks: "local-disk 100 SSD"
	}

	command <<<
		/R/RunR.sh -f main /src/GDAC_CNMF.R -m${expfile} -u${k_int} -v${k_final} -o${outputprefix} -s/src
	>>>

	output {
		File plotk2="${outputprefix}.consensus.plot.k2.png"
		File plotk3="${outputprefix}.consensus.plot.k3.png"
		File plotk4="${outputprefix}.consensus.plot.k4.png"
		File plotk5="${outputprefix}.consensus.plot.k5.png"
		File plotk6="${outputprefix}.consensus.plot.k6.png"
		File plotk7="${outputprefix}.consensus.plot.k7.png"
		File plotk8="${outputprefix}.consensus.plot.k8.png"

		File coefficientfile="${outputprefix}.cophenetic.coefficient.txt"
		File membershipfile="${outputprefix}.membership.txt"

	}
}

task cnmf_select_cluster {
	String measure
	File inputexp
	String outputprefix 
	File inputallexp
	File clumembership
	File cophenetic

	command <<<
		/R/RunR.sh -f main /src/select_best_cluster_cnmf_v1.R -m${measure} -u${inputexp} \
		-v${outputprefix} -w${clumembership} -p${inputallexp} -c${cophenetic}
	>>>

	runtime { 
		docker: "broadgdac/cnmf_select_cluster:62"
    	memory: "12 GB" 
    	disks: "local-disk 100 SSD"
	}

	output {
		File sigkclus="${outputprefix}.silfig.png"
		File bestcluster="${outputprefix}.bestclus.txt"
		File markerfile="${outputprefix}.subclassmarkers.txt"
		File cormatrixpng="${outputprefix}.cormatrix.png"
		File selectmarkers="${outputprefix}.selectmarker.txt"
		File geneheatmap="${outputprefix}.geneheatmap.png"
		File geneheatmapall="${outputprefix}.geneheatmaptopgenes.png"
		
	}

}




task cnmf_report {
	File expdata
	File kclus
	File bestclu
	File allcluster
	File markers
	File file_gif_2
	File file_gif_3
	File file_gif_4
	File file_gif_5
	File file_gif_6
	File file_gif_7
	File file_gif_8
	File markersP
	File heatmap
	File heatmapall
	String report
	String outputprefix

	command <<<

		#increase verbosity
		set -x

		/R/RunR.sh -f main /src/ReportGenerator.R -hNULL -o${expdata} -v${kclus} -s${bestclu} -u${allcluster} \
		-w${markers} -a${file_gif_2} -b${file_gif_3} -c${file_gif_4} -d${file_gif_5} -e${file_gif_6} \
		-f${file_gif_7} -g${file_gif_8} -q${markersP} -r${heatmap} -t${heatmapall} -y/src -z${report}

		#bring in the PNGs to that they can be de-localized
		find `pwd` -iname "*.png" -exec cp -vf {} . \;

	>>>

	runtime { 
		docker: "broadgdac/report_cluster_cnmf:35"
		memory: "12 GB"
    	disks: "local-disk 100 SSD"
	}

	output {
		File nozzle_report="nozzle.html"
		File rData="nozzle.RData"
		File heatmap_out="${outputprefix}.geneheatmap.png"
		File heatmapall_out="${outputprefix}.geneheatmaptopgenes.png"
		File plotk2="${outputprefix}.silfig.png"
		File kclus_out="${outputprefix}.consensus.plot.k2.png"

	}

}



workflow clustering_cnmf {
	call top_genes_for_cluster

	call nmf_consensus_clustering {
		input: expfile=top_genes_for_cluster.expression_file,
		outputprefix=top_genes_for_cluster.outoutputprefix
	}

	call cnmf_select_cluster {
		input: inputexp=top_genes_for_cluster.expression_file, 
		clumembership=nmf_consensus_clustering.membershipfile, 
		cophenetic=nmf_consensus_clustering.coefficientfile, 
		outputprefix=top_genes_for_cluster.outoutputprefix
	}

	call cnmf_report {
		input: expdata=top_genes_for_cluster.expression_file, 
		kclus=cnmf_select_cluster.sigkclus, 
		bestclu=cnmf_select_cluster.bestcluster,
		allcluster=nmf_consensus_clustering.membershipfile, 
		markers=cnmf_select_cluster.markerfile,
		file_gif_2=nmf_consensus_clustering.plotk2, 
		file_gif_3=nmf_consensus_clustering.plotk3, 
		file_gif_4=nmf_consensus_clustering.plotk4,
		file_gif_5=nmf_consensus_clustering.plotk5, 
		file_gif_6=nmf_consensus_clustering.plotk6, 
		file_gif_7=nmf_consensus_clustering.plotk7,
		file_gif_8=nmf_consensus_clustering.plotk8, 
		markersP=cnmf_select_cluster.selectmarkers, 
		heatmap=cnmf_select_cluster.geneheatmap,
		heatmapall=cnmf_select_cluster.geneheatmapall,
		outputprefix=top_genes_for_cluster.outoutputprefix
	}


}