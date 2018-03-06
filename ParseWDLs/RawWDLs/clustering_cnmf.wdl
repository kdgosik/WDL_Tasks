task cluster_cnmf {
	File expfile
	Int k_int
	Int k_final
	String outputprefix

	runtime { docker: "broadgdac/cluster_cnmf:51"}

	command {
		/R/RunR.sh -f main /src/GDAC_CNMF.R -m${expfile} -u${k_int} -v${k_final} -o${outputprefix} -s/src
	}

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

task cluster_cnmf_select {
	String measure
	File inputexp
	String outputprefix 
	File inputallexp
	File clumembership
	File cophenetic

	command {
		/R/RunR.sh -f main /src/select_best_cluster_cnmf_v1.R -m${measure} -u${inputexp} -v${outputprefix} -w${clumembership} -p${inputallexp} -c${cophenetic}
	}

	runtime { docker: "broadgdac/cluster_cnmf_select:62"}

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

task preprocess_expression_for_cluster {
	File expfile
	String? selectedgenes ##Optional in FH, but not really optional
	String outputprefix

	command {
		/R/RunR.sh -f main /src/Topgenes_v1.R -s/src -m${expfile} -u${selectedgenes} -o${outputprefix}
	}

	runtime { 
		docker: "broadgdac/preprocess_expression_for_cluster:47"
	}

	output {
		File expression_file="${outputprefix}.expclu.gct"
	}

}

task report_cluster_cnmf {
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
	String report_type

	command {
		/R/RunR.sh -f main /src/ReportGenerator.R -hNULL -o${expdata} -v${kclus} -s${bestclu} -u${allcluster} -w${markers} -a${file_gif_2} -b${file_gif_3} -c${file_gif_4} -d${file_gif_5} -e${file_gif_6} -f${file_gif_7} -g${file_gif_8} -q${markersP} -r${heatmap} -t${heatmapall} -y/src -z${report_type}

	}

	runtime { docker: "broadgdac/report_cluster_cnmf:35"}

	output {
		File report="nozzle.html"
		File rData="nozzle.RData"

	}

}

workflow clustering_cnmf {
	call preprocess_expression_for_cluster
	call cluster_cnmf {
		input: expfile=preprocess_expression_for_cluster.expression_file
	}
	call cluster_cnmf_select {
		input: inputexp=preprocess_expression_for_cluster.expression_file, clumembership=cluster_cnmf.membershipfile, cophenetic=cluster_cnmf.coefficientfile
	}

	call report_cluster_cnmf {
		input: expdata=preprocess_expression_for_cluster.expression_file, kclus=cluster_cnmf_select.sigkclus, bestclu=cluster_cnmf_select.bestcluster,
		allcluster=cluster_cnmf.membershipfile, markers=cluster_cnmf_select.markerfile,
		file_gif_2=cluster_cnmf.plotk2, file_gif_3=cluster_cnmf.plotk3, file_gif_4=cluster_cnmf.plotk4,
		file_gif_5=cluster_cnmf.plotk5, file_gif_6=cluster_cnmf.plotk6, file_gif_7=cluster_cnmf.plotk7,
		file_gif_8=cluster_cnmf.plotk8, markersP=cluster_cnmf_select.selectmarkers, heatmap=cluster_cnmf_select.geneheatmap,
		heatmapall=cluster_cnmf_select.geneheatmapall
	}

}

