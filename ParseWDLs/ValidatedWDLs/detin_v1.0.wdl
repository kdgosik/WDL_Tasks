task detin_release {
  File call_stats_file
  File copy_number_data
  File tumor_hets
  File normal_hets
  File exac_pickle
  File indels_file
  String MuTect2_or_Strelka
  String TiN_prior
  String Mutation_prior
  String pair_id
  String preemptible_limit
  String output_disk_gb
  String boot_disk_gb
  String ram_gb
  String cpu_cores
  String release_version
  
      command <<<

      #increase verbosity
      set -x
      pwd
      ls -latrh 
      # generate null outputs
      echo null > ${pair_id}.TiN_estimate.txt
      echo null > ${pair_id}.TiN_estimate_CI.txt
      echo null > ${pair_id}.deTiN_SSNVs.txt
      echo null > ${pair_id}_SSNVs_plot.png
      echo null > ${pair_id}_KmeansEval_plot.png
      echo null > ${pair_id}_TiN_models_plot.png
      echo null > ${pair_id}_KmeansEval_scatter_plot.png 
      
      git clone https://github.com/broadinstitute/deTiN.git
      git -C deTiN/ checkout tags/${release_version}
      python deTiN/deTiN.py --mutation_data_path ${call_stats_file} --cn_data_path ${copy_number_data} --tumor_het_data ${tumor_hets} --normal_het_data ${normal_hets} --exac_data_path ${exac_pickle} --output_name ${pair_id} --TiN_prior ${TiN_prior} --mutation_prior ${Mutation_prior} --indel_data_path ${indels_file} --indel_data_type ${MuTect2_or_Strelka}
      >>>
      
        output {
        String TiN=read_string("${pair_id}.TiN_estimate.txt")
		String TiN_CI=read_string("${pair_id}.TiN_estimate_CI.txt")
      	File deTiN_call_stats="${pair_id}.deTiN_SSNVs.txt"
   		File deTiN_SSNVs_plot = "${pair_id}_SSNVs_plot.png"
        File deTiN_aSCNA_kmeans_RSS_plot = "${pair_id}_KmeansEval_plot.png"
        File deTiN_aSCNA_scatter_plot = "${pair_id}_KmeansEval_scatter_plot.png"
        File deTiN_TiN_modes_plot = "${pair_id}_TiN_models_plot.png"
  }
  runtime {
    memory: "${ram_gb}GB"
    cpu: "${cpu_cores}"
    disks: "local-disk ${output_disk_gb} HDD"
    bootDiskSizeGb: "${boot_disk_gb}"
    preemptible: "${preemptible_limit}"
    docker: "amarotaylor/detin_python:latest"

  }
 meta {
        author : "Amaro Taylor-Weiner"
        email : "amaro@broadinstitute.org"
    }
}
workflow detin_release_wkflow{

  call detin_release 
    
}
 

