workflow picardTargetMapperWorkflow {
	String capture_type

	Array[String] ice_files = [
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/CRSP_ICE_hg19_wex_illumina_v1.no_X_Y_MT.bed",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/ICE_Combined_EEW_EWS_v1.pon",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/ICE_406_Normal_samples_PoN.vcf",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/whole_exome_illumina_coding_v1_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/ice_targets.bed",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/ice_rcs_eval.v1.pd250.spark.pon"
	]

	Array[String] agilent_files = [
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/whole_exome_agilent_1.1_refseq_plus_3_boosters_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets_germline_copy_number_variants_X_Y_removed.bed",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/BI_with_FFPE_Agilent_collapsed_v1.pon",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/gaf_20111020+broad_wex_1.1_hg19.bed",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/agilent_targets.bed",
	"gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/picard_target_mapper/param_files/agilent_rcs_eval_FFPE.v1.pd250.spark.pon"
	]

	Array[String] files_to_use = if capture_type == "ICE" then ice_files else agilent_files

	output {
		File contest_capture_targets_interval_list = files_to_use[0]
		File ReCapSeg_target_bed = files_to_use[1]
		File ReCapSeg_PON = files_to_use[2]
		File mutect_panel_of_normals_vcf_capture = files_to_use[3]
		File mutect_capture_targets_interval_list = files_to_use[4]
		File gatk4cnv_target_bed_capture = files_to_use[5]
		File gatk4cnv_pon_capture = files_to_use[6]
	}
}