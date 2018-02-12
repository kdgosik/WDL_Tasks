task aggregate_clinical {
	Array[String] sample_ids
	Array[File] clinical_primary
	Array[File] clinical_biospecimen
	String outPrefix

	command {
	python <<CODE
NULL_SENTINEL = "GDAC_FC_NULL"

# Parse array inputs
keys = '${sep="," sample_ids}'.split(",")
prim_values = '${sep="," clinical_primary}'.split(",")
bio_values = '${sep="," clinical_biospecimen}'.split(",")

# Intermediate tsv files
prim_tsv = open("clin_datapaths.tsv", "w")
bio_tsv = open("bio_datapaths.tsv", "w")

# Put in datapaths for clinical and bio values
for i in range(len(keys)):
	if not prim_values[i].endswith(NULL_SENTINEL):
		prim_tsv.write(keys[i] + "\t" + prim_values[i] + "\n")

	if not bio_values[i].endswith(NULL_SENTINEL):
		bio_tsv.write(keys[i] + "\t" + bio_values[i] + "\n")

bio_tsv.close()
prim_tsv.close()
CODE

python /src/clinical_data_merger.py --clinical-inputs "clin_datapaths.tsv" --biospecimen-inputs "bio_datapaths.tsv" ${outPrefix}

	}

	runtime {
		docker: "broadgdac/aggregate_clinical:73"
	}

	meta {
		author : "Tim DeFreitas"
		email : "timdef@broadinstitute.org"
	}

	output {
		File clinical__merged="${outPrefix}.clin.merged.txt"
	}

}

workflow aggregate_clinical_workflow {
	call aggregate_clinical
}
