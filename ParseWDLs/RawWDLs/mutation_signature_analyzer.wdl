task mutation_signature_analyzer {
	Boolean package
	String null_file = "gs://broad-institute-gdac/GDAC_FC_NULL"
	String libdir = "/src"
	Int n_iter
	Int Kcol
	String tol
	Int a0
	String tumor_type
	String hyper
	File maffile
	String genomeRefFolder
	Array[File] genomeRefFiles = [genomeRefFolder + "chr1.txt", genomeRefFolder + "chr2.txt", genomeRefFolder + "chr3.txt", genomeRefFolder + "chr4.txt", genomeRefFolder + "chr5.txt", genomeRefFolder + "chr6.txt", genomeRefFolder + "chr7.txt", genomeRefFolder + "chr8.txt", genomeRefFolder + "chr9.txt", genomeRefFolder + "chr10.txt", genomeRefFolder + "chr11.txt", genomeRefFolder + "chr12.txt", genomeRefFolder + "chr13.txt", genomeRefFolder + "chr14.txt", genomeRefFolder + "chr15.txt", genomeRefFolder + "chr16.txt", genomeRefFolder + "chr17.txt", genomeRefFolder + "chr18.txt", genomeRefFolder + "chr19.txt", genomeRefFolder + "chr20.txt", genomeRefFolder + "chr21.txt", genomeRefFolder + "chr22.txt", genomeRefFolder + "chrX.txt", genomeRefFolder + "chrY.txt"]
	String prior

	command {
		set -euo pipefail
		genomeRefFolder=`dirname ${genomeRefFiles[0]}`
		Rscript /src/SignatureAnalyzer.Broad.R \
		--n.iter ${n_iter} \
		--Kcol ${Kcol} \
		--tol ${tol} \
		--a0 ${a0} \
		--tumor.type ${tumor_type} \
		--hyper ${hyper} \
		--maffile ${maffile} \
		--genomeRefFolder $genomeRefFolder \
		--prior ${prior} 


		if ${package}
			then zip -r mutation_signature_analyzer . -x \
				"fc-[a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9]/*" \
				lost+found/\* \
				broad-institute-gdac/\* \
				"tmp.[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]/*" \
				exec.sh
		fi

	}

	output {
		File mutation_signature_analyzer_pkg="${if package then 'mutation_signature_analyzer.zip' else null_file}"

	}

	runtime {
		docker : "broadgdac/mutation_signature_analyzer:1.2"
	}

	meta {
		author : "Vicky Horst"
		email : "gdac@broadinstitute.org"
	}

}

workflow mutation_signature_analyzer_workflow {
	call mutation_signature_analyzer
	meta {
		version : "1.0"

	}
}
