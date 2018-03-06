

workflow QCAndPrePhase {
	# sample plink files

	File dataset_bim
	File dataset_bed
	File dataset_fam

#	- dbsnp reference vcf (referenced in this documentation: imputation_datafiles/dbSNP)
	
	File dbSnp_vcf
	File dbSnp_vcf_index

#- human reference fasta, indexed (same used in your panel calling) (imputation_datafiles/fasta)
	File reference_fasta
	File reference_dict
	File reference_index

#- plink --freq output for all variants in panel(s) for all chromosomes in one single text file (referenced here as: imputation_datafiles/refPanel/WGS_WES_union_reference.frq). We did this for both panels separately and then created an union file of these two, prioritizing WGS frequency for overlapping sites.
	File plink_frq

#- genetic maps in i) eagle2 format, ii) impute2 format (imputation_datafiles/geneticMaps)
	File genetic_maps_eagle

	Int dbSnpSize=sub(size(dbSnp_vcf,"GB"), "\\..*", "")
	Int bedSize=sub(size(dataset_bed,"GB") * 1.3 , "\\..*", "") 

	Int extractRsIdFromVCF_size=dbSnpSize + 100

	call ExtractRsIdFromVCF{
		input:
			vcf=dbSnp_vcf,
			disk_size=extractRsIdFromVCF_size
	}

	call SubTask as dbRsIDsSize_string{input: string=size(ExtractRsIdFromVCF.rsIDs,"GB"), pattern= "\\..*", replacement=""}
	Int dbRsIDsSize=dbRsIDsSize_string.out
	
 	Int checkBuild_size=dbRsIDsSize + bedSize + 100

 	## instead of "checking" we need to compare the number of mismatches against several dbSNPs (from different builds) and see if the one corresponding 
 	## to hg37 is the one with the smallest number of mismatches...this will need to take chr1/1 differences into account....

 	call CheckBuild {
 		input:
 			dataset_bim=dataset_bim,
 			dataset_bed=dataset_bed,
 			dataset_fam=dataset_fam,
 			dbSnp_vcf  =dbSnp_vcf,
 			dbSnp_RSIDs=ExtractRsIdFromVCF.rsIDs,
 			disk_size = checkBuild_size 
 	}

 	if ( !CheckBuild.OK ) {

		Int liftOverPlinkFiles_size = 2*bedSize + dbRsIDsSize + 30
		call LiftOverPlinkFiles{
			input:
				dataset_bim=dataset_bim,
				dataset_bed=dataset_bed,
				dataset_fam=dataset_fam,
				dbSnp_RSIDs=ExtractRsIdFromVCF.rsIDs,
				update_map =ExtractRsIdFromVCF.update_map,
				disk_size  =liftOverPlinkFiles_size
		}
	}

	File dataset_bim_post_check = select_first([LiftOverPlinkFiles.dataset_liftedover_bim, dataset_bim])
	File dataset_bed_post_check = select_first([LiftOverPlinkFiles.dataset_liftedover_bed, dataset_bed])
	File dataset_fam_post_check = select_first([LiftOverPlinkFiles.dataset_liftedover_fam, dataset_fam])


	Int removeHighMafAndChrXParVariants_size=2*bedSize + 10
	call RemoveHighMafAndChrXParVariants{
		input:
			dataset_bim=dataset_bim_post_check,
			dataset_bed=dataset_bed_post_check,
			dataset_fam=dataset_fam_post_check,
			disk_size  =removeHighMafAndChrXParVariants_size
	}

	Int changeCodingToChrPos=2*bedSize + 10
	call ChangeCodingToChrPos {
			input:
			dataset_bim=RemoveHighMafAndChrXParVariants.dataset_filtered_bim,
			dataset_bed=RemoveHighMafAndChrXParVariants.dataset_filtered_bed,
			dataset_fam=RemoveHighMafAndChrXParVariants.dataset_filtered_fam,
			disk_size = changeCodingToChrPos
	}

	call ExcludeDupsAndIndels{
		input:
			dataset_bim=ChangeCodingToChrPos.dataset_recoded_bim,
			dataset_bed=ChangeCodingToChrPos.dataset_recoded_bed,
			dataset_fam=ChangeCodingToChrPos.dataset_recoded_fam,
			disk_size = changeCodingToChrPos
	}

	call StrandAmbiguousHarmonization {
		input:
			dataset_bim=ExcludeDupsAndIndels.dataset_no_dups_bim,
			dataset_bed=ExcludeDupsAndIndels.dataset_no_dups_bed,
			dataset_fam=ExcludeDupsAndIndels.dataset_no_dups_fam,

			panel_freq=plink_frq,
			disk_size = changeCodingToChrPos
	}

	call PadhraigsScript {
		input:
			freq1=StrandAmbiguousHarmonization.dataset_strandAmbig,
			freq2=StrandAmbiguousHarmonization.panel_strandAmbig,
			disk_size = 100
	}

	call ExcludeFlipAndGetFreq {
		input:
			dataset_bim=ExcludeDupsAndIndels.dataset_no_dups_bim,
			dataset_bed=ExcludeDupsAndIndels.dataset_no_dups_bed,
			dataset_fam=ExcludeDupsAndIndels.dataset_no_dups_fam,

			ambiguous_flip=PadhraigsScript.ambiguous_flip,
			ambiguous_excl=PadhraigsScript.ambiguous_excl,
			disk_size = changeCodingToChrPos
	}

	call GetMoreExcludedAndPlot {
		input:
	 		dataset_freq=ExcludeFlipAndGetFreq.dataset_excluded_freq,
		 	panel_freq=plink_frq,
		 	disk_size = 10
	}

	call ExcludeFlipAndGetFreq as ExcludeFlipRestAndGetFreq {
		input:
			dataset_bim=ExcludeFlipAndGetFreq.dataset_excluded_bim,
			dataset_bed=ExcludeFlipAndGetFreq.dataset_excluded_bed,
			dataset_fam=ExcludeFlipAndGetFreq.dataset_excluded_fam,

			ambiguous_flip=GetMoreExcludedAndPlot.ambiguous_flip,
			ambiguous_excl=GetMoreExcludedAndPlot.ambiguous_excl,
			disk_size = changeCodingToChrPos
	}

	
	Int fasta_size=sub(size(reference_fasta,"GB") , "\\..*", "")  

	Int fixMaleHetSitesConvertToBCF_size=2 * bedSize + fasta_size + 100
	call FixMaleHetSitesConvertToBCF {
		input:
			dataset_bim=ExcludeFlipRestAndGetFreq.dataset_excluded_bim,
			dataset_bed=ExcludeFlipRestAndGetFreq.dataset_excluded_bed,
			dataset_fam=ExcludeFlipRestAndGetFreq.dataset_excluded_fam,
			dataset_hh=ExcludeFlipRestAndGetFreq.dataset_excluded_hh,

			reference_fasta=reference_fasta,
			reference_dict=reference_dict,
			reference_index=reference_index,

			disk_size = fixMaleHetSitesConvertToBCF_size
	}

	call SubTask as bcf_size_string{input: string=size(FixMaleHetSitesConvertToBCF.dataset_bcf,"GB"), pattern= "\\..*", replacement=""}
	Int bcf_size = bcf_size_string.out

	Int annotateVariants_size = 2 * bcf_size + dbSnpSize + 100
	call AnnotateVariants {
		input:
			dataset_bcf=FixMaleHetSitesConvertToBCF.dataset_bcf,
			dbSnp_vcf=dbSnp_vcf,
			dbSnp_vcf_index=dbSnp_vcf_index,
			disk_size = annotateVariants_size
	}		
	Int eagle_size = 2 * bcf_size + 20 

	scatter( chrom in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]) { 
		call PrePhaseVariantsEagle {
			input: 
				chrom=chrom,
			 	dataset_bcf=AnnotateVariants.dataset_annotated_bcf,
	 			genetic_map_file=genetic_maps_eagle,
	 			disk_size = eagle_size
		}
	}

	call PrePhaseVariantsEagle as PrePhaseVariantsEagleChrX  {
			input: 
				chrom="23",
			 	dataset_bcf=AnnotateVariants.dataset_annotated_bcf,
	 			genetic_map_file=genetic_maps_eagle,
	 			disk_size = eagle_size
	}

	call SubTask as chrx_vcf_size_string{input: string=size(PrePhaseVariantsEagleChrX.dataset_prephased_vcf, "GB") * 1.5, pattern= "\\..*", replacement=""}
	
	Int chrx_vcf_size = chrx_vcf_size_string.out + 100


	call ConcatenateStringToArray{ input: array=PrePhaseVariantsEagle.dataset_prephased_vcf,string=PrePhaseVariantsEagleChrX.dataset_prephased_vcf}

	call ConvertToHap {
		input:
			dataset_vcf=PrePhaseVariantsEagleChrX.dataset_prephased_vcf,
			disk_size = chrx_vcf_size
	}	
	
	call SampleLeftJoinFamGender {
		input:
			sample_file=ConvertToHap.output_hap,
			fam_file=ExcludeFlipRestAndGetFreq.dataset_excluded_fam,
			disk_size=100
	}

	output {
		File MAF_plot=GetMoreExcludedAndPlot.MAFs_plot
		File MAF_hist=GetMoreExcludedAndPlot.MAF_hist
		File sample_file=SampleLeftJoinFamGender.sample_file_out
		Boolean build_OK=CheckBuild.OK
		Array[File] dataset_bcf=ConcatenateStringToArray.array_out
	}	
}


############ 
### Tasks
############

task CheckBuild {
    File dataset_bim
	File dataset_bed
	File dataset_fam

	File dbSnp_RSIDs
	File dbSnp_vcf

	Int disk_size

	command <<<
		set -xeuo pipefail

		# extract dbsnp variants from dataset
		/usr/impute/plink/plink \
			--bim ${dataset_bim} \
			--bed ${dataset_bed} \
			--fam ${dataset_fam} \
			--extract ${dbSnp_RSIDs} \
			--make-bed \
			--out "dataset_b37"

		cut -f2 dataset_b37.bim > rsIds_subset 

		# extract the corresponding variants from dbsnp
		# vcfeval doesn't return error code different from 0
		# on error..thus the need to parse the logs and look
		# for an error string

		vcftools \
		--snps rsIds_subset \
		--out dbSnp_vcf.subset \
		--gzvcf ${dbSnp_vcf} \
		--recode \
		--recode-INFO-all 2>&1 |
			tee vcftools.log

		set +e
		grep -i "error" vcftools.log && echo "error in vcftools" && exit 2
		set -e

		# compare the locations of the variants from dbsnp and the dataset.
		# removing the set -e is required since diff will "fail" when the files differ.

		awk '{OFS="\t";print $1,$4,$2}' dataset_b37.bim | sed 's/^23\t/X\t/; s/^24\t/Y\t/'> subset_from_dataset.txt

		grep -v '^#' dbSnp_vcf.subset.recode.vcf | 
			cut -f1,2,3 > subset_from_dbSNP.txt


		set +e
		differentlines=$(diff --side-by-side --suppress-common-lines subset_from_dataset.txt subset_from_dbSNP.txt | wc -l)
		set -e 
		totallines=$(wc -l <subset_from_dbSNP.txt)
		
		# if differences are less than 1% --> OK
		if [[ $(echo "print(int($differentlines / ($totallines + 0.0) < 0.01)) " | python ) == 1 ]]; then 
			echo "true" > build_OK
		else
			echo "false" > build_OK
		fi		
	>>>
	output {
		Boolean OK = read_boolean("build_OK")
		File subset = "rsIds_subset"
		File subset_from_dataset = "subset_from_dataset.txt"
		File subset_from_dbSNP = "subset_from_dbSNP.txt"
		File dbSNP_subset_vcf = "dbSnp_vcf.subset.recode.vcf"
		File vcftools_log = "vcftools.log"
	}
	runtime {
		docker: "farjoun/impute:0.0.4-1506086533"
    	memory: "15 GB"
    	cpu: "3"
        disks: "local-disk " + disk_size + " HDD"
  }
}


#A) Get QC’d chip array data with denials and duplicated samples removed, sample IDs in correct format, and check that it is in build 37 (or respective to your panel build)
task LiftOverPlinkFiles {
	File dataset_bim
	File dataset_bed
	File dataset_fam

	File dbSnp_RSIDs
	File update_map

	Int disk_size

	command {
		set -e 
		/usr/impute/plink/plink --bim ${dataset_bim} --bed ${dataset_bed} --fam ${dataset_fam} \
			--extract ${dbSnp_RSIDs} \
			--update-map ${update_map} \
			--make-bed --out "dataset_lifted"
	}
	output {
		File dataset_liftedover_bim="dataset_lifted.bim"
		File dataset_liftedover_bed="dataset_lifted.bed"
		File dataset_liftedover_fam="dataset_lifted.fam"
	}
	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "32 GB"
    	cpu: "5"
        disks: "local-disk " + disk_size + " HDD"
  }
}

# B) Remove variants with PLINK MAF > 0.35 and remove PAR regions from chrX

task RemoveHighMafAndChrXParVariants{
	File dataset_bim
	File dataset_bed
	File dataset_fam

	Int disk_size

# 1) Move chrX PAR regions to chr code 25 with plink --split-x and remove variants with MAF>0.35 on the same go. If you did not do remapping in step A:
# 2) Take only chromosomes 1-23 (leave out chrX PAR regions with chr code 25 and chrY with chr code 24)
	command {
		set -euxo pipefail

		/usr/impute/plink/plink --bim ${dataset_bim} --bed ${dataset_bed} --fam ${dataset_fam} \
 			  --split-x b37 no-fail --max-maf 0.35 --make-bed --out dataset_mafover0.35_removed_X_split

		/usr/impute/plink/plink --bfile dataset_mafover0.35_removed_X_split --chr 1-23 --make-bed --out dataset_filtered

	}
	output {
		File dataset_filtered_bim="dataset_filtered.bim"
		File dataset_filtered_bed="dataset_filtered.bed"
		File dataset_filtered_fam="dataset_filtered.fam"
	}
	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "1 GB"
    	cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
  }
}


# C) Change variant coding to CHR_POS

task ChangeCodingToChrPos {
	File dataset_bim
	File dataset_bed
	File dataset_fam

	Int disk_size

	command <<<
		set -euxo pipefail

		# Create a bim file with new variant IDs (BIM files has the variant meta-information in the /usr/impute/plink/plink file format)
		awk -vOFS="\t" '{print $1, $1"_"$4,$3,$4,$5,$6}' ${dataset_bim} > dataset_recoded.bim
		cp ${dataset_bed} dataset_recoded.bed
		cp ${dataset_fam} dataset_recoded.fam
	>>>

	output {
		File dataset_recoded_bim="dataset_recoded.bim"
		File dataset_recoded_bed="dataset_recoded.bed"
		File dataset_recoded_fam="dataset_recoded.fam"
	}

	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "1 GB"
    	cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
  }
}


# D) Just to make sure: Exclude variants with duplicate IDs or I/D-coded indel variants

task ExcludeDupsAndIndels {
	File dataset_bim
	File dataset_bed
	File dataset_fam

	Int disk_size

	command <<<
		set -euxo pipefail
		
		touch dataset.dups
		cut -f 2 ${dataset_bim} | sort | uniq -d > dataset.dups

# 2) Remove these duplicates and DI-variants
		/usr/impute/plink/plink --bim ${dataset_bim} --bed ${dataset_bed} --fam ${dataset_fam} \
 			--exclude dataset.dups \
			--snps-only just-acgt \
			--make-bed \
			--out dataset_no_dups

# 3) Check (grep . will fail (i.e. return error code 1 if empty stream))

		## no I/D's
		[ $(cat dataset_no_dups.bim | awk '($5 =="I" || $5 =="D" || $6=="D" || $6=="I")' | wc -l) == 0 ] || ( echo "found I or D in bim file " && exit 2 ) 

		## no dups
		[ $(cat dataset_no_dups.bim | cut -f2 | sort | uniq -d | wc -l ) == 0 ] || ( echo "found duplicates entries in bim file" && exit 3 )
	>>>

	output {
		File dataset_no_dups_bim="dataset_no_dups.bim"
		File dataset_no_dups_bed="dataset_no_dups.bed"
		File dataset_no_dups_fam="dataset_no_dups.fam"
	}

	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "1 GB"
    	cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
  }
}


# E) Strand ambiguous harmonization (Padhraig Gormley's script)

task StrandAmbiguousHarmonization {
	File dataset_bim
	File dataset_bed
	File dataset_fam

	File panel_freq

	Int disk_size

	command <<<
		set -e 

		/usr/impute/plink/plink --bim ${dataset_bim} --bed ${dataset_bed} --fam ${dataset_fam} \
 		 --freq --out dataset_freq

		# Pick strand ambiguous SNP frequencies for chip and panel (please see the pre-requirements part in the beginning of this document)
		awk '$3=="A" && $4=="T" || $3=="T" && $4=="A" || $3=="G" && $4=="C" || $3=="C" && $4=="G"' dataset_freq.frq > dataset_freq.strandAmbig.frq
		awk  '$3=="A" && $4=="T" || $3=="T" && $4=="A" || $3=="G" && $4=="C" || $3=="C" && $4=="G"' ${panel_freq} > panel.strandAmbig.frq
		
	>>>

	output {
		File dataset_strandAmbig="dataset_freq.strandAmbig.frq"
		File panel_strandAmbig="panel.strandAmbig.frq" 
	}

	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "1 GB"
    	cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
  }
}


task PadhraigsScript {
	File freq1
	File freq2

	Int disk_size

	command <<<

		set -e
# Strand ambiguous variant flipping and exclusion
# from Padhraig Gormley (pgormley@broadinstitute.org)

# $1 - plink .frq from fileset1 with strand-ambiguous SNPs only
# $2 - plink .frq from fileset2 with strand-ambiguous SNPs only
		awk -v OFS="\t" 'FNR==NR{ maf[$2]=$5; minor[$2]=$3; major[$2]=$4; next} \
    	($2 in maf){ \
            	if( $5 > 0.35 ){ \
                    	print $2,minor[$2],major[$2],maf[$2],$3,$4,$5,"EXCLUDE" \
            	} \
            	else if( minor[$2]==$3 && sqrt((maf[$2]-$5)^2)>0.05 ){ \
                    	print $2,minor[$2],major[$2],maf[$2],$3,$4,$5,"EXCLUDE" \
            	} \
            	else if( minor[$2]==$3 && sqrt((maf[$2]-$5)^2)<=0.05 ){ \
                    	print $2,minor[$2],major[$2],maf[$2],$3,$4,$5,"OK" \
            	} \
            	else if( minor[$2]==$4 && sqrt((maf[$2]-$5)^2)>0.05 ){ \
                    	print $2,minor[$2],major[$2],maf[$2],$3,$4,$5,"EXCLUDE" \
            	} \
            	else if( minor[$2]==$4 && sqrt((maf[$2]-$5)^2)<=0.05 ){ \
                    	print $2,minor[$2],major[$2],maf[$2],$3,$4,$5,"FLIP" \
            	} \
            	else{ \
                    	print $2,minor[$2],major[$2],maf[$2],$3,$4,$5,"EXCLUDE" \
            	} \
    	}' ${freq1} ${freq2} > out.frq

    	# Gather list of flippable and excludable variants
		cat out.frq | grep FLIP | awk '{print $1}' > strandAmb.res.FLIP
		cat out.frq | grep EXCL | awk '{print $1}' > strandAmb.res.EXCL

	>>>
	output {
		File ambiguous_flip="strandAmb.res.FLIP"
		File ambiguous_excl="strandAmb.res.EXCL"
	}
	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "1 GB"
    	cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
  }
}

 # F) Flip flippable, exclude non-flippable SNVs and calculate freqs



task ExcludeFlipAndGetFreq {
	File dataset_bim
	File dataset_bed
	File dataset_fam
	
	File ambiguous_flip
	File ambiguous_excl

	Int disk_size

	command <<<

		set -euxo 

		touch dataset_excluded_and_flipped.hh
	 	/usr/impute/plink/plink --bim ${dataset_bim} --bed ${dataset_bed} --fam ${dataset_fam} \
			--exclude ${ambiguous_excl} \
			--flip  ${ambiguous_flip} \
			--make-bed \
			--freq \
			--out dataset_excluded_and_flipped | tee flip.log

		# check that
		# - Flipped plink log file reports same number as in ../E_harmonization_script/strandAmb.res.FLIP
		# - Excluded number (variants in - variants out) same as in ../E_harmonization_script/strandAmb.res.EXCL

		in=$(wc -l <${dataset_bim})
		out=$(wc -l <dataset_excluded_and_flipped.bim)
		excl=$(wc -l <${ambiguous_excl})
		freq=$(wc -l <dataset_excluded_and_flipped.frq)

		# verify that freq count is the same as bim count (subtract 1 for the header line in the frq file):
		[[ $(( $freq - 1 )) == $out ]]  || 
			(echo "ERROR: Number of lines in frq file isn't equal to the number of lines in bim file" && exit 1)

		# will return non zero error code if false
		[[ $(( $in - $out )) == $excl ]] || 
			(echo "ERROR: Number of excluded snps not matching difference between before vs. after." && exit 2)

		# Check that number of variants flipped matches with NORM_multiallelic_flip.txt. 
		
		flipped=$(sed -n 's/^--flip: \(.*\) SNP.* flipped\..*$/\1/p' flip.log)

		# will return non zero error code if false
		[[ $(wc -l <${ambiguous_flip}) == $flipped ]] || 
			(echo "ERROR: Number of flipped SNPs reported differs from number of lines in input file." && exit 3)

	>>>
	output {
		File dataset_excluded_bim ="dataset_excluded_and_flipped.bim"
		File dataset_excluded_bed ="dataset_excluded_and_flipped.bed"
		File dataset_excluded_fam ="dataset_excluded_and_flipped.fam"
		File dataset_excluded_hh  ="dataset_excluded_and_flipped.hh"
		File dataset_excluded_freq="dataset_excluded_and_flipped.frq"
	}
	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "1 GB"
    	cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
  }
}


# G) Create a list of list of additional flips and exclusions + few diagnostic plots (MAF histogram and scatter plot)

task GetMoreExcludedAndPlot {
	File dataset_freq
	File panel_freq

	Int disk_size

	String MAF="MAFs.jpg"
	String Hist="MAF_hist.jpg"
	String Flip="flip.txt"
	String Excl="exclude.txt"


	command <<<

		# Then run
		Rscript --no-save -<<'RCODE'

		library(data.table) # For fast fread

		# Create a helper function to flip alleles
		flip<-function(x) ifelse(x=="A","T", ifelse(x=="T","A", ifelse(x=="G","C", ifelse(x=="C","G","Error"))))

		cat("Defined flip function\n")

		# Read and merge frequency files (take intersection)
		# DEFINE FREQUENCY FILES HERE (please check pre-requirements from the beginning of this document)
		chip<-fread("${dataset_freq}", header=T)

		cat("read dataset frequency file\n")

		panel<-fread("${panel_freq}", header=T)

		cat("read panel frequency file\n")

		isec<-merge(panel, chip, by="SNP")

		cat("merged files\n")

		# Future proof: check if indels present
		indel<-nchar(isec$A1.x) != 1 | nchar(isec$A2.x) != 1 | nchar(isec$A1.y) != 1 | nchar(isec$A2.y) != 1

		# A1 and A2 alleles match exactly
		concordant<-(isec$A1.x == isec$A1.y & isec$A2.x == isec$A2.y)

		# A1 and A2 alleles match exactly after flipping (exclude indels from this check)
		flippable<-!indel & (isec$A1.x == flip(isec$A1.y) & isec$A2.x == flip(isec$A2.y))

		# Check that MAFs is within range of 5 pp in both datasets
		maf_ok<-abs(isec$MAF.x - isec$MAF.y)<0.05

		# Exclude those that are neither concordant nor flippable
		exclude<-!(concordant | flippable) | !maf_ok

		cat("preparing to plot MAFs\n")

		jpeg("${MAF}")
		# Plot first concordant variants, then flippable and finally excludable
		plot(isec[concordant]$MAF.x, isec[concordant]$MAF.y, col=1,pch=20, main="Chip data MAF x WES+WGS panel MAF (intersecting positions)", xlab="WGS+WES MAF", ylab="Chip MAF")
		points(isec[flippable]$MAF.x, isec[flippable]$MAF.y, col=4,pch=20)
		points(isec[exclude]$MAF.x, isec[exclude]$MAF.y, col=2,pch=20)

		# Draw a legend
		legend("topleft", legend=c("Concordant","Flippable","Mismatching alleles or MAF"), col=c(1,4,2), pch=20, cex=0.9)
		dev.off()

		cat("preparing to plot histogram.\n")

		jpeg("${Hist}")
		# Chip MAF histogram for intersecting variants
		hist(isec[!exclude]$MAF.y, breaks=100, main="Chip MAF for isec variants")
		dev.off()

		cat("plotted.\n")

		# Write flipping and exclusion lists
		write.table(isec[flippable]$SNP, "${Flip}", quote=F, row.names=F, col.names=F)

		cat("wrote flip.txt\n")

		# IF YOU DIDN’T DO A BUILD CONVERSION IN THE STEP A, use following code:
	    write.table(isec[exclude]$SNP, "exclude.txt", quote=F, row.names=F, col.names=F) 

		# IF YOU DID A BUILD CONVERSION IN THE STEP A: 
		# Please comment out red write.table above and uncomment following three yellow commands. We want to remove all remapped variants that are not in the panel and cannot be checked for frequency (most 
		# will get imputed back) 

		# non_panel<-!(chip$SNP %in% isec$SNP)
		# exclude_w_nonpanel<-union(isec[exclude]$SNP, chip[non_panel]$SNP)
		# write.table(exclude_w_nonpanel, "${Excl}", quote=F, row.names=F, col.names=F)

		cat("wrote exclude.txt\n")

		RCODE
	>>>
	output {
		File MAFs_plot="${MAF}"
		File MAF_hist="${Hist}"

		File ambiguous_flip="${Flip}"
		File ambiguous_excl="${Excl}"
	}
	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "10 GB"
    	cpu: "2"
        disks: "local-disk " + disk_size + " HDD"
  }
}

# H) Flip rest flippable, exclude non-flippable
# same task as F


# I) Extract gender before conversion, fix het sites, convert chip to VCF/BCF and check against human reference


task FixMaleHetSitesConvertToBCF {
	File dataset_bim
	File dataset_bed
	File dataset_fam
	File dataset_hh

	File reference_fasta
	File reference_dict
	File reference_index

	Int disk_size

	command <<<
		set -euxo pipefail

		# 1) Extract gender information
		cat ${dataset_fam} | awk '($5 == 1) {print $1,$2}' > men_samples.txt

		# 2) Check sites with excess heterozygosity:
		# First cat command may give warning / error if .hh file doesn’t exist - do not worry, this is a good sign (= there are no heterozygosity errors) -> ignore and continue normally
		cat ${dataset_hh} | cut -f3 | sort | uniq -c > het_counts.txt
		cat het_counts.txt | awk '($1>5) {print $2}' > het_over5males.EXCL

		# 3) Exclude identified sites (with over 5 male het calls) and set rest heterozygous male chrX calls missing:
		/usr/impute/plink/plink --bim ${dataset_bim} --bed ${dataset_bed} --fam ${dataset_fam} \
		--exclude het_over5males.EXCL --set-hh-missing --make-bed --out dataset_het_fixed

		# 4) VCF/BCF CONVERSION and reference assembly check
		# 4a) RUN I: Convert plink file into a VCF (bgzipped)
		# flag “--output-chr MT” codes chrX as “X”, not “23”. This is necessary for checking alleles against human reference assembly (VCF uses “X” not “23”)
		/usr/impute/plink/plink --bfile dataset_het_fixed --recode vcf-iid bgz --output-chr MT --out dataset_vcf

		# 4b) Check data against human reference assembly (used in panels).
		# 	On the same run generate a list of mismatching variants (will end up multiallelic due to an allelic mismatch with the assembly)
		bcftools norm -f ${reference_fasta} -c ws dataset_vcf.vcf.gz -Ou --threads 6 | 
			bcftools view -m 3 --threads 6 -Ou | bcftools query -f '%ID\n' > NORM_multiallelic_flip.txt


		# 4c) RUN II: 3) Convert plink file into a VCF (bgzipped) with NORM_multiallelic_flip.txt sites flipped
		/usr/impute/plink/plink --bfile dataset_het_fixed --recode vcf-iid bgz --output-chr MT \
		--flip NORM_multiallelic_flip.txt \
		--out dataset_normalized_vcf

		# Check that number of variants flipped matches with NORM_multiallelic_flip.txt. 
		# TODO: HOW?

		# 4d) Normalize again and remove possible remaining multiallelic sites (= sites not fixed by flipping)
		bcftools norm -f ${reference_fasta} -c ws dataset_normalized_vcf.vcf.gz  -Ou --threads 6 | bcftools view -m2 -M2 -Ob -o dataset_bcf.bcf.gz --threads 6


	>>>
	output {
		File dataset_bcf="dataset_bcf.bcf.gz"
	}
	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "7 GB"
    	cpu: "10"
        disks: "local-disk " + disk_size + " HDD"
  }
}


# J -- Annotate variant IDs with VAR_CHROM_POS_REF_ALT and dbSNP (RSID ANNOTATION MAY NOT BE NECESSARY, IF IMPUTATION IS THE ONLY PURPOSE)

task AnnotateVariants {
	File dataset_bcf
	File dbSnp_vcf
	File dbSnp_vcf_index
	

	Int disk_size
	command <<<
		set -euxo pipefail

		# 1) Replace all variant IDs with CHROM_POS_REF_ALT 
		bcftools annotate ${dataset_bcf} \
		--set-id 'VAR\_%CHROM\_%POS\_%REF\_%FIRST_ALT' -Ob --threads 6 \
		-o dataset_VARid.bcf.gz

		# 2) Tabix resulting file (needed for dbSNP annotation)
		tabix dataset_VARid.bcf.gz

		# 3) Annotate variants IDs from dbSNP for those sites that have a matching record
		bcftools annotate -a ${dbSnp_vcf} -c ID \
		 dataset_VARid.bcf.gz \
		-Ob --threads 6 \
		-o dataset_dbSNP.bcf.gz

	>>>
	output {
		File dataset_annotated_bcf="dataset_dbSNP.bcf.gz"
	}
	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "7 GB"
    	cpu: "7"
  	    disks: "local-disk " + disk_size + " HDD"
  }
}


# K) Pre-Phasing
# ----------------

task PrePhaseVariantsEagle {
	File dataset_bcf
	String chrom
	String genetic_map_file

	Int disk_size

	command <<<
		  /usr/impute/Eagle_v2.3.5/eagle  \
           	--vcf ${dataset_bcf} \
           	--chrom ${chrom} \
           	--geneticMapFile ${genetic_map_file} \
           	--numThreads=8 \
           	--Kpbwt=20000 \
           	--outPrefix pre_phased_${chrom}
	>>>
	output {
		File dataset_prephased_vcf="pre_phased_${chrom}.vcf.gz"
	}
	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "32 GB"
    	cpu: "8"
        disks: "local-disk " + disk_size + " HDD"
  }
}

# K) Convert to Hap
# ----------------

task ConvertToHap {
	File dataset_vcf

	Int disk_size
	
	command {
		bcftools convert ${dataset_vcf} --hapsample file --vcf-ids
	}
	
	output {
		File output_hap="file.sample"
	}
	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "7 GB"	
    	cpu: "8"
        disks: "local-disk " + disk_size + " HDD"
  }
}


task SampleLeftJoinFamGender {
	File sample_file
	File fam_file
	
	String dollar_sign="$"

	Int disk_size
	command <<<

		set -euxo pipefail

		cat <<'SCALADOC' > script.scala
			
			import scala.io._

			def require0(cond: Boolean, msg: String) = {
			  if (cond == false) {
				println(s"error ::  $msg")
				sys.exit(1)
			  }
			}

			case class Sample(id1: String, id2: String, missing: String) {
			  override def toString() = this.productIterator.toArray.mkString(" ")
			}

			object Sample {
			  def apply(xs: Array[String]): Sample = Sample(xs(0), xs(1), xs(2))
			}

			case class Fam(id1: String, id2: String, sex: String)

			object Fam {
			  def apply(xs: Array[String]): Fam = Fam(xs(0), xs(1), xs(4))
			}

			val samples = Source.fromFile(args(0)).getLines.drop(2).flatMap { line =>
			            	val xs = line.trim.split("\\s+")
			            	if (xs.size > 2) Some(Sample(xs)) else None
			              }.toArray

			val fams = Source.fromFile(args(1)).getLines.flatMap { line =>
			             val xs = line.trim.split("\\s+")
			             if (xs.size > 4) Some(Fam(xs)) else None
			     	   }.toArray

			val samplesMap = samples.map ( s => ((s.id1, s.id2)) -> s ).toMap
			val famsMap    = fams.map ( f => ((f.id1, f.id2)) -> f ).toMap

			val samplesIds = samplesMap.keySet
			val famsIds    = famsMap.keySet

			require0( samples.size == samplesIds.size, "duplicate ids in samples" )
			require0( fams.size    == famsIds.size,	"duplicate ids in fam" )

			if ( samplesIds != famsIds ) {
			  val msg1 = "sample ids not in fam ids " + samplesIds.diff(famsIds)
			  val msg2 = "fam ids not in sample ids " + famsIds.diff(samplesIds)
			  require0( false, s"sample ids set differs from fam ids set\n  $msg1\n  $msg2" )
			}

			println("ID_1 ID_2 missing sex")
			println("0 0 0 D")
			samples.foreach { s => println( s"${dollar_sign}s ${dollar_sign}{famsMap((s.id1, s.id2)).sex}" ) }

		SCALADOC

		/root/scala-2.11.6/bin/scala script.scala ${sample_file} ${fam_file} > new_sample_file.sample

	>>>
	output {
		File sample_file_out="new_sample_file.sample"
	}
	runtime {
		docker: "farjoun/impute:0.0.3-1504715575"
    	memory: "1 GB"
    	cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
  }
}


#→ Now we should be ready to impute!


############### helper tasks #######################

task ExtractRsIdFromVCF {
	File vcf

	Int disk_size
	command <<<
		set -exuo pipefail

		gunzip -cd ${vcf} | grep -v "^#" | cut -f 2,3 | tee >(awk '$2!="."{print $2" "$1}' > update_map.txt ) | cut -f 2  | grep -v '\.' > rsIDs.txt
	>>>
	runtime {
      docker: "python:2.7"
      memory: "1 GB"
      disks: "local-disk " + disk_size + " HDD"
    }
    output {
        File rsIDs="rsIDs.txt"
        File update_map="update_map.txt"
    }
}

task ConcatenateStringToArray{
	Array[String] array
	String string

	command <<<
		set -exuo pipefail

		echo -n ${sep='\t' array}
		echo -e '\t'${string}
	>>>
	runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }
    output {
        Array[String] array_out = read_lines(stdout())
    }
}


 task SubTask{

 	String string
 	String pattern
 	String replacement

 	command{}

 	output {
 		String out=sub(string, pattern, replacement)
 	}

 	runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }
 }


