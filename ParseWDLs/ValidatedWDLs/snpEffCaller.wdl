workflow SnpEff {
		File snpEff_jar
        File inputVCF
        File configFile
        File databaseFile
        String outputName
        
        call snpEffCaller {
         input:
           snpEff_jar = snpEff_jar,
           inputVCF = inputVCF,
           configFile = configFile,
           databaseFile = databaseFile,
           outputName = outputName,
        }
}

task snpEffCaller {
        File snpEff_jar
        File inputVCF
        File databaseFile
        String outputName
        File configFile

        command {
        		subdir=$(dirname ${snpEff_jar})
        		unzip ${databaseFile} -d $subdir
                java -Xmx4g -jar ${snpEff_jar} -v -config ${configFile} -canon -nodownload -hgvs1LetterAa hg19 ${inputVCF} > ${outputName}_snpEff.vcf
        }

        output {
                File vcfOut = "${outputName}_snpEff.vcf"
        }
        
		runtime {
        	docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
            memory: "4GB"
            disks: "local-disk 200 HDD"
    }
}