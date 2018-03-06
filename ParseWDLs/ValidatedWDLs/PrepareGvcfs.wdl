
task GetUUID {
    command <<<
        python <<CODE
import uuid; 
print(str(uuid.uuid4()))
CODE
     >>>

      runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }

    output {
        String uuid = read_string(stdout())
    }
}

 task UniquifyGvcf {
      File file
      String uuid
      Int disk_size
      command <<<
        out=${uuid}.g.vcf.gz &&
        cp ${file} $out
     >>>

      runtime {
      docker: "python:2.7"
      memory: "1 GB"
      disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_gvcf = "${uuid}.g.vcf.gz" 
    }
}

task Uniquify{
      File file
      Int disk_size
      command <<<
        out=$(python <(echo "print(str(uuid.uuid4()))").g.vcf.gz &&
        cp ${file} $out
        echo $out > "file_name.txt"
     >>>

      runtime {
      docker: "python:2.7"
      memory: "1 GB"
      disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_gvcf = read_string("file_name.txt") 
    }
}



task FileToStringArray {
    File input_file

    command {

    }

    runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }

    output {
        Array[String] array = read_lines(input_file)
    }
}

task MakeSampleMap {
  Array[String] samples
  Array[String] gvcfs

  command {

    python<<CODE

samples=[s.strip() for s in open("${write_lines(samples)}","r").readlines()]
gvcfs = [s.strip() for s in open("${write_lines(gvcfs)}","r").readlines()]

assert len(samples)==len(gvcfs) , "lengths of gvcfs ands samples arrays are not equal"
assert len(set(samples))==len(samples), "samples list is not unique, please fix"
assert len(set(gvcfs))==len(gvcfs), "gvcfs list is not unique, please fix"

w = open("sample_map.txt","wt")
for (sample,gvcf) in zip(samples,gvcfs):
    w.write( "%s\t%s\n"%(gvcf.split("/")[-1], sample))
w.close()
CODE

  }
  output {
    File sample_map = "sample_map.txt"
  }

  runtime {
      docker: "python:2.7"
      memory: "1 GB"
  }
}

#for 901 intervals only
task SplitGvcf {
  File gvcf
  File interval_list
  String sample_name
  Int disk_size

  command <<<
    # cut -f1-3 returns <chromosome> <start> <stop>
    set -o pipefail


    /usr/gitc/tabix ${gvcf} &&
    cat ${interval_list} | grep -v "@" | cut -f1-3 > regions.txt &&
    mkdir split_gvcfs &&
    piece=0 &&
    while read -r chrom start stop; do
      # might need to change %04 if we decide to have num_partitions greater than 9999
      OUT="printf ${sample_name}.%04d.g.vcf.gz $piece"
      /usr/gitc/tabix -h ${gvcf} $chrom:$start-$stop | /usr/gitc/bgzip > split_gvcfs/$($OUT)
      piece=$(($piece+1))
    done < regions.txt

  >>>
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.1-1465310914"
    memory: "3 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Array[File] gvcf_list = ["split_gvcfs/${sample_name}.0000.g.vcf.gz", "split_gvcfs/${sample_name}.0001.g.vcf.gz", "split_gvcfs/${sample_name}.0002.g.vcf.gz", "split_gvcfs/${sample_name}.0003.g.vcf.gz", "split_gvcfs/${sample_name}.0004.g.vcf.gz", "split_gvcfs/${sample_name}.0005.g.vcf.gz", "split_gvcfs/${sample_name}.0006.g.vcf.gz", "split_gvcfs/${sample_name}.0007.g.vcf.gz", "split_gvcfs/${sample_name}.0008.g.vcf.gz", "split_gvcfs/${sample_name}.0009.g.vcf.gz", "split_gvcfs/${sample_name}.0010.g.vcf.gz", "split_gvcfs/${sample_name}.0011.g.vcf.gz", "split_gvcfs/${sample_name}.0012.g.vcf.gz", "split_gvcfs/${sample_name}.0013.g.vcf.gz", "split_gvcfs/${sample_name}.0014.g.vcf.gz", "split_gvcfs/${sample_name}.0015.g.vcf.gz", "split_gvcfs/${sample_name}.0016.g.vcf.gz", "split_gvcfs/${sample_name}.0017.g.vcf.gz", "split_gvcfs/${sample_name}.0018.g.vcf.gz", "split_gvcfs/${sample_name}.0019.g.vcf.gz", "split_gvcfs/${sample_name}.0020.g.vcf.gz", "split_gvcfs/${sample_name}.0021.g.vcf.gz", "split_gvcfs/${sample_name}.0022.g.vcf.gz", "split_gvcfs/${sample_name}.0023.g.vcf.gz", "split_gvcfs/${sample_name}.0024.g.vcf.gz", "split_gvcfs/${sample_name}.0025.g.vcf.gz", "split_gvcfs/${sample_name}.0026.g.vcf.gz", "split_gvcfs/${sample_name}.0027.g.vcf.gz", "split_gvcfs/${sample_name}.0028.g.vcf.gz", "split_gvcfs/${sample_name}.0029.g.vcf.gz", "split_gvcfs/${sample_name}.0030.g.vcf.gz", "split_gvcfs/${sample_name}.0031.g.vcf.gz", "split_gvcfs/${sample_name}.0032.g.vcf.gz", "split_gvcfs/${sample_name}.0033.g.vcf.gz", "split_gvcfs/${sample_name}.0034.g.vcf.gz", "split_gvcfs/${sample_name}.0035.g.vcf.gz", "split_gvcfs/${sample_name}.0036.g.vcf.gz", "split_gvcfs/${sample_name}.0037.g.vcf.gz", "split_gvcfs/${sample_name}.0038.g.vcf.gz", "split_gvcfs/${sample_name}.0039.g.vcf.gz", "split_gvcfs/${sample_name}.0040.g.vcf.gz", "split_gvcfs/${sample_name}.0041.g.vcf.gz", "split_gvcfs/${sample_name}.0042.g.vcf.gz", "split_gvcfs/${sample_name}.0043.g.vcf.gz", "split_gvcfs/${sample_name}.0044.g.vcf.gz", "split_gvcfs/${sample_name}.0045.g.vcf.gz", "split_gvcfs/${sample_name}.0046.g.vcf.gz", "split_gvcfs/${sample_name}.0047.g.vcf.gz", "split_gvcfs/${sample_name}.0048.g.vcf.gz", "split_gvcfs/${sample_name}.0049.g.vcf.gz", "split_gvcfs/${sample_name}.0050.g.vcf.gz", "split_gvcfs/${sample_name}.0051.g.vcf.gz", "split_gvcfs/${sample_name}.0052.g.vcf.gz", "split_gvcfs/${sample_name}.0053.g.vcf.gz", "split_gvcfs/${sample_name}.0054.g.vcf.gz", "split_gvcfs/${sample_name}.0055.g.vcf.gz", "split_gvcfs/${sample_name}.0056.g.vcf.gz", "split_gvcfs/${sample_name}.0057.g.vcf.gz", "split_gvcfs/${sample_name}.0058.g.vcf.gz", "split_gvcfs/${sample_name}.0059.g.vcf.gz", "split_gvcfs/${sample_name}.0060.g.vcf.gz", "split_gvcfs/${sample_name}.0061.g.vcf.gz", "split_gvcfs/${sample_name}.0062.g.vcf.gz", "split_gvcfs/${sample_name}.0063.g.vcf.gz", "split_gvcfs/${sample_name}.0064.g.vcf.gz", "split_gvcfs/${sample_name}.0065.g.vcf.gz", "split_gvcfs/${sample_name}.0066.g.vcf.gz", "split_gvcfs/${sample_name}.0067.g.vcf.gz", "split_gvcfs/${sample_name}.0068.g.vcf.gz", "split_gvcfs/${sample_name}.0069.g.vcf.gz", "split_gvcfs/${sample_name}.0070.g.vcf.gz", "split_gvcfs/${sample_name}.0071.g.vcf.gz", "split_gvcfs/${sample_name}.0072.g.vcf.gz", "split_gvcfs/${sample_name}.0073.g.vcf.gz", "split_gvcfs/${sample_name}.0074.g.vcf.gz", "split_gvcfs/${sample_name}.0075.g.vcf.gz", "split_gvcfs/${sample_name}.0076.g.vcf.gz", "split_gvcfs/${sample_name}.0077.g.vcf.gz", "split_gvcfs/${sample_name}.0078.g.vcf.gz", "split_gvcfs/${sample_name}.0079.g.vcf.gz", "split_gvcfs/${sample_name}.0080.g.vcf.gz", "split_gvcfs/${sample_name}.0081.g.vcf.gz", "split_gvcfs/${sample_name}.0082.g.vcf.gz", "split_gvcfs/${sample_name}.0083.g.vcf.gz", "split_gvcfs/${sample_name}.0084.g.vcf.gz", "split_gvcfs/${sample_name}.0085.g.vcf.gz", "split_gvcfs/${sample_name}.0086.g.vcf.gz", "split_gvcfs/${sample_name}.0087.g.vcf.gz", "split_gvcfs/${sample_name}.0088.g.vcf.gz", "split_gvcfs/${sample_name}.0089.g.vcf.gz", "split_gvcfs/${sample_name}.0090.g.vcf.gz", "split_gvcfs/${sample_name}.0091.g.vcf.gz", "split_gvcfs/${sample_name}.0092.g.vcf.gz", "split_gvcfs/${sample_name}.0093.g.vcf.gz", "split_gvcfs/${sample_name}.0094.g.vcf.gz", "split_gvcfs/${sample_name}.0095.g.vcf.gz", "split_gvcfs/${sample_name}.0096.g.vcf.gz", "split_gvcfs/${sample_name}.0097.g.vcf.gz", "split_gvcfs/${sample_name}.0098.g.vcf.gz", "split_gvcfs/${sample_name}.0099.g.vcf.gz", "split_gvcfs/${sample_name}.0100.g.vcf.gz", "split_gvcfs/${sample_name}.0101.g.vcf.gz", "split_gvcfs/${sample_name}.0102.g.vcf.gz", "split_gvcfs/${sample_name}.0103.g.vcf.gz", "split_gvcfs/${sample_name}.0104.g.vcf.gz", "split_gvcfs/${sample_name}.0105.g.vcf.gz", "split_gvcfs/${sample_name}.0106.g.vcf.gz", "split_gvcfs/${sample_name}.0107.g.vcf.gz", "split_gvcfs/${sample_name}.0108.g.vcf.gz", "split_gvcfs/${sample_name}.0109.g.vcf.gz", "split_gvcfs/${sample_name}.0110.g.vcf.gz", "split_gvcfs/${sample_name}.0111.g.vcf.gz", "split_gvcfs/${sample_name}.0112.g.vcf.gz", "split_gvcfs/${sample_name}.0113.g.vcf.gz", "split_gvcfs/${sample_name}.0114.g.vcf.gz", "split_gvcfs/${sample_name}.0115.g.vcf.gz", "split_gvcfs/${sample_name}.0116.g.vcf.gz", "split_gvcfs/${sample_name}.0117.g.vcf.gz", "split_gvcfs/${sample_name}.0118.g.vcf.gz", "split_gvcfs/${sample_name}.0119.g.vcf.gz", "split_gvcfs/${sample_name}.0120.g.vcf.gz", "split_gvcfs/${sample_name}.0121.g.vcf.gz", "split_gvcfs/${sample_name}.0122.g.vcf.gz", "split_gvcfs/${sample_name}.0123.g.vcf.gz", "split_gvcfs/${sample_name}.0124.g.vcf.gz", "split_gvcfs/${sample_name}.0125.g.vcf.gz", "split_gvcfs/${sample_name}.0126.g.vcf.gz", "split_gvcfs/${sample_name}.0127.g.vcf.gz", "split_gvcfs/${sample_name}.0128.g.vcf.gz", "split_gvcfs/${sample_name}.0129.g.vcf.gz", "split_gvcfs/${sample_name}.0130.g.vcf.gz", "split_gvcfs/${sample_name}.0131.g.vcf.gz", "split_gvcfs/${sample_name}.0132.g.vcf.gz", "split_gvcfs/${sample_name}.0133.g.vcf.gz", "split_gvcfs/${sample_name}.0134.g.vcf.gz", "split_gvcfs/${sample_name}.0135.g.vcf.gz", "split_gvcfs/${sample_name}.0136.g.vcf.gz", "split_gvcfs/${sample_name}.0137.g.vcf.gz", "split_gvcfs/${sample_name}.0138.g.vcf.gz", "split_gvcfs/${sample_name}.0139.g.vcf.gz", "split_gvcfs/${sample_name}.0140.g.vcf.gz", "split_gvcfs/${sample_name}.0141.g.vcf.gz", "split_gvcfs/${sample_name}.0142.g.vcf.gz", "split_gvcfs/${sample_name}.0143.g.vcf.gz", "split_gvcfs/${sample_name}.0144.g.vcf.gz", "split_gvcfs/${sample_name}.0145.g.vcf.gz", "split_gvcfs/${sample_name}.0146.g.vcf.gz", "split_gvcfs/${sample_name}.0147.g.vcf.gz", "split_gvcfs/${sample_name}.0148.g.vcf.gz", "split_gvcfs/${sample_name}.0149.g.vcf.gz", "split_gvcfs/${sample_name}.0150.g.vcf.gz", "split_gvcfs/${sample_name}.0151.g.vcf.gz", "split_gvcfs/${sample_name}.0152.g.vcf.gz", "split_gvcfs/${sample_name}.0153.g.vcf.gz", "split_gvcfs/${sample_name}.0154.g.vcf.gz", "split_gvcfs/${sample_name}.0155.g.vcf.gz", "split_gvcfs/${sample_name}.0156.g.vcf.gz", "split_gvcfs/${sample_name}.0157.g.vcf.gz", "split_gvcfs/${sample_name}.0158.g.vcf.gz", "split_gvcfs/${sample_name}.0159.g.vcf.gz", "split_gvcfs/${sample_name}.0160.g.vcf.gz", "split_gvcfs/${sample_name}.0161.g.vcf.gz", "split_gvcfs/${sample_name}.0162.g.vcf.gz", "split_gvcfs/${sample_name}.0163.g.vcf.gz", "split_gvcfs/${sample_name}.0164.g.vcf.gz", "split_gvcfs/${sample_name}.0165.g.vcf.gz", "split_gvcfs/${sample_name}.0166.g.vcf.gz", "split_gvcfs/${sample_name}.0167.g.vcf.gz", "split_gvcfs/${sample_name}.0168.g.vcf.gz", "split_gvcfs/${sample_name}.0169.g.vcf.gz", "split_gvcfs/${sample_name}.0170.g.vcf.gz", "split_gvcfs/${sample_name}.0171.g.vcf.gz", "split_gvcfs/${sample_name}.0172.g.vcf.gz", "split_gvcfs/${sample_name}.0173.g.vcf.gz", "split_gvcfs/${sample_name}.0174.g.vcf.gz", "split_gvcfs/${sample_name}.0175.g.vcf.gz", "split_gvcfs/${sample_name}.0176.g.vcf.gz", "split_gvcfs/${sample_name}.0177.g.vcf.gz", "split_gvcfs/${sample_name}.0178.g.vcf.gz", "split_gvcfs/${sample_name}.0179.g.vcf.gz", "split_gvcfs/${sample_name}.0180.g.vcf.gz", "split_gvcfs/${sample_name}.0181.g.vcf.gz", "split_gvcfs/${sample_name}.0182.g.vcf.gz", "split_gvcfs/${sample_name}.0183.g.vcf.gz", "split_gvcfs/${sample_name}.0184.g.vcf.gz", "split_gvcfs/${sample_name}.0185.g.vcf.gz", "split_gvcfs/${sample_name}.0186.g.vcf.gz", "split_gvcfs/${sample_name}.0187.g.vcf.gz", "split_gvcfs/${sample_name}.0188.g.vcf.gz", "split_gvcfs/${sample_name}.0189.g.vcf.gz", "split_gvcfs/${sample_name}.0190.g.vcf.gz", "split_gvcfs/${sample_name}.0191.g.vcf.gz", "split_gvcfs/${sample_name}.0192.g.vcf.gz", "split_gvcfs/${sample_name}.0193.g.vcf.gz", "split_gvcfs/${sample_name}.0194.g.vcf.gz", "split_gvcfs/${sample_name}.0195.g.vcf.gz", "split_gvcfs/${sample_name}.0196.g.vcf.gz", "split_gvcfs/${sample_name}.0197.g.vcf.gz", "split_gvcfs/${sample_name}.0198.g.vcf.gz", "split_gvcfs/${sample_name}.0199.g.vcf.gz", "split_gvcfs/${sample_name}.0200.g.vcf.gz", "split_gvcfs/${sample_name}.0201.g.vcf.gz", "split_gvcfs/${sample_name}.0202.g.vcf.gz", "split_gvcfs/${sample_name}.0203.g.vcf.gz", "split_gvcfs/${sample_name}.0204.g.vcf.gz", "split_gvcfs/${sample_name}.0205.g.vcf.gz", "split_gvcfs/${sample_name}.0206.g.vcf.gz", "split_gvcfs/${sample_name}.0207.g.vcf.gz", "split_gvcfs/${sample_name}.0208.g.vcf.gz", "split_gvcfs/${sample_name}.0209.g.vcf.gz", "split_gvcfs/${sample_name}.0210.g.vcf.gz", "split_gvcfs/${sample_name}.0211.g.vcf.gz", "split_gvcfs/${sample_name}.0212.g.vcf.gz", "split_gvcfs/${sample_name}.0213.g.vcf.gz", "split_gvcfs/${sample_name}.0214.g.vcf.gz", "split_gvcfs/${sample_name}.0215.g.vcf.gz", "split_gvcfs/${sample_name}.0216.g.vcf.gz", "split_gvcfs/${sample_name}.0217.g.vcf.gz", "split_gvcfs/${sample_name}.0218.g.vcf.gz", "split_gvcfs/${sample_name}.0219.g.vcf.gz", "split_gvcfs/${sample_name}.0220.g.vcf.gz", "split_gvcfs/${sample_name}.0221.g.vcf.gz", "split_gvcfs/${sample_name}.0222.g.vcf.gz", "split_gvcfs/${sample_name}.0223.g.vcf.gz", "split_gvcfs/${sample_name}.0224.g.vcf.gz", "split_gvcfs/${sample_name}.0225.g.vcf.gz", "split_gvcfs/${sample_name}.0226.g.vcf.gz", "split_gvcfs/${sample_name}.0227.g.vcf.gz", "split_gvcfs/${sample_name}.0228.g.vcf.gz", "split_gvcfs/${sample_name}.0229.g.vcf.gz", "split_gvcfs/${sample_name}.0230.g.vcf.gz", "split_gvcfs/${sample_name}.0231.g.vcf.gz", "split_gvcfs/${sample_name}.0232.g.vcf.gz", "split_gvcfs/${sample_name}.0233.g.vcf.gz", "split_gvcfs/${sample_name}.0234.g.vcf.gz", "split_gvcfs/${sample_name}.0235.g.vcf.gz", "split_gvcfs/${sample_name}.0236.g.vcf.gz", "split_gvcfs/${sample_name}.0237.g.vcf.gz", "split_gvcfs/${sample_name}.0238.g.vcf.gz", "split_gvcfs/${sample_name}.0239.g.vcf.gz", "split_gvcfs/${sample_name}.0240.g.vcf.gz", "split_gvcfs/${sample_name}.0241.g.vcf.gz", "split_gvcfs/${sample_name}.0242.g.vcf.gz", "split_gvcfs/${sample_name}.0243.g.vcf.gz", "split_gvcfs/${sample_name}.0244.g.vcf.gz", "split_gvcfs/${sample_name}.0245.g.vcf.gz", "split_gvcfs/${sample_name}.0246.g.vcf.gz", "split_gvcfs/${sample_name}.0247.g.vcf.gz", "split_gvcfs/${sample_name}.0248.g.vcf.gz", "split_gvcfs/${sample_name}.0249.g.vcf.gz", "split_gvcfs/${sample_name}.0250.g.vcf.gz", "split_gvcfs/${sample_name}.0251.g.vcf.gz", "split_gvcfs/${sample_name}.0252.g.vcf.gz", "split_gvcfs/${sample_name}.0253.g.vcf.gz", "split_gvcfs/${sample_name}.0254.g.vcf.gz", "split_gvcfs/${sample_name}.0255.g.vcf.gz", "split_gvcfs/${sample_name}.0256.g.vcf.gz", "split_gvcfs/${sample_name}.0257.g.vcf.gz", "split_gvcfs/${sample_name}.0258.g.vcf.gz", "split_gvcfs/${sample_name}.0259.g.vcf.gz", "split_gvcfs/${sample_name}.0260.g.vcf.gz", "split_gvcfs/${sample_name}.0261.g.vcf.gz", "split_gvcfs/${sample_name}.0262.g.vcf.gz", "split_gvcfs/${sample_name}.0263.g.vcf.gz", "split_gvcfs/${sample_name}.0264.g.vcf.gz", "split_gvcfs/${sample_name}.0265.g.vcf.gz", "split_gvcfs/${sample_name}.0266.g.vcf.gz", "split_gvcfs/${sample_name}.0267.g.vcf.gz", "split_gvcfs/${sample_name}.0268.g.vcf.gz", "split_gvcfs/${sample_name}.0269.g.vcf.gz", "split_gvcfs/${sample_name}.0270.g.vcf.gz", "split_gvcfs/${sample_name}.0271.g.vcf.gz", "split_gvcfs/${sample_name}.0272.g.vcf.gz", "split_gvcfs/${sample_name}.0273.g.vcf.gz", "split_gvcfs/${sample_name}.0274.g.vcf.gz", "split_gvcfs/${sample_name}.0275.g.vcf.gz", "split_gvcfs/${sample_name}.0276.g.vcf.gz", "split_gvcfs/${sample_name}.0277.g.vcf.gz", "split_gvcfs/${sample_name}.0278.g.vcf.gz", "split_gvcfs/${sample_name}.0279.g.vcf.gz", "split_gvcfs/${sample_name}.0280.g.vcf.gz", "split_gvcfs/${sample_name}.0281.g.vcf.gz", "split_gvcfs/${sample_name}.0282.g.vcf.gz", "split_gvcfs/${sample_name}.0283.g.vcf.gz", "split_gvcfs/${sample_name}.0284.g.vcf.gz", "split_gvcfs/${sample_name}.0285.g.vcf.gz", "split_gvcfs/${sample_name}.0286.g.vcf.gz", "split_gvcfs/${sample_name}.0287.g.vcf.gz", "split_gvcfs/${sample_name}.0288.g.vcf.gz", "split_gvcfs/${sample_name}.0289.g.vcf.gz", "split_gvcfs/${sample_name}.0290.g.vcf.gz", "split_gvcfs/${sample_name}.0291.g.vcf.gz", "split_gvcfs/${sample_name}.0292.g.vcf.gz", "split_gvcfs/${sample_name}.0293.g.vcf.gz", "split_gvcfs/${sample_name}.0294.g.vcf.gz", "split_gvcfs/${sample_name}.0295.g.vcf.gz", "split_gvcfs/${sample_name}.0296.g.vcf.gz", "split_gvcfs/${sample_name}.0297.g.vcf.gz", "split_gvcfs/${sample_name}.0298.g.vcf.gz", "split_gvcfs/${sample_name}.0299.g.vcf.gz", "split_gvcfs/${sample_name}.0300.g.vcf.gz", "split_gvcfs/${sample_name}.0301.g.vcf.gz", "split_gvcfs/${sample_name}.0302.g.vcf.gz", "split_gvcfs/${sample_name}.0303.g.vcf.gz", "split_gvcfs/${sample_name}.0304.g.vcf.gz", "split_gvcfs/${sample_name}.0305.g.vcf.gz", "split_gvcfs/${sample_name}.0306.g.vcf.gz", "split_gvcfs/${sample_name}.0307.g.vcf.gz", "split_gvcfs/${sample_name}.0308.g.vcf.gz", "split_gvcfs/${sample_name}.0309.g.vcf.gz", "split_gvcfs/${sample_name}.0310.g.vcf.gz", "split_gvcfs/${sample_name}.0311.g.vcf.gz", "split_gvcfs/${sample_name}.0312.g.vcf.gz", "split_gvcfs/${sample_name}.0313.g.vcf.gz", "split_gvcfs/${sample_name}.0314.g.vcf.gz", "split_gvcfs/${sample_name}.0315.g.vcf.gz", "split_gvcfs/${sample_name}.0316.g.vcf.gz", "split_gvcfs/${sample_name}.0317.g.vcf.gz", "split_gvcfs/${sample_name}.0318.g.vcf.gz", "split_gvcfs/${sample_name}.0319.g.vcf.gz", "split_gvcfs/${sample_name}.0320.g.vcf.gz", "split_gvcfs/${sample_name}.0321.g.vcf.gz", "split_gvcfs/${sample_name}.0322.g.vcf.gz", "split_gvcfs/${sample_name}.0323.g.vcf.gz", "split_gvcfs/${sample_name}.0324.g.vcf.gz", "split_gvcfs/${sample_name}.0325.g.vcf.gz", "split_gvcfs/${sample_name}.0326.g.vcf.gz", "split_gvcfs/${sample_name}.0327.g.vcf.gz", "split_gvcfs/${sample_name}.0328.g.vcf.gz", "split_gvcfs/${sample_name}.0329.g.vcf.gz", "split_gvcfs/${sample_name}.0330.g.vcf.gz", "split_gvcfs/${sample_name}.0331.g.vcf.gz", "split_gvcfs/${sample_name}.0332.g.vcf.gz", "split_gvcfs/${sample_name}.0333.g.vcf.gz", "split_gvcfs/${sample_name}.0334.g.vcf.gz", "split_gvcfs/${sample_name}.0335.g.vcf.gz", "split_gvcfs/${sample_name}.0336.g.vcf.gz", "split_gvcfs/${sample_name}.0337.g.vcf.gz", "split_gvcfs/${sample_name}.0338.g.vcf.gz", "split_gvcfs/${sample_name}.0339.g.vcf.gz", "split_gvcfs/${sample_name}.0340.g.vcf.gz", "split_gvcfs/${sample_name}.0341.g.vcf.gz", "split_gvcfs/${sample_name}.0342.g.vcf.gz", "split_gvcfs/${sample_name}.0343.g.vcf.gz", "split_gvcfs/${sample_name}.0344.g.vcf.gz", "split_gvcfs/${sample_name}.0345.g.vcf.gz", "split_gvcfs/${sample_name}.0346.g.vcf.gz", "split_gvcfs/${sample_name}.0347.g.vcf.gz", "split_gvcfs/${sample_name}.0348.g.vcf.gz", "split_gvcfs/${sample_name}.0349.g.vcf.gz", "split_gvcfs/${sample_name}.0350.g.vcf.gz", "split_gvcfs/${sample_name}.0351.g.vcf.gz", "split_gvcfs/${sample_name}.0352.g.vcf.gz", "split_gvcfs/${sample_name}.0353.g.vcf.gz", "split_gvcfs/${sample_name}.0354.g.vcf.gz", "split_gvcfs/${sample_name}.0355.g.vcf.gz", "split_gvcfs/${sample_name}.0356.g.vcf.gz", "split_gvcfs/${sample_name}.0357.g.vcf.gz", "split_gvcfs/${sample_name}.0358.g.vcf.gz", "split_gvcfs/${sample_name}.0359.g.vcf.gz", "split_gvcfs/${sample_name}.0360.g.vcf.gz", "split_gvcfs/${sample_name}.0361.g.vcf.gz", "split_gvcfs/${sample_name}.0362.g.vcf.gz", "split_gvcfs/${sample_name}.0363.g.vcf.gz", "split_gvcfs/${sample_name}.0364.g.vcf.gz", "split_gvcfs/${sample_name}.0365.g.vcf.gz", "split_gvcfs/${sample_name}.0366.g.vcf.gz", "split_gvcfs/${sample_name}.0367.g.vcf.gz", "split_gvcfs/${sample_name}.0368.g.vcf.gz", "split_gvcfs/${sample_name}.0369.g.vcf.gz", "split_gvcfs/${sample_name}.0370.g.vcf.gz", "split_gvcfs/${sample_name}.0371.g.vcf.gz", "split_gvcfs/${sample_name}.0372.g.vcf.gz", "split_gvcfs/${sample_name}.0373.g.vcf.gz", "split_gvcfs/${sample_name}.0374.g.vcf.gz", "split_gvcfs/${sample_name}.0375.g.vcf.gz", "split_gvcfs/${sample_name}.0376.g.vcf.gz", "split_gvcfs/${sample_name}.0377.g.vcf.gz", "split_gvcfs/${sample_name}.0378.g.vcf.gz", "split_gvcfs/${sample_name}.0379.g.vcf.gz", "split_gvcfs/${sample_name}.0380.g.vcf.gz", "split_gvcfs/${sample_name}.0381.g.vcf.gz", "split_gvcfs/${sample_name}.0382.g.vcf.gz", "split_gvcfs/${sample_name}.0383.g.vcf.gz", "split_gvcfs/${sample_name}.0384.g.vcf.gz", "split_gvcfs/${sample_name}.0385.g.vcf.gz", "split_gvcfs/${sample_name}.0386.g.vcf.gz", "split_gvcfs/${sample_name}.0387.g.vcf.gz", "split_gvcfs/${sample_name}.0388.g.vcf.gz", "split_gvcfs/${sample_name}.0389.g.vcf.gz", "split_gvcfs/${sample_name}.0390.g.vcf.gz", "split_gvcfs/${sample_name}.0391.g.vcf.gz", "split_gvcfs/${sample_name}.0392.g.vcf.gz", "split_gvcfs/${sample_name}.0393.g.vcf.gz", "split_gvcfs/${sample_name}.0394.g.vcf.gz", "split_gvcfs/${sample_name}.0395.g.vcf.gz", "split_gvcfs/${sample_name}.0396.g.vcf.gz", "split_gvcfs/${sample_name}.0397.g.vcf.gz", "split_gvcfs/${sample_name}.0398.g.vcf.gz", "split_gvcfs/${sample_name}.0399.g.vcf.gz", "split_gvcfs/${sample_name}.0400.g.vcf.gz", "split_gvcfs/${sample_name}.0401.g.vcf.gz", "split_gvcfs/${sample_name}.0402.g.vcf.gz", "split_gvcfs/${sample_name}.0403.g.vcf.gz", "split_gvcfs/${sample_name}.0404.g.vcf.gz", "split_gvcfs/${sample_name}.0405.g.vcf.gz", "split_gvcfs/${sample_name}.0406.g.vcf.gz", "split_gvcfs/${sample_name}.0407.g.vcf.gz", "split_gvcfs/${sample_name}.0408.g.vcf.gz", "split_gvcfs/${sample_name}.0409.g.vcf.gz", "split_gvcfs/${sample_name}.0410.g.vcf.gz", "split_gvcfs/${sample_name}.0411.g.vcf.gz", "split_gvcfs/${sample_name}.0412.g.vcf.gz", "split_gvcfs/${sample_name}.0413.g.vcf.gz", "split_gvcfs/${sample_name}.0414.g.vcf.gz", "split_gvcfs/${sample_name}.0415.g.vcf.gz", "split_gvcfs/${sample_name}.0416.g.vcf.gz", "split_gvcfs/${sample_name}.0417.g.vcf.gz", "split_gvcfs/${sample_name}.0418.g.vcf.gz", "split_gvcfs/${sample_name}.0419.g.vcf.gz", "split_gvcfs/${sample_name}.0420.g.vcf.gz", "split_gvcfs/${sample_name}.0421.g.vcf.gz", "split_gvcfs/${sample_name}.0422.g.vcf.gz", "split_gvcfs/${sample_name}.0423.g.vcf.gz", "split_gvcfs/${sample_name}.0424.g.vcf.gz", "split_gvcfs/${sample_name}.0425.g.vcf.gz", "split_gvcfs/${sample_name}.0426.g.vcf.gz", "split_gvcfs/${sample_name}.0427.g.vcf.gz", "split_gvcfs/${sample_name}.0428.g.vcf.gz", "split_gvcfs/${sample_name}.0429.g.vcf.gz", "split_gvcfs/${sample_name}.0430.g.vcf.gz", "split_gvcfs/${sample_name}.0431.g.vcf.gz", "split_gvcfs/${sample_name}.0432.g.vcf.gz", "split_gvcfs/${sample_name}.0433.g.vcf.gz", "split_gvcfs/${sample_name}.0434.g.vcf.gz", "split_gvcfs/${sample_name}.0435.g.vcf.gz", "split_gvcfs/${sample_name}.0436.g.vcf.gz", "split_gvcfs/${sample_name}.0437.g.vcf.gz", "split_gvcfs/${sample_name}.0438.g.vcf.gz", "split_gvcfs/${sample_name}.0439.g.vcf.gz", "split_gvcfs/${sample_name}.0440.g.vcf.gz", "split_gvcfs/${sample_name}.0441.g.vcf.gz", "split_gvcfs/${sample_name}.0442.g.vcf.gz", "split_gvcfs/${sample_name}.0443.g.vcf.gz", "split_gvcfs/${sample_name}.0444.g.vcf.gz", "split_gvcfs/${sample_name}.0445.g.vcf.gz", "split_gvcfs/${sample_name}.0446.g.vcf.gz", "split_gvcfs/${sample_name}.0447.g.vcf.gz", "split_gvcfs/${sample_name}.0448.g.vcf.gz", "split_gvcfs/${sample_name}.0449.g.vcf.gz", "split_gvcfs/${sample_name}.0450.g.vcf.gz", "split_gvcfs/${sample_name}.0451.g.vcf.gz", "split_gvcfs/${sample_name}.0452.g.vcf.gz", "split_gvcfs/${sample_name}.0453.g.vcf.gz", "split_gvcfs/${sample_name}.0454.g.vcf.gz", "split_gvcfs/${sample_name}.0455.g.vcf.gz", "split_gvcfs/${sample_name}.0456.g.vcf.gz", "split_gvcfs/${sample_name}.0457.g.vcf.gz", "split_gvcfs/${sample_name}.0458.g.vcf.gz", "split_gvcfs/${sample_name}.0459.g.vcf.gz", "split_gvcfs/${sample_name}.0460.g.vcf.gz", "split_gvcfs/${sample_name}.0461.g.vcf.gz", "split_gvcfs/${sample_name}.0462.g.vcf.gz", "split_gvcfs/${sample_name}.0463.g.vcf.gz", "split_gvcfs/${sample_name}.0464.g.vcf.gz", "split_gvcfs/${sample_name}.0465.g.vcf.gz", "split_gvcfs/${sample_name}.0466.g.vcf.gz", "split_gvcfs/${sample_name}.0467.g.vcf.gz", "split_gvcfs/${sample_name}.0468.g.vcf.gz", "split_gvcfs/${sample_name}.0469.g.vcf.gz", "split_gvcfs/${sample_name}.0470.g.vcf.gz", "split_gvcfs/${sample_name}.0471.g.vcf.gz", "split_gvcfs/${sample_name}.0472.g.vcf.gz", "split_gvcfs/${sample_name}.0473.g.vcf.gz", "split_gvcfs/${sample_name}.0474.g.vcf.gz", "split_gvcfs/${sample_name}.0475.g.vcf.gz", "split_gvcfs/${sample_name}.0476.g.vcf.gz", "split_gvcfs/${sample_name}.0477.g.vcf.gz", "split_gvcfs/${sample_name}.0478.g.vcf.gz", "split_gvcfs/${sample_name}.0479.g.vcf.gz", "split_gvcfs/${sample_name}.0480.g.vcf.gz", "split_gvcfs/${sample_name}.0481.g.vcf.gz", "split_gvcfs/${sample_name}.0482.g.vcf.gz", "split_gvcfs/${sample_name}.0483.g.vcf.gz", "split_gvcfs/${sample_name}.0484.g.vcf.gz", "split_gvcfs/${sample_name}.0485.g.vcf.gz", "split_gvcfs/${sample_name}.0486.g.vcf.gz", "split_gvcfs/${sample_name}.0487.g.vcf.gz", "split_gvcfs/${sample_name}.0488.g.vcf.gz", "split_gvcfs/${sample_name}.0489.g.vcf.gz", "split_gvcfs/${sample_name}.0490.g.vcf.gz", "split_gvcfs/${sample_name}.0491.g.vcf.gz", "split_gvcfs/${sample_name}.0492.g.vcf.gz", "split_gvcfs/${sample_name}.0493.g.vcf.gz", "split_gvcfs/${sample_name}.0494.g.vcf.gz", "split_gvcfs/${sample_name}.0495.g.vcf.gz", "split_gvcfs/${sample_name}.0496.g.vcf.gz", "split_gvcfs/${sample_name}.0497.g.vcf.gz", "split_gvcfs/${sample_name}.0498.g.vcf.gz", "split_gvcfs/${sample_name}.0499.g.vcf.gz", "split_gvcfs/${sample_name}.0500.g.vcf.gz", "split_gvcfs/${sample_name}.0501.g.vcf.gz", "split_gvcfs/${sample_name}.0502.g.vcf.gz", "split_gvcfs/${sample_name}.0503.g.vcf.gz", "split_gvcfs/${sample_name}.0504.g.vcf.gz", "split_gvcfs/${sample_name}.0505.g.vcf.gz", "split_gvcfs/${sample_name}.0506.g.vcf.gz", "split_gvcfs/${sample_name}.0507.g.vcf.gz", "split_gvcfs/${sample_name}.0508.g.vcf.gz", "split_gvcfs/${sample_name}.0509.g.vcf.gz", "split_gvcfs/${sample_name}.0510.g.vcf.gz", "split_gvcfs/${sample_name}.0511.g.vcf.gz", "split_gvcfs/${sample_name}.0512.g.vcf.gz", "split_gvcfs/${sample_name}.0513.g.vcf.gz", "split_gvcfs/${sample_name}.0514.g.vcf.gz", "split_gvcfs/${sample_name}.0515.g.vcf.gz", "split_gvcfs/${sample_name}.0516.g.vcf.gz", "split_gvcfs/${sample_name}.0517.g.vcf.gz", "split_gvcfs/${sample_name}.0518.g.vcf.gz", "split_gvcfs/${sample_name}.0519.g.vcf.gz", "split_gvcfs/${sample_name}.0520.g.vcf.gz", "split_gvcfs/${sample_name}.0521.g.vcf.gz", "split_gvcfs/${sample_name}.0522.g.vcf.gz", "split_gvcfs/${sample_name}.0523.g.vcf.gz", "split_gvcfs/${sample_name}.0524.g.vcf.gz", "split_gvcfs/${sample_name}.0525.g.vcf.gz", "split_gvcfs/${sample_name}.0526.g.vcf.gz", "split_gvcfs/${sample_name}.0527.g.vcf.gz", "split_gvcfs/${sample_name}.0528.g.vcf.gz", "split_gvcfs/${sample_name}.0529.g.vcf.gz", "split_gvcfs/${sample_name}.0530.g.vcf.gz", "split_gvcfs/${sample_name}.0531.g.vcf.gz", "split_gvcfs/${sample_name}.0532.g.vcf.gz", "split_gvcfs/${sample_name}.0533.g.vcf.gz", "split_gvcfs/${sample_name}.0534.g.vcf.gz", "split_gvcfs/${sample_name}.0535.g.vcf.gz", "split_gvcfs/${sample_name}.0536.g.vcf.gz", "split_gvcfs/${sample_name}.0537.g.vcf.gz", "split_gvcfs/${sample_name}.0538.g.vcf.gz", "split_gvcfs/${sample_name}.0539.g.vcf.gz", "split_gvcfs/${sample_name}.0540.g.vcf.gz", "split_gvcfs/${sample_name}.0541.g.vcf.gz", "split_gvcfs/${sample_name}.0542.g.vcf.gz", "split_gvcfs/${sample_name}.0543.g.vcf.gz", "split_gvcfs/${sample_name}.0544.g.vcf.gz", "split_gvcfs/${sample_name}.0545.g.vcf.gz", "split_gvcfs/${sample_name}.0546.g.vcf.gz", "split_gvcfs/${sample_name}.0547.g.vcf.gz", "split_gvcfs/${sample_name}.0548.g.vcf.gz", "split_gvcfs/${sample_name}.0549.g.vcf.gz", "split_gvcfs/${sample_name}.0550.g.vcf.gz", "split_gvcfs/${sample_name}.0551.g.vcf.gz", "split_gvcfs/${sample_name}.0552.g.vcf.gz", "split_gvcfs/${sample_name}.0553.g.vcf.gz", "split_gvcfs/${sample_name}.0554.g.vcf.gz", "split_gvcfs/${sample_name}.0555.g.vcf.gz", "split_gvcfs/${sample_name}.0556.g.vcf.gz", "split_gvcfs/${sample_name}.0557.g.vcf.gz", "split_gvcfs/${sample_name}.0558.g.vcf.gz", "split_gvcfs/${sample_name}.0559.g.vcf.gz", "split_gvcfs/${sample_name}.0560.g.vcf.gz", "split_gvcfs/${sample_name}.0561.g.vcf.gz", "split_gvcfs/${sample_name}.0562.g.vcf.gz", "split_gvcfs/${sample_name}.0563.g.vcf.gz", "split_gvcfs/${sample_name}.0564.g.vcf.gz", "split_gvcfs/${sample_name}.0565.g.vcf.gz", "split_gvcfs/${sample_name}.0566.g.vcf.gz", "split_gvcfs/${sample_name}.0567.g.vcf.gz", "split_gvcfs/${sample_name}.0568.g.vcf.gz", "split_gvcfs/${sample_name}.0569.g.vcf.gz", "split_gvcfs/${sample_name}.0570.g.vcf.gz", "split_gvcfs/${sample_name}.0571.g.vcf.gz", "split_gvcfs/${sample_name}.0572.g.vcf.gz", "split_gvcfs/${sample_name}.0573.g.vcf.gz", "split_gvcfs/${sample_name}.0574.g.vcf.gz", "split_gvcfs/${sample_name}.0575.g.vcf.gz", "split_gvcfs/${sample_name}.0576.g.vcf.gz", "split_gvcfs/${sample_name}.0577.g.vcf.gz", "split_gvcfs/${sample_name}.0578.g.vcf.gz", "split_gvcfs/${sample_name}.0579.g.vcf.gz", "split_gvcfs/${sample_name}.0580.g.vcf.gz", "split_gvcfs/${sample_name}.0581.g.vcf.gz", "split_gvcfs/${sample_name}.0582.g.vcf.gz", "split_gvcfs/${sample_name}.0583.g.vcf.gz", "split_gvcfs/${sample_name}.0584.g.vcf.gz", "split_gvcfs/${sample_name}.0585.g.vcf.gz", "split_gvcfs/${sample_name}.0586.g.vcf.gz", "split_gvcfs/${sample_name}.0587.g.vcf.gz", "split_gvcfs/${sample_name}.0588.g.vcf.gz", "split_gvcfs/${sample_name}.0589.g.vcf.gz", "split_gvcfs/${sample_name}.0590.g.vcf.gz", "split_gvcfs/${sample_name}.0591.g.vcf.gz", "split_gvcfs/${sample_name}.0592.g.vcf.gz", "split_gvcfs/${sample_name}.0593.g.vcf.gz", "split_gvcfs/${sample_name}.0594.g.vcf.gz", "split_gvcfs/${sample_name}.0595.g.vcf.gz", "split_gvcfs/${sample_name}.0596.g.vcf.gz", "split_gvcfs/${sample_name}.0597.g.vcf.gz", "split_gvcfs/${sample_name}.0598.g.vcf.gz", "split_gvcfs/${sample_name}.0599.g.vcf.gz", "split_gvcfs/${sample_name}.0600.g.vcf.gz", "split_gvcfs/${sample_name}.0601.g.vcf.gz", "split_gvcfs/${sample_name}.0602.g.vcf.gz", "split_gvcfs/${sample_name}.0603.g.vcf.gz", "split_gvcfs/${sample_name}.0604.g.vcf.gz", "split_gvcfs/${sample_name}.0605.g.vcf.gz", "split_gvcfs/${sample_name}.0606.g.vcf.gz", "split_gvcfs/${sample_name}.0607.g.vcf.gz", "split_gvcfs/${sample_name}.0608.g.vcf.gz", "split_gvcfs/${sample_name}.0609.g.vcf.gz", "split_gvcfs/${sample_name}.0610.g.vcf.gz", "split_gvcfs/${sample_name}.0611.g.vcf.gz", "split_gvcfs/${sample_name}.0612.g.vcf.gz", "split_gvcfs/${sample_name}.0613.g.vcf.gz", "split_gvcfs/${sample_name}.0614.g.vcf.gz", "split_gvcfs/${sample_name}.0615.g.vcf.gz", "split_gvcfs/${sample_name}.0616.g.vcf.gz", "split_gvcfs/${sample_name}.0617.g.vcf.gz", "split_gvcfs/${sample_name}.0618.g.vcf.gz", "split_gvcfs/${sample_name}.0619.g.vcf.gz", "split_gvcfs/${sample_name}.0620.g.vcf.gz", "split_gvcfs/${sample_name}.0621.g.vcf.gz", "split_gvcfs/${sample_name}.0622.g.vcf.gz", "split_gvcfs/${sample_name}.0623.g.vcf.gz", "split_gvcfs/${sample_name}.0624.g.vcf.gz", "split_gvcfs/${sample_name}.0625.g.vcf.gz", "split_gvcfs/${sample_name}.0626.g.vcf.gz", "split_gvcfs/${sample_name}.0627.g.vcf.gz", "split_gvcfs/${sample_name}.0628.g.vcf.gz", "split_gvcfs/${sample_name}.0629.g.vcf.gz", "split_gvcfs/${sample_name}.0630.g.vcf.gz", "split_gvcfs/${sample_name}.0631.g.vcf.gz", "split_gvcfs/${sample_name}.0632.g.vcf.gz", "split_gvcfs/${sample_name}.0633.g.vcf.gz", "split_gvcfs/${sample_name}.0634.g.vcf.gz", "split_gvcfs/${sample_name}.0635.g.vcf.gz", "split_gvcfs/${sample_name}.0636.g.vcf.gz", "split_gvcfs/${sample_name}.0637.g.vcf.gz", "split_gvcfs/${sample_name}.0638.g.vcf.gz", "split_gvcfs/${sample_name}.0639.g.vcf.gz", "split_gvcfs/${sample_name}.0640.g.vcf.gz", "split_gvcfs/${sample_name}.0641.g.vcf.gz", "split_gvcfs/${sample_name}.0642.g.vcf.gz", "split_gvcfs/${sample_name}.0643.g.vcf.gz", "split_gvcfs/${sample_name}.0644.g.vcf.gz", "split_gvcfs/${sample_name}.0645.g.vcf.gz", "split_gvcfs/${sample_name}.0646.g.vcf.gz", "split_gvcfs/${sample_name}.0647.g.vcf.gz", "split_gvcfs/${sample_name}.0648.g.vcf.gz", "split_gvcfs/${sample_name}.0649.g.vcf.gz", "split_gvcfs/${sample_name}.0650.g.vcf.gz", "split_gvcfs/${sample_name}.0651.g.vcf.gz", "split_gvcfs/${sample_name}.0652.g.vcf.gz", "split_gvcfs/${sample_name}.0653.g.vcf.gz", "split_gvcfs/${sample_name}.0654.g.vcf.gz", "split_gvcfs/${sample_name}.0655.g.vcf.gz", "split_gvcfs/${sample_name}.0656.g.vcf.gz", "split_gvcfs/${sample_name}.0657.g.vcf.gz", "split_gvcfs/${sample_name}.0658.g.vcf.gz", "split_gvcfs/${sample_name}.0659.g.vcf.gz", "split_gvcfs/${sample_name}.0660.g.vcf.gz", "split_gvcfs/${sample_name}.0661.g.vcf.gz", "split_gvcfs/${sample_name}.0662.g.vcf.gz", "split_gvcfs/${sample_name}.0663.g.vcf.gz", "split_gvcfs/${sample_name}.0664.g.vcf.gz", "split_gvcfs/${sample_name}.0665.g.vcf.gz", "split_gvcfs/${sample_name}.0666.g.vcf.gz", "split_gvcfs/${sample_name}.0667.g.vcf.gz", "split_gvcfs/${sample_name}.0668.g.vcf.gz", "split_gvcfs/${sample_name}.0669.g.vcf.gz", "split_gvcfs/${sample_name}.0670.g.vcf.gz", "split_gvcfs/${sample_name}.0671.g.vcf.gz", "split_gvcfs/${sample_name}.0672.g.vcf.gz", "split_gvcfs/${sample_name}.0673.g.vcf.gz", "split_gvcfs/${sample_name}.0674.g.vcf.gz", "split_gvcfs/${sample_name}.0675.g.vcf.gz", "split_gvcfs/${sample_name}.0676.g.vcf.gz", "split_gvcfs/${sample_name}.0677.g.vcf.gz", "split_gvcfs/${sample_name}.0678.g.vcf.gz", "split_gvcfs/${sample_name}.0679.g.vcf.gz", "split_gvcfs/${sample_name}.0680.g.vcf.gz", "split_gvcfs/${sample_name}.0681.g.vcf.gz", "split_gvcfs/${sample_name}.0682.g.vcf.gz", "split_gvcfs/${sample_name}.0683.g.vcf.gz", "split_gvcfs/${sample_name}.0684.g.vcf.gz", "split_gvcfs/${sample_name}.0685.g.vcf.gz", "split_gvcfs/${sample_name}.0686.g.vcf.gz", "split_gvcfs/${sample_name}.0687.g.vcf.gz", "split_gvcfs/${sample_name}.0688.g.vcf.gz", "split_gvcfs/${sample_name}.0689.g.vcf.gz", "split_gvcfs/${sample_name}.0690.g.vcf.gz", "split_gvcfs/${sample_name}.0691.g.vcf.gz", "split_gvcfs/${sample_name}.0692.g.vcf.gz", "split_gvcfs/${sample_name}.0693.g.vcf.gz", "split_gvcfs/${sample_name}.0694.g.vcf.gz", "split_gvcfs/${sample_name}.0695.g.vcf.gz", "split_gvcfs/${sample_name}.0696.g.vcf.gz", "split_gvcfs/${sample_name}.0697.g.vcf.gz", "split_gvcfs/${sample_name}.0698.g.vcf.gz", "split_gvcfs/${sample_name}.0699.g.vcf.gz", "split_gvcfs/${sample_name}.0700.g.vcf.gz", "split_gvcfs/${sample_name}.0701.g.vcf.gz", "split_gvcfs/${sample_name}.0702.g.vcf.gz", "split_gvcfs/${sample_name}.0703.g.vcf.gz", "split_gvcfs/${sample_name}.0704.g.vcf.gz", "split_gvcfs/${sample_name}.0705.g.vcf.gz", "split_gvcfs/${sample_name}.0706.g.vcf.gz", "split_gvcfs/${sample_name}.0707.g.vcf.gz", "split_gvcfs/${sample_name}.0708.g.vcf.gz", "split_gvcfs/${sample_name}.0709.g.vcf.gz", "split_gvcfs/${sample_name}.0710.g.vcf.gz", "split_gvcfs/${sample_name}.0711.g.vcf.gz", "split_gvcfs/${sample_name}.0712.g.vcf.gz", "split_gvcfs/${sample_name}.0713.g.vcf.gz", "split_gvcfs/${sample_name}.0714.g.vcf.gz", "split_gvcfs/${sample_name}.0715.g.vcf.gz", "split_gvcfs/${sample_name}.0716.g.vcf.gz", "split_gvcfs/${sample_name}.0717.g.vcf.gz", "split_gvcfs/${sample_name}.0718.g.vcf.gz", "split_gvcfs/${sample_name}.0719.g.vcf.gz", "split_gvcfs/${sample_name}.0720.g.vcf.gz", "split_gvcfs/${sample_name}.0721.g.vcf.gz", "split_gvcfs/${sample_name}.0722.g.vcf.gz", "split_gvcfs/${sample_name}.0723.g.vcf.gz", "split_gvcfs/${sample_name}.0724.g.vcf.gz", "split_gvcfs/${sample_name}.0725.g.vcf.gz", "split_gvcfs/${sample_name}.0726.g.vcf.gz", "split_gvcfs/${sample_name}.0727.g.vcf.gz", "split_gvcfs/${sample_name}.0728.g.vcf.gz", "split_gvcfs/${sample_name}.0729.g.vcf.gz", "split_gvcfs/${sample_name}.0730.g.vcf.gz", "split_gvcfs/${sample_name}.0731.g.vcf.gz", "split_gvcfs/${sample_name}.0732.g.vcf.gz", "split_gvcfs/${sample_name}.0733.g.vcf.gz", "split_gvcfs/${sample_name}.0734.g.vcf.gz", "split_gvcfs/${sample_name}.0735.g.vcf.gz", "split_gvcfs/${sample_name}.0736.g.vcf.gz", "split_gvcfs/${sample_name}.0737.g.vcf.gz", "split_gvcfs/${sample_name}.0738.g.vcf.gz", "split_gvcfs/${sample_name}.0739.g.vcf.gz", "split_gvcfs/${sample_name}.0740.g.vcf.gz", "split_gvcfs/${sample_name}.0741.g.vcf.gz", "split_gvcfs/${sample_name}.0742.g.vcf.gz", "split_gvcfs/${sample_name}.0743.g.vcf.gz", "split_gvcfs/${sample_name}.0744.g.vcf.gz", "split_gvcfs/${sample_name}.0745.g.vcf.gz", "split_gvcfs/${sample_name}.0746.g.vcf.gz", "split_gvcfs/${sample_name}.0747.g.vcf.gz", "split_gvcfs/${sample_name}.0748.g.vcf.gz", "split_gvcfs/${sample_name}.0749.g.vcf.gz", "split_gvcfs/${sample_name}.0750.g.vcf.gz", "split_gvcfs/${sample_name}.0751.g.vcf.gz", "split_gvcfs/${sample_name}.0752.g.vcf.gz", "split_gvcfs/${sample_name}.0753.g.vcf.gz", "split_gvcfs/${sample_name}.0754.g.vcf.gz", "split_gvcfs/${sample_name}.0755.g.vcf.gz", "split_gvcfs/${sample_name}.0756.g.vcf.gz", "split_gvcfs/${sample_name}.0757.g.vcf.gz", "split_gvcfs/${sample_name}.0758.g.vcf.gz", "split_gvcfs/${sample_name}.0759.g.vcf.gz", "split_gvcfs/${sample_name}.0760.g.vcf.gz", "split_gvcfs/${sample_name}.0761.g.vcf.gz", "split_gvcfs/${sample_name}.0762.g.vcf.gz", "split_gvcfs/${sample_name}.0763.g.vcf.gz", "split_gvcfs/${sample_name}.0764.g.vcf.gz", "split_gvcfs/${sample_name}.0765.g.vcf.gz", "split_gvcfs/${sample_name}.0766.g.vcf.gz", "split_gvcfs/${sample_name}.0767.g.vcf.gz", "split_gvcfs/${sample_name}.0768.g.vcf.gz", "split_gvcfs/${sample_name}.0769.g.vcf.gz", "split_gvcfs/${sample_name}.0770.g.vcf.gz", "split_gvcfs/${sample_name}.0771.g.vcf.gz", "split_gvcfs/${sample_name}.0772.g.vcf.gz", "split_gvcfs/${sample_name}.0773.g.vcf.gz", "split_gvcfs/${sample_name}.0774.g.vcf.gz", "split_gvcfs/${sample_name}.0775.g.vcf.gz", "split_gvcfs/${sample_name}.0776.g.vcf.gz", "split_gvcfs/${sample_name}.0777.g.vcf.gz", "split_gvcfs/${sample_name}.0778.g.vcf.gz", "split_gvcfs/${sample_name}.0779.g.vcf.gz", "split_gvcfs/${sample_name}.0780.g.vcf.gz", "split_gvcfs/${sample_name}.0781.g.vcf.gz", "split_gvcfs/${sample_name}.0782.g.vcf.gz", "split_gvcfs/${sample_name}.0783.g.vcf.gz", "split_gvcfs/${sample_name}.0784.g.vcf.gz", "split_gvcfs/${sample_name}.0785.g.vcf.gz", "split_gvcfs/${sample_name}.0786.g.vcf.gz", "split_gvcfs/${sample_name}.0787.g.vcf.gz", "split_gvcfs/${sample_name}.0788.g.vcf.gz", "split_gvcfs/${sample_name}.0789.g.vcf.gz", "split_gvcfs/${sample_name}.0790.g.vcf.gz", "split_gvcfs/${sample_name}.0791.g.vcf.gz", "split_gvcfs/${sample_name}.0792.g.vcf.gz", "split_gvcfs/${sample_name}.0793.g.vcf.gz", "split_gvcfs/${sample_name}.0794.g.vcf.gz", "split_gvcfs/${sample_name}.0795.g.vcf.gz", "split_gvcfs/${sample_name}.0796.g.vcf.gz", "split_gvcfs/${sample_name}.0797.g.vcf.gz", "split_gvcfs/${sample_name}.0798.g.vcf.gz", "split_gvcfs/${sample_name}.0799.g.vcf.gz", "split_gvcfs/${sample_name}.0800.g.vcf.gz", "split_gvcfs/${sample_name}.0801.g.vcf.gz", "split_gvcfs/${sample_name}.0802.g.vcf.gz", "split_gvcfs/${sample_name}.0803.g.vcf.gz", "split_gvcfs/${sample_name}.0804.g.vcf.gz", "split_gvcfs/${sample_name}.0805.g.vcf.gz", "split_gvcfs/${sample_name}.0806.g.vcf.gz", "split_gvcfs/${sample_name}.0807.g.vcf.gz", "split_gvcfs/${sample_name}.0808.g.vcf.gz", "split_gvcfs/${sample_name}.0809.g.vcf.gz", "split_gvcfs/${sample_name}.0810.g.vcf.gz", "split_gvcfs/${sample_name}.0811.g.vcf.gz", "split_gvcfs/${sample_name}.0812.g.vcf.gz", "split_gvcfs/${sample_name}.0813.g.vcf.gz", "split_gvcfs/${sample_name}.0814.g.vcf.gz", "split_gvcfs/${sample_name}.0815.g.vcf.gz", "split_gvcfs/${sample_name}.0816.g.vcf.gz", "split_gvcfs/${sample_name}.0817.g.vcf.gz", "split_gvcfs/${sample_name}.0818.g.vcf.gz", "split_gvcfs/${sample_name}.0819.g.vcf.gz", "split_gvcfs/${sample_name}.0820.g.vcf.gz", "split_gvcfs/${sample_name}.0821.g.vcf.gz", "split_gvcfs/${sample_name}.0822.g.vcf.gz", "split_gvcfs/${sample_name}.0823.g.vcf.gz", "split_gvcfs/${sample_name}.0824.g.vcf.gz", "split_gvcfs/${sample_name}.0825.g.vcf.gz", "split_gvcfs/${sample_name}.0826.g.vcf.gz", "split_gvcfs/${sample_name}.0827.g.vcf.gz", "split_gvcfs/${sample_name}.0828.g.vcf.gz", "split_gvcfs/${sample_name}.0829.g.vcf.gz", "split_gvcfs/${sample_name}.0830.g.vcf.gz", "split_gvcfs/${sample_name}.0831.g.vcf.gz", "split_gvcfs/${sample_name}.0832.g.vcf.gz", "split_gvcfs/${sample_name}.0833.g.vcf.gz", "split_gvcfs/${sample_name}.0834.g.vcf.gz", "split_gvcfs/${sample_name}.0835.g.vcf.gz", "split_gvcfs/${sample_name}.0836.g.vcf.gz", "split_gvcfs/${sample_name}.0837.g.vcf.gz", "split_gvcfs/${sample_name}.0838.g.vcf.gz", "split_gvcfs/${sample_name}.0839.g.vcf.gz", "split_gvcfs/${sample_name}.0840.g.vcf.gz", "split_gvcfs/${sample_name}.0841.g.vcf.gz", "split_gvcfs/${sample_name}.0842.g.vcf.gz", "split_gvcfs/${sample_name}.0843.g.vcf.gz", "split_gvcfs/${sample_name}.0844.g.vcf.gz", "split_gvcfs/${sample_name}.0845.g.vcf.gz", "split_gvcfs/${sample_name}.0846.g.vcf.gz", "split_gvcfs/${sample_name}.0847.g.vcf.gz", "split_gvcfs/${sample_name}.0848.g.vcf.gz", "split_gvcfs/${sample_name}.0849.g.vcf.gz", "split_gvcfs/${sample_name}.0850.g.vcf.gz", "split_gvcfs/${sample_name}.0851.g.vcf.gz", "split_gvcfs/${sample_name}.0852.g.vcf.gz", "split_gvcfs/${sample_name}.0853.g.vcf.gz", "split_gvcfs/${sample_name}.0854.g.vcf.gz", "split_gvcfs/${sample_name}.0855.g.vcf.gz", "split_gvcfs/${sample_name}.0856.g.vcf.gz", "split_gvcfs/${sample_name}.0857.g.vcf.gz", "split_gvcfs/${sample_name}.0858.g.vcf.gz", "split_gvcfs/${sample_name}.0859.g.vcf.gz", "split_gvcfs/${sample_name}.0860.g.vcf.gz", "split_gvcfs/${sample_name}.0861.g.vcf.gz", "split_gvcfs/${sample_name}.0862.g.vcf.gz", "split_gvcfs/${sample_name}.0863.g.vcf.gz", "split_gvcfs/${sample_name}.0864.g.vcf.gz", "split_gvcfs/${sample_name}.0865.g.vcf.gz", "split_gvcfs/${sample_name}.0866.g.vcf.gz", "split_gvcfs/${sample_name}.0867.g.vcf.gz", "split_gvcfs/${sample_name}.0868.g.vcf.gz", "split_gvcfs/${sample_name}.0869.g.vcf.gz", "split_gvcfs/${sample_name}.0870.g.vcf.gz", "split_gvcfs/${sample_name}.0871.g.vcf.gz", "split_gvcfs/${sample_name}.0872.g.vcf.gz", "split_gvcfs/${sample_name}.0873.g.vcf.gz", "split_gvcfs/${sample_name}.0874.g.vcf.gz", "split_gvcfs/${sample_name}.0875.g.vcf.gz", "split_gvcfs/${sample_name}.0876.g.vcf.gz", "split_gvcfs/${sample_name}.0877.g.vcf.gz", "split_gvcfs/${sample_name}.0878.g.vcf.gz", "split_gvcfs/${sample_name}.0879.g.vcf.gz", "split_gvcfs/${sample_name}.0880.g.vcf.gz", "split_gvcfs/${sample_name}.0881.g.vcf.gz", "split_gvcfs/${sample_name}.0882.g.vcf.gz", "split_gvcfs/${sample_name}.0883.g.vcf.gz", "split_gvcfs/${sample_name}.0884.g.vcf.gz", "split_gvcfs/${sample_name}.0885.g.vcf.gz", "split_gvcfs/${sample_name}.0886.g.vcf.gz", "split_gvcfs/${sample_name}.0887.g.vcf.gz", "split_gvcfs/${sample_name}.0888.g.vcf.gz", "split_gvcfs/${sample_name}.0889.g.vcf.gz", "split_gvcfs/${sample_name}.0890.g.vcf.gz", "split_gvcfs/${sample_name}.0891.g.vcf.gz", "split_gvcfs/${sample_name}.0892.g.vcf.gz", "split_gvcfs/${sample_name}.0893.g.vcf.gz", "split_gvcfs/${sample_name}.0894.g.vcf.gz", "split_gvcfs/${sample_name}.0895.g.vcf.gz", "split_gvcfs/${sample_name}.0896.g.vcf.gz", "split_gvcfs/${sample_name}.0897.g.vcf.gz", "split_gvcfs/${sample_name}.0898.g.vcf.gz", "split_gvcfs/${sample_name}.0899.g.vcf.gz", "split_gvcfs/${sample_name}.0900.g.vcf.gz"]

  }
}

task MatrixRotation {
  Array[Array[String]] input_matrix

    command <<<
      python <<CODE
      import csv
      with open("${write_tsv(input_matrix)}") as tsv_in:
        input_matrix = [line.strip().split("\t") for line in tsv_in]
        final_matrix = [["" for x in range(len(input_matrix))] for y in range(len(input_matrix[0]))]
        for x in range(len(input_matrix)):
          for y in range(len(input_matrix[0])):
            final_matrix[y][x] = input_matrix[x][y]
      with open("output.tsv", "w") as tsv_out:
        # lineterminator is a workaround for a cromwell bug that doesnt allow for "\r\n" line endings which this outputs by default
        writer = csv.writer(tsv_out, delimiter="\t", lineterminator="\n")
        writer.writerows(final_matrix)
      CODE
      >>>

    runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }

    output {
        File output_tsv = "output.tsv"
    }
}

task makefofns {
    Array[Array[String]] tsv

    command <<<
      python <<CODE
      i=0
      for line in tsv:
        for l in line:
          with open("fofo." + str(i) , "w+"):
            file.write(l + "\n")
        i+=1
        CODE
    >>>

    output {
        Array[File] outp = glob("fofo.*")
    }
}

task ArrayToFile {
    Array[String] array

    command {
    }

    runtime {
      docker: "python:2.7"
      memory: "1 GB"
    }

    output {
        File write_output = write_lines(array)
    }
}

workflow JointGenotyping {
  Array[String] gvcf_string_array
 
  Array[String] samples

  File split_interval_list
  Int small_disk


  scatter (gvcf in gvcf_string_array) {

    call GetUUID 

    call UniquifyGvcf {
      input: 
      file=gvcf,
      uuid = GetUUID.uuid,
      disk_size = small_disk
    }

    call SplitGvcf {
      input:
        gvcf = UniquifyGvcf.output_gvcf,
        sample_name = sub(sub(UniquifyGvcf.output_gvcf, "gs://.*/",""), ".g.vcf.gz", ""),
        interval_list = split_interval_list,
        disk_size = small_disk
    }
  }

  call MakeSampleMap {
    input:
      gvcfs=UniquifyGvcf.output_gvcf,
      samples=samples
  }

  call MatrixRotation as RotateGVCF {
    input:
      input_matrix = SplitGvcf.gvcf_list
  }

  output {
    RotateGVCF.output_tsv
    MakeSampleMap.sample_map
  }
}
