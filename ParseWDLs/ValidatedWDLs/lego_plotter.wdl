task lego_plotter_task {

#maf file for plotting
File maf
#reference data
File categs
#label name
String labName
#parameters for the platter
String covString
# 'exome', 'genome', or 'unit' to pick typical exome coverage, genome coverage, or "unit" as every base counts equally in the normalization. 


#runtime parameter for memory/disk
String memGB
String diskGB

command <<<

#increase verbosity to inform of commands run
set -x

#create directory where category data is expected
mkdir -pv /xchip/cga/reference/hg19/annotation/db/tracks/hg19/c65e/

#copy the category data
cp -vf ${categs} /xchip/cga/reference/hg19/annotation/db/tracks/hg19/c65e/

#link for the name
ln -vs ${maf} ${labName}.maf

#run the plotter and acquire the exit code (plot rates and zscale)
/usr/local/bin/run_call_lego_plotter.sh ${labName}.maf  . ${covString} yes yes
PLOTTER_EXIT_CODE=$? ;

#encode the PNGS
for PNG in `find *.png`; do 
	uuencode  -m $PNG  /dev/stdout | grep -Pv '^begin'|grep -Pv '^====$'|tr -d "\n" > $PNG.enc ; 
	UUE=$? ; 
	if [ "$UUE" -ne "0" ] ; then exit 2 ; fi ;
	done 

#update the HTML in memory/place
python_cmd="
import glob
import re
for H in glob.glob('*_legos.html'):
	reader=open(H,'r')
	#store the HTML into a single string
	all_lines=''
	for line in reader:
		all_lines=all_lines+line
	reader.close()
	#iterate over the encs
	for ENC in glob.glob('*.enc'):
		print 'got enc = ',ENC
		PNG=re.sub('\.enc$','',ENC)
		enc_reader=open(ENC,'r')
		enc_data=''
		for enc_line in enc_reader:
			enc_data=enc_data+enc_line.strip()
		enc_reader.close()
		#embed the img, close the src, and add a title
		enc_replacement='data:image/png;base64,'+str(enc_data)+'\" title=\"'+str(PNG)+'\" '
		all_lines=all_lines.replace(PNG,enc_replacement)
	writer=open(H,'w')
	writer.write(all_lines)
	writer.close()
"
python -c "$python_cmd" 
PY_EXIT=$? ;
if [ "$PY_EXIT" -ne "0" ] ; then exit 3 ; fi ;


#propagate exit code
bash -c "exit $PLOTTER_EXIT_CODE" 

>>>

runtime {
	memory: "${memGB} GB"
	cpu : "1"
	docker: "broadinstitute/lego_plotter@sha256:d8cc0eb021aed63134619d3f9c3169436484579476940b85f93ac3df71f82908"
	disks: "local-disk ${diskGB} HDD"
	}

output {
	Array[File] ais=glob("*.ai")
	Array[File] pngs=glob("*.png")
	Array[File] figs=glob("*.fig")
	Array[File] pss=glob("*.ps")
	File mut_legos_html=		"${labName}.maf.mutation_legos.html"
	}
}

workflow lego_plotter_workflow {
 
    call lego_plotter_task

	  meta {
	    author: "Eddie Salinas"
	    email: "esalinas@broadinstitute.org"
	  }

	output {
		lego_plotter_task.mut_legos_html		
		}

	}