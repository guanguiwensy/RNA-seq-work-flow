#!/bin/bash

#The function of help document.
func(){
    echo "Usage:"
	echo "RNA-seq.sh [-t threads] [-f directory of reference of hisat2] [-g directory of gtf] [-n derectory of feature of fastq files' name]"
    echo "Example : for 'aa.R1.fastq.gz aa.R2.fastq.gz'. The -n could be aa.R"
    exit 1
}

#Getting arguments or print help document.
while getopts 't:f:g:n:h' OPT;do
    case $OPT in
	t) threads="$OPTARG";;
	f) ref="$OPTARG";;
	g) gtf="$OPTARG";;
	n) name="$OPTARG";;
	h) func;;
	?) func;;
	esac
done

#If one of the arguments is missing, then print help document or the command of RNA-seq will going.
if [[ ! -n "$ref" || ! -n "$threads" || ! -n "$gtf" || ! -n "$name" ]];then
    func
    exit 1
else

#Test of the function (deleted)
#    echo "ref is" $ref
#    echo "threads is" $threads
#    echo "gtf is" $gtf
#    echo "name is" $name

     ls *$name.1* >1
     ls *$name.2* >2
     paste 1 2 > 3
     cat 3 | while read id
     do
	    arr=($id)
	    R1=${arr[0]}
	    R2=${arr[1]}
	
        hisat2 -p $threads  -x $ref -1 $R1 -2 $R2 -S ${R1}.sam

        samtools view -bS ${R1}.sam > ${R1}.bam

        samtools sort -l 9 -m 2G -o ${R1}.sort.bam -@ $threads ${R1}.bam

        samtools index -@ $threads ${R1}.sort.bam

done
#Since htseq-count not supports muti-threads, we use GNU-parallel to solve this issue.
ls *.sort.bam | parallel -j $threads htseq-count -f bam -i gene_id -m union {} $gtf '>' {}.txt
Rscript exon_length.R $gtf
Rscript merge_count_to_fpkm.R
fi
