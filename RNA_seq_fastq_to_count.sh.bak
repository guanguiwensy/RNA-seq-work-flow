#¡/bin/bash
t=44
ref=/home/guanguiwen/data1/reference/human/hg19/hisat_hg19_HBV/ht2/ucsc.hg19.HBV.fasta
gtf=/home/guanguiwen/data1/reference/human/hg19/gtf/hg19.ensGene.HBV.gtf

#$1 is the id¡¡contain file
cat $1 |while read id
do
	arr=($id)
	R1=${arr[0]}
	R2=${arr[1]}
hisat2 -p $t -x $ref -1 $R1 -2 $R2 -S $R1.sam

samtools view -bS -1 $R1.sam > $R1.sam.bam

samtools sort -l 9 -m 1G -o $R1.sam.bam.sort.bam -T sorted -@ $t $R1.sam.bam

samtools index $R1.sam.bam.sort.bam -@ $t

rm $R1.sam $R1.sam.bam 

featureCounts -T $t -p -t exon -g gene_id  -a $gtf -o $R1.count.txt $R1.sam.bam.sort.bam; done

done
mkdir count
mv *.count.txt count
cp ~/biotools/RNA-seq/*.csv count
cd count
Rscript ~/biotoosl/RNA-seq/RNA_seq_merge.R

