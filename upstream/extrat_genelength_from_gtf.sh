#!/usr/bin/bash
#Calculate gene length from gtf file. Depends on perl awk and bedtools. Usage: extrat_genelength_from_gtf.sh exanple.gtf
gtf=$1

awk -v FS="\t" -v OFS=" " '{if($3=="exon")print $4,$5,$9}' $gtf | \
awk -v FS=" " -v OFS="\t" '{gene=substr($4,2,15);print gene,$1,$2}'>$gtf.bed

bedtools sort -i $gtf.bed > $gtf.sort.bed

bedtools merge -i $gtf.sort.bed > $gtf.sort.merge.bed

perl extract_chr_cover_length.from.bed.pl $gtf.sort.merge.bed > genelength.txt

rm $gtf.bed $gtf.sort.bed $gtf.sort.merge.bed
