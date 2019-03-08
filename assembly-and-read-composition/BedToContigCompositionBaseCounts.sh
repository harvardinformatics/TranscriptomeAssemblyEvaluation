#!/bin/bash
# $1 is a bed file of GMAP contig alignments to the genome, created by bamToBed -splitD -i <bamfile>.bam > <bamfile>.be
d
intersectBed -wo -a $1 -b merge_sorted_included_CDS_Mus_musculus.GRCm38.83.sortedGtfTo.bed > cds.bed

subtractBed -a $1 -b Mus_musculus.GRCm38.83.sortedGtfTo.bed> intergenic_external.bed

intersectBed -wo -a $1 -b merge_sorted_intergenic_within_ts_boundaries.bed > intra_intergenic.bed

intersectBed -wo -a $1 -b  merge_sorted_putative_intronic.bed > intron.bed

intersectBed -wo -a $1 -b merge_sorted_noncoding_tsrm_rmCdsUtrStopStart.bed > noncoding.bed

intersectBed -wo -a $1 -b merge_sorted_filt_noncodingintrons.bed > noncoding_introns.bed

intersectBed -wo -a $1 -b merge_sorted_startAndstopCodonsProteinCoding_CDSsubtract.bed > startstop.bed

intersectBed -wo -a $1 -b merge_sorted_utr_proteincoding_Mus_musculus.GRCm38.83.sortedGtfTo.bed > utr.bed 

python GetTotalBasesMappedPerContig.py intergenic_external.bed intergenic_external_contigtallies.txt inter_ext 

python GetTotalBasesMappedPerContig.py $1 mappedbasecount_per_contig.txt contigcount

for i in cds.bed intra_intergenic.bed intron.bed noncoding.bed noncoding_introns.bed startstop.bed utr.bed
do
functionclass=`echo $i |sed 's/.bed//g'`

python TallyBedIntersect.py $i mapped_contigs.txt ${functionclass}_contigtallies.txt ${functionclass} 

done
