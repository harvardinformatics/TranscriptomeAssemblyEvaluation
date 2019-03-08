# Functional Composition of Transcriptome Contigs and Sequence Reads

### Hierarchical classification of annotated nucleotides
We created bed files of reference genome intervals from annotations in gtf format, following a hierarchical classification of unique nucleotides described in the Supplementary Text of Freedman et al. 2019 (bioRxiv). Specific details, with respect to how nucleotides were classifed in Mus, are provided in MethodsToProduceAnnotationBedfiles_Mus.txt, including bedtools execution commands, and the use of FilterBedWithElementList.py in building annotation bed files for particular functional element classes. The resulting annotation bed files are provided in mouse_annotations.tar.gz We modified these methods slightly in the event that a gtf file did not contain the same annotation nomenclature, e.g. the Drosophila annotation does not contain the full complement of biotypes present in the Mus annotation.

### Functional classification of assembly contigs
We mapped assembly contigs to the genome with GMAP, with an example command line as follows:

    gmap -n 0 -b -B 5 -t 15  -d no_patches_Mus_musculus.GRCm38.dna_sm.toplevel -D /PATH/TO/gmap_index/ --failed-input failed.fasta -f samse transcriptome.fasta > GmapToGenome.sam

We then extract the successfully aligned contigs with samtools to a bam file, then convert this bam to bed with bedtools using the -splitD flag that splits on both Ns and Ds in the CIGAR string. We then execute a shell script, BedToContigCompositionBaseCounts.sh, which takes as its one command line argument the bed file of contig alignments to the genome. This script calls both GetTotalBasesMappedPerContig.py and TallyBedIntersect.py. Once the "contig tallies" files have been created for each functional class of nucleotides, we aggregate them into a summary table with WriteFunctCompPerRead.py. 
 
### Functional classification of sequence reads
Reads are aligned to the genome with HiSat using the following generic command line:

    hisat2 -p 8 -x /PATH/TO/HISAT_INDEX/ -q --phred33 --min-intronlen 20 --max-intronlen 500000 -1 R1.fq -2 R2.fq -S alignment.sam  

For stranded libraries, we also use the --rna-strandness switch, and the FR or RF argument, if the library is ligation-stranded or dUTP, respectively.

Quantification of the fraction of reads assigned to different functional classes of nucleotides are conducted in the same was as for contigs, using the same scripts.

 
 
