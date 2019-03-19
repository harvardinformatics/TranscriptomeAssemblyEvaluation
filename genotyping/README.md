# SNP-based analysis of assembly composition

Our comparision of map-to-reference (MR) and SuperTranscript (ST) genotypes involves four main steps:

    * Inferring genotypes from alignments of RNA-seq reads to the reference genome
    * Collapsing TRINITY assemblies into STs and inferring genotypes from alignments to STs
    * Projecting ST genotypes into genomic coordinates
    * Calculating various performance metrics on ST genotypes using MR as the gold standard

## MR genotyping
We implement the Broad Institute's [best practices](https://gatkforums.broadinstitute.org/gatk/discussion/3892/the-gatk-best-practices-for-variant-calling-on-rnaseq-in-full-detail) for genotyping with RNA-seq data. The only step we do not utilize is base quality recalibration (BQR). This step requires known SNPs, such as those in dbSNP, as a truth set. BQR is not part of the TRINITY developers' best practice for genotyping for the obvious reason that non-model organisms do not have dbSNP resources. Therefore, to make genotypes derived from the two methods more comparable, we do not do BQR.

Genotyping involves the following steps:  

### Build a STAR genome index
A 1st pass genome index for [STAR](https://github.com/alexdobin/STAR) is built. A generic example command line is as follows:  

    STAR --runMode genomeGenerate --runThreadN 8 --genomeDir PATH/TO/OUTPUT/INDEX/ --genomeFastaFiles /PATH/TO/Mus_musculus.GRCm38.dna_sm.toplevel.fa   --sjdbGTFfile /PATH/TO/Mus_musculus.GRCm38.83.gtf   --sjdbOverhang 75 

The --sjdbOverhang argument is set to one minus the read length. In our case, it was 75 for the MDC sample, and 99 for all other *Mus* libraries. 

### First-pass read alignment

We do this as follows:

    STAR --genomeDir PATH/TO/GENOME/INDEX --readFilesIn PATH/TO/R1.fa /PATH/TO/R2.fq  --runThreadN 12 
 
### Re-build genome index
This step uses the initial read-mapping to update information on splice junctions. We create a new directory for the updated index, create a symlink of the genome fasta to that new directory, and from the directory where read mapping was done (that contains the splice-junction file SJ.out.tab), do:

    STAR --runMode genomeGenerate --genomeDir  /PATH/TO/NEW/INDEX --genomeFastaFiles /PATH/TO/NEW/INDEX/Mus_musculus.GRCm38.dna_sm.toplevel.fa --sjdbFileChrStartEnd `pwd`/SJ.out.tab --sjdbOverhang 75 --runThreadN 12

### 2nd-pass read alignment

    STAR --genomeDir /PATH/TO/NEW/INDEX/ --readFilesIn PATH/TO/R1.fa /PATH/TO/R2.fq --runThreadN 16

### Add read groups
You will need an updated jdk development version of java for this to work, and [picard](https://broadinstitute.github.io/picard/). The command line below uses our specified settings for the MDC library, which is a dUTP, stranded library, therere we set RGLB to RF. Other specified arguments are RGID and RGSM.

    java -Xmx32g -jar /PATH/TO/picard.jar AddOrReplaceReadGroups I=Aligned.out.sam O=rg_added_sorted_SRR203276.bam SO=coordinate RGID=SRR203276 RGLB=RF RGPL=Illumina RGPU=Unknown RGSM=Mus_MDC TMP_DIR=temp 
    
### Mark duplicates
Similarly, using picard we mark duplicates:

    java -Xmx32g -jar /PATH/TO/picard.jar MarkDuplicates I=rg_added_sorted_SRR203276.bam O=dedupped_rg_added_sorted_SRR203276.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=markdup.metrics TMP_DIR=temp

### Split Reads on CIGAR string Ns
This is the first step that uses GATK, so be sure it is in your path to do the following:

    gatk SplitNCigarReads --java-options "-Xmx32g" -R /PATH/TO/NEW/INDEX/Mus_musculus.GRCm38.dna_sm.toplevel.fa -I dedupped_rg_added_sorted_SRR203276.bam -O split.dedupped_rg_added_sorted_SRR203276.bam

### Variant calling with GATK HaplotypeCaller

    gatk HaplotypeCaller --java-options "-Xmx32g"  -R /PATH/TO/NEW/INDEX/Mus_musculus.GRCm38.dna_sm.toplevel.fa -I split.dedupped_rg_added_sorted_SRR203276.bam --dont-use-soft-clipped-bases -stand-call-conf 20.0 -O MDC_maptoref_gatk4.vcf

For samples that are know to be pools of separate individuals--BALB/c and wild pools analyzed in this study are comprised of 6 and eight individuals, respectively-- we also included the "--sample-ploidy" with our execution of HaplotypeCaller, specifying the diploid numbers of 12 and 16.

### Filtering variant calls
We employ variant call filters recommended by the Broad Institute for genotypes derived from RNA-seq data. Note, identical filters are applied to both the MR and ST-derived genotypes.

    gatk VariantFiltration --java-options "-Xmx16g" -R /PATH/TO/NEW/INDEX/Mus_musculus.GRCm38.dna_sm.toplevel.fa -V MDC_maptoref_gatk4.vcf -window 35 -cluster 3 --genotype-filter-name FS --genotype-filter-expression "FS > 30.0" --genotype-filter-name QD --genotype-filter-expression "QD < 2.0" -O filtered_MDC_maptoref_gatk4.vcf

## Collapsing Trinity Assembly into SuperTranscripts

To construct putative non-redundant assemblies for genotyping, we collapsed Trinity assemblies into SuperTranscripts using the python script provided with Trinity version 2.6.5, Trinity_gene_splice_modeler.py. This approach collapses sets of contigs belonging to the same Trinity component (i.e. belonging to the same Trinity "gene"). As described by [Davidson et al. 2017, Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1284-1), and to directly quote the [Trinity developers](https://github.com/trinityrnaseq/trinityrnaseq/wiki/SuperTranscripts), "A SuperTranscript is constructed by collapsing unique and common sequence regions among splicing isoforms into a single linear sequence."

## ST Genotyping
We inferred genotypes using custom modifications run_variant_calling.py, a wrapper provided with Trinity. This script implements the Broad Institute best practices pipeline that we have documented above. Because we used a more recent version of GATK (version 4) than was initially available since this script was made available, For all analyses, we modified the script to point to a local install of version 4. Parameters were provided to match those used in MR genotyping with respect to sjdb overhang and ploidy. As already noted, the same set of filters were imposed on the variant calls to produce a final, filtered set.

## Coverage Depth Calculations

The default mode of HaplotypeCaller, both as implemented directly for MR and for ST via run_variant_calling.py, is to only emit genotype calls at variant sites. Information on the total number of sites with sufficient depth of coverage for genotyping is necessary both for diagnosing genotyping errors, and for calculating genetic diversity measures such as heterozygosity per base. 

For MR, To determine which sites had sufficient coverage depth for genotyping, we converted the STAR output bam file to bed with [bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html), using the bamToBed utility and the -splitD flag to split on N and D CIGAR string flags. In order to create a target bed file on which to calulate coverage depth, We then sorted and merged this bed file with the bedtools sortBed and mergeBed utilities. Depth of coverage was then calculated with the coverageBed utility using the -d flag to output base-level coverage, and supplying the sorted, merged bed to the -a flag, and the sorted, unmerged bed to -b.The resulting coverage file is not in standard bed format, so we used basic awk operations to convert this to proper format for downstream analyses.  

For ST, we first converted the gtf file output by Trinity_gene_splice_modeler.py to bed format by extracting the relevant columns with simple awk commands, converting SuperTranscript intervals to zero start (i.e. subtracting 1 from the beginning of each interval). Next, we sorted and merged this bed file into unique SuperTranscript intervals, and calculated coverage on these intervals with bedtools in the same manner as we did for MR.  For ST, these operations are carried out from within [ConvertSuperTranscriptDataToGenomicCoordinates.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/ConvertSuperTranscriptDataToGenomicCoordinates.py) (see below).  
For both MR and ST, we produce filtered callable sites bed files by filtering out all sites with depth < 5. In additionl, for analyses in which we only consider sites overlapping reference exons, we further filter the callable sites bed file by using the bedtools intersectBed utility with a bed file of *Mus* exons which we extracted from the original annotation file in gtf format.

