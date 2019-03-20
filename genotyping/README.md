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

## Projection of trancriptome into genomic coordinates
In order to evaluate the ability of genotypes derived from transcriptome assemblies to accurately reconstruct patterns of genomic diversity, one must compare these genotypes to those derived from the map-to-reference genome approach. Doing so requires projecting SuperTranscripts into genomic coordinate space. To do this, we mapped *Mus* SuperTranscripts to the reference genome using GMAP, version 2018-07-04. Note: this version represents an update for that used to assess functional composition, that includes updates to handling CIGAR strings,so as to avoid any possible errors in the conversion of ST genotypes to genomic coordinates.  

A generic example of our GMAP command line is as follows:

    gmap -n 0 -b -B 5 -t 15  -d no_patches_Mus_musculus.GRCm38.dna_sm.toplevel -D path/to/mus/GRCm38/gmap_index --failed-input Trinity_failed.fasta -f samse trinity_genes.fasta > Trinitysupertscipts.GmapToGenome.sam

We then use samtools to extract only the aligned SuperTranscripts and convert the output to bam:


    samtools view -b -F 4 Trinitysupertscipts.GmapToGenome.sam > alignedonly_Trinitysupertscipts.GmapToGenome.bam

Next, we convert this bamfile to bed, retaining the full cigar string for each alignment, using bedtools bamToBed utility and supplying the -cigar argument. Conversion of genotypes to genomic coordinate space relies on parsing the ST alignment CIGAR strings. This, filtering on exons, annotation with coverage depth, and other related operations are carried out with [ConvertSuperTranscriptDataToGenomicCoordinates.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/ConvertSuperTranscriptDataToGenomicCoordinates.py), which takes as input:  
    
* the bed file of ST genomic alignments with CIGAR strings, via -i
* the ST fasta, via -f
* the sorted, merged bed file of genomic exon intervals, via -ex
* the output of coverageBed on the STs, i.e. prior to conversion to proper bed format, via -sd
* and the filtered ST vcf-format genotypes file, via -v 

## Intersection of MR and ST genotypes

Prior to intersecting the genotypes, we first convert the MR genotypes to bed format with [ConvertMaprefVcfToBed.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/ConvertMaprefVcfToBed.py), and only retain genotypes overlapping exons. This script takes as its arguments:  

* the MR filtered vcf file, via -rvcf
* the merged bedfile of annotated exons, via -e

Next, we perform the intersection with [IntersectMapRefandSuperTranscriptGenotypes.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/IntersectMapRefandSuperTranscriptGenotypes.py), which takes as its arguments:  

* the bed file of filtered, exonic MR genotypes, via -m
* a bed file representing a column reorganization and truncation of the exon-filtered ST genotypes bed that begins with "wGenomePosAndSuperTsDepth", and that was generated by ConvertSuperTranscriptDataToGenomicCoordinates.py, via -s
  * the reorganized file is produced by *awk '{print $15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' wGenomePosAndSuperTsDepthFile.bed > OUTFILE*, such that it is back in genomic coordinate space.
* the raw bedtools output from coverageBed for MR exon-filtered genomic intervals, via -mrc
* the output from coverageBed for ST,in proper bed format, i.e. the coverage bed with "reformatted" as prefix generated by ConvertSuperTranscriptDataToGenomicCoordinates.py, via -scov
* the 1-bp resolution bedfile of exonic sites in ST coordinate space, with genomic coordinates as values in additonal columns) produced previously by ConvertSuperTranscriptDataToGenomicCoordinates.py, which has the prefix "supertscoords", via -sexons
* the desired name of the tab-separated output file for the intersection, via -o
* the genome fasta file, via -gf

As part of the intersection, and because we were focused on which alleles were recovered, and whether they matched between methods, this script collapse polyploid genotyes (e.g. with GATK ploidy argument set to 12, we convert 0/0/0/0/0/0/1/1/1/1/1/1 to 0/1) such that every recorded allele is only recorded once.

## Intersection table filtering and modification

We perform an additional filtering step such that we only consider intersections where MR genotypes are bi-allelic. We do this with [FilterGtypeIntersectionsForSnvBiallelic.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/FilterGtypeIntersectionsForSnvBiallelic.py). This script takes as arguments:  

* the intersection file, via -gt
* an option flag to indicate whether intersections will also filter on bi-allelic sites in STs.

Because we were interested in the ability of ST to recover "true" bi-allelic genotypes as called by MR, we did not use the -sf flag in our analyses of genotype concordance between the two methods.

In our original development of IntersectMapRefandSuperTranscriptGenotypes.py, we only reported ST positions if SNPs were called in them. As a result, false negative sites would have no position recorded. Because we were interested in whether coverage might influence false negatives, we modified the bi-allelic filtered intersection file, appending, for each ST that mapped to a genomic position, the contig id, the coverage depth, and the strand of the mapping as an additional table field. We do that with [UpdateContigOverlapsSnpIntersectTable.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/UpdateContigOverlapsSnpIntersectTable.py), which takes as ordered command line arguments (via sys.argv specification):
  
* 1st argument, the intersection file
* 2nd argument the exon filtered ST coverage bed for the sample, that ends in "exons_coverage.bed", produced by ConvertSuperTranscriptDataToGenomicCoordinates.py

In examing outputs, we realized that there were certain cases where the genotypes from the two methods were effective concordant were being called as errors, specifically:  

* cases where there was no ST genotype call, and the ST base corresponds to a fixed alternative variant in MR, are incorrectly called as FN
* cases where a SNP is called in ST, and it is a fixed alternative (relative to the ST reference) that matches the reference genome base, are incorrectly called false positives. This can result from a misassembly in the ST.

In a less than optimal way of recifying this, when we calculate error/condcordance metrics based upon the genotype intersection file updated with UpdateContigOverlapsSnpIntersectTable.py, we correct these two errors on the fly, calculating metrics on the corrected row entries of the table, and writing a new table that begins with the prefix "corrected", for use in other downstream analyses. The script that does this, [CalculateConcordanceMetrics_modular.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/CalculateConcordanceMetrics_modular.py), takes as arguments:  

* the genotype intersection table, via -gt
* the sample id, for the purpose of adding it to the concordance table name, via -sid

In this script, each metric is calculated with its own function. While this leads to some inefficiency in processing from a coding perspective, our hope is that it will increase clarity with respect to how each metric is calculated.

Because we were interested in the effect on concordance metrics of annotated vs. unannotated contigs, for a subset of our samples we ran concordance metric analyses in which, if the Trinity ST had no contig with a hit to Uniref90, we set all values used to evaluate concordance to "NA", i.e. to missing, and filtered out entirely any entries where the original entry had no genotype called in MR. To do this, we use, [NoProteinHitSuperTsGenotypesToNA.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/NoProteinHitSuperTsGenotypesToNA.py), which takes as arguments:  

* the genotype intersection table, via -gt
* a table of BLASTP or BLASTX hits, via -bhits

To summarize the locations of genotyping errors in STs, and to evaluate whether such errors are enriched at ends of STs, we analyse the intersection tables with [GetSuperTsErrorCount.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/GetSuperTsErrorCount.py), which takes as ordered command line arguments (via sys.argv):  

* 1st argument, the intersection table
* 2nd argument, a sample id
* 3rd argument, the ST fasta.

## Heterozygosity calculations

We provide three python scripts to obtain counts of heterozygous genotype calls. These counts are effectively the numerator, and the total number of callable sites (with depth greater than or equal to 5) are the denominator in making heterozygosity-per-base estimates presented in our manuscript.

[GetHetcountFromGenomeCoordSupertsGtypeBed.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/GetHetcountFromGenomeCoordSupertsGtypeBed.py) obtains het counts from the "wGenomePosAndSuperTsDepth" file reorganized back into genomic coordinate space (via the rather painful awk command provided above), taking this file as a single command line argument.  

[ExtractHetCountFromMaprefGtypeBedfile.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/ExtractHetCountFromMaprefGtypeBedfile.py), which takes as arguments:  
* an exon-filtered MR genotype bed file, via -b
* a flag '-p' that, if called, only includes genotypes with the GATK "PASS" filter flag, e.g excluding "SnpCluster" entries.

To determine whether focusing on exons might change the relative difference in estimates of heterozygosity for MR and ST, we also generate counts of heterozygous sites prior to exon filtering, by extracting these directly from the filtered vcf files. We do this with [ExtractHetCountsFromVcfs.py](https://github.com/harvardinformatics/TranscriptomeAssemblyEvaluation/blob/master/genotyping/python_code/ExtractHetCountsFromVcfs.py), which takes as arguments:  

* the input vcf file, via -vvcf
* a flat text outfile containing the number of counted heterozygous sites, via -o
* the name of an outfile incluing in bed format, the recorded heterzygous sites, via -b
* the ploidy of the sample, used in collapsing a polyploidy genome into the uniquely observed alleles via -p

In both exon-filtered and unfiltered estimates of heterozygosity, we only consider sites with the GATK "PASS" filter, and the scripts are written to this effect.  
