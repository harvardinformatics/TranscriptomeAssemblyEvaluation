# SuperTranscripts

## Construction of non-redundant transcriptome reference
To construct putative non-redundant assemblies for genotyping, we collapsed Trinity assemblies into SuperTranscripts using the python script provided with Trinity version 2.6.5, Trinity_gene_splice_modeler.py. This approach collapses sets of contigs belonging to the same Trinity component (i.e. belonging to the same Trinity "gene"). 

## Genotyping with GATK
After constructing SuperTranscripts, We inferred genotypes using custom modifications run_variant_calling.py, a wrapper provided with Trinity. This script implements the Broad Institute best practices pipeline for genoytping from RNA-seq data, that aligns reads in a splice-aware fashion using STAR, then calls genotypes using GATK. For all analyses, we modified the script to point to a more current version of GATK (version 4). Because some of our *Mus* samples are effectively polyploid -- BALB/c and wild *Mus* samples are pools of 6 and 8 individuals, with ploidies of 12 and 16, respectively -- for these samples we added a hard-coded value for the "--sample-ploidy" optional argument. We also supplied a value for the sjdb overhang parameter for STAR to reflect the RNA-seq read lengths, and increased allocated memory using the --maxram argument when necessary. After the initial round of genotyping, genotypes were filtered using defaults, and in the exact same manner as with genotyping from the reference genome: FS > 30.0, and QD < 2.0.

## Coverage depth calculations
The default mode of run_variant_calling.py is to only emit genotype calls at variant sites. Information on the total number of sites with sufficient depth of coverage for genotyping is necessary both for diagnosing genotyping errors, and for calculating genetic diversity measures such as heterozygosity. To determine which sites had sufficient coverage depth for genotyping we converted the gtf file output by Trinity_gene_splice_modeler.py to bed format, collapsed this bed file into supertranscript intervals with the bedtools mergeBed utility, and calculated depth of coverage of STAR-aligned RNA-seq reads on the supertranscripts. For this last step, we converted the STAR alignments split on N cigar string flags (resulting from splicing) to bed files using the bedtools bamToBed utility, then calculated depth of coverage for each SuperTranscript nucleotide position using bedtools coverageBed with the "-d" flag. The resulting coverage file is not in bed format, so we used basic awk operations to convert this to proper format for downstream analyses. We created a filtered coverage bed of callable sites, defined as those positions with coverage depth \ge 5. 

## Projection of trancriptome into genomic coordinates
In order to evaluate the ability of genotypes derived from transcriptome assemblies to accurately reconstruct patterns of genomic diversity, one must compare these genotypes to those derived from the map-to-reference genome approach. Doing so requires projecting SuperTranscripts into genomic coordinate space. To do this, we mapped *Mus* SuperTranscripts to the reference genome using GMAP, version 2016-06-30. A generic example of our GMAP command line is as follows:

    gmap -n 0 -b -B 5 -t 15  -d no_patches_Mus_musculus.GRCm38.dna_sm.toplevel -D path/to/mus/GRCm38/gmap_index --failed-input Trinity_failed.fasta -f samse trinity_genes.fasta > Trinitysupertscipts.GmapToGenome.sam

where trinity_genes.fasta is the SuperTranscript sequences. We filtered out SuperTranscripts that did not map to the *Mus* genome using samtools.

 
