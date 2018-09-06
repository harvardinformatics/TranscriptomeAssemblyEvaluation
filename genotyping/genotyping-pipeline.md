# SuperTranscripts

## Construction of non-redundant transcriptome reference
To construct putative non-redundant assemblies for genotyping, we collapsed Trinity assemblies into SuperTranscripts using the python script provided with Trinity version 2.6.5, Trinity_gene_splice_modeler.py. This approach collapses sets of contigs belonging to the same Trinity component (i.e. belonging to the same Trinity "gene"). 

## Genotyping with GATK
After constructing SuperTranscripts, We inferred genotypes using run_variant_calling_gatk4.py, a wrapper script provided with Trinity that implements the Broad Institute best practices pipeline for genoytping from RNA-seq data, that aligns reads in a splice-aware fashion using STAR, then calls genotypes using GATK version 4. Other than setting the sjdb overhang parameter for STAR to reflect the RNA-seq read lengths, and increasing allocated memory using the --maxram argument, default settings were used.

## Projection of trancriptome into genomic coordinates
In order to evaluate the ability of genotypes derived from transcriptome assemblies to accurately reconstruct patterns of genomic diversity, one must compare these genotypes to those derived from the map-to-reference genome approach. Doing so requires projecting SuperTranscripts into genomic coordinate space. To do this, we mapped *Mus* SuperTranscripts to the reference genome using GMAP, version 2016-06-30. A generic example of our GMAP command line is as follows:

    gmap -n 0 -b -B 5 -t 15  -d no_patches_Mus_musculus.GRCm38.dna_sm.toplevel -D path/to/mus/GRCm38/gmap_index --failed-input Trinity_failed.fasta -f samse trinity_genes.fasta > Trinitysupertscipts.GmapToGenome.sam

 
