# SuperTranscripts

To construct putative non-redundant assemblies for genotyping, we collapsed Trinity assemblies into SuperTranscripts using the python script provided with Trinity version 2.6.5, Trinity_gene_splice_modeler.py. We then performed genotyping using run_variant_calling_gatk4.py, a wrapper script provided with Trinity that implements the Broad Institute best practices pipeline for genoytping from RNA-seq data, that aligns reads in a splice-aware fashion using STAR, then calls genotypes using GATK version 4. Other than setting the sjdb overhang parameter for STAR to reflect the RNA-seq read lengths, and increasing allocated memory using the --maxram argument, default settings were used.

 
