# Intra-assembly redundancy

To assess the redundancy within assemblies across all assembled contigs, we perform a BLAT all-by-all alignment. We used a similar command line to that used to align contigs to reference transcripts for the analysis of transcript coverage, i.e.

    blat assembly.fasta assembly.fasta -t=dna -q=dna -fastMap allbyall.psl

with the only difference that we partitioned the assembly fasta into pieces for alignment to the whole transcriptome so that we could parallelize alignment across a job array.

Because in an all-by-all analysis, every pair of aligned contigs will show up twice (one will be a query the other a target, and vice versa), we removed duplicate entires with an rscript, 
