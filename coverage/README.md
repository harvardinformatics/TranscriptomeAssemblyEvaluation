# Calculating Assembly Coverage of Reference Genes and Transcripts

For a given de novo transcriptome assembly, to calculate reference transcript and gene-level coverage, we first use [BLAT](https://genome.ucsc.edu/FAQ/FAQblat) to align transcriptome contigs to reference transcripts. In order to relate coverage to expression level of reference transcripts and genes, we then quantify abundance with [RSEM](https://deweylab.github.io/RSEM/), based upon alignments to transcripts using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

A generic example of our BLAT command line is as follows:

    blat /PATH/TO/TRANSCRIPTS/FASTA  assembly.fasta -t=dna -q=dna -fastMap assembly.ENSEBML_RSEMtscripts.blat.psl

A generic example of our RSEM command line for mapping to reference transcripts is as follows:

    rsem-calculate-expression --bowtie2 -p 16  --output-genome-bam --time --paired-end R2.fq R1.fq /PATH/TO/RSEM/INDEX output_prefix

For dUTP and ligation-stranded libraries, we added --forward-prob 0 and --forward-prob 1, respectively to the command line.  
We then calculate transcript and gene-level coverage by the de novo assembly with OmnibusAssemblyCoverageAnalysis.py, supplying it via command line arguments the BLAT output, the map-to-reference gene and isoform expression matrices, and a tab-separated file where the first column contains a reference transcript id, and the second column cotains the reference gene it is associated with. The expression matrices are in transcripts-per-million (TPM) units, and have a header. The transcript-gene map file does not have a header.  

This script outputs coverage in terms of:  

  * coverage of individual transcripts by its best-hit contig
  * summed coverage of transcripts by all contigs that are hits, and
  *  summed coverage of genes by all contigs that are hits to that gene



**It is important to note** that this code was written to handle assemblies build from a single sample. To adapt this pipeline to an assembly derived from multiple samples, one could either generate a new expression matrix file that sums TPMs from individual-sample outputs, or modify the python code to handle a multi-sample fastq file and generate a single TPM value for each transcript or gene.  

We examine the combined effects of:

  * TPM filtering on the transcriptome assembly by filtering out contigs with TPM < 1,
  * Only considering protein-coding reference annotations (of primary interest for users of transcriptome assembly), and
  * Only considering reference transcripts whose minimum length is 200bp, because this is the minimum default length for a TRINITY contig.

To do this, we must filter the BLAT output to address the TPM filter, and use FilterTableEntries.py to do so. Because OmnibusAssemblyCoverageAnalysis.py iterates over the reference transcript and gene ids collected from the expression matrix, we also must filter the expression matrices so we only consider protein-coding annotations and transcripts with a miminum length of 200 bp. We do this using FilterMapRefExprMatrices.py. This script takes as in files the map-to-reference expression matrices, lists of protein-coding gene and transcript ids, the transcript-to-gene map described above, and a list of reference transcripts whose lengths are greater or equal to 200bp. We obtained this last file using the [BioPython](https://biopython.org/) SeqIO module for fasta processing, which can return sequence lengths and ids for easy parsing. 

    
    
    
