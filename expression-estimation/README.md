# Expression Estimation

## Map-to-reference

To generate our truth-set benchmark for RNA abundance estimation we used [RSEM](https://deweylab.github.io/RSEM/). Specifically, we used RSEM to build an index using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as an aligner with the following generic command line:  

    rsem-prepare-reference -p 12 --bowtie --bowtie2 --gtf Mus_musculus.GRCm38.83.gtf Mus_musculus.GRCm38.dna_sm.toplevel.fa rsem_ensembl_indices/Mus_musculus.GRCm38.83

Then, for each sample we used RSEM to estimate RNA abundance from the bowtie2 alignments:

    rsem-calculate-expression --bowtie2 -p 16 --output-genome-bam --time --paired-end R1.fq R2.fq /PATH/TO/RSEM/INDEX/DIR/Mus_musculus.GRCm38.83 OUTFILE_PREFIX

## De novo transcriptomes

For each sample-assembler combination, we did two rounds of expression estimation. First, we did a "naive" estimation where we did not attempt to group contigs by their BLAT hits to annotated *Mus* transcripts

### TRINITY assemblies

We used a perl script provided by the developers to wrap RSEM, that builds the bowtie2 index and generates "gene" and "isoform" level abundance estimates for Trinity components and contigs, respectively:  

    align_and_estimate_abundance.pl --prep_reference --thread_count 16 --transcripts Trinity.fasta --seqType fq --left R1.fq --right R2.fq --est_method RSEM --aln_method bowtie2 --trinity_mode --output_dir OUTDIR --output_prefix PREFIX


For dUTP and ligation stranded libraries, we supplied the --SS_lib_type flag values of RF and FR, respectively. 

We then performed a second round, where we grouped TRINITY contigs into genes according to their top-scoring BLAT hits. In doing this, we did not use align_and_estimate_abundance.pl, as we discovered that the argument in which one supplies it the map table DOES NOT WORK!! It simply returns abundance estimates based upon the TRINITY component groupings. Thus,and, per RSEM's requirements, we generated a table that maps contig ids to gene ids, and build the index as follows:

    rsem-prepare-reference --transcript-to-gene-map MAP_TABLE --bowtie2 -p 8 /PATH/TO/Trinity.fasta /PATH/TO/INDEX/Trinity

We then estimated abundance:

    rsem-calculate-expression  --bowtie2 -p 16 --time --paired-end R1.fq R2.fq PATH/TO/INDEX/ OUTFILE_PREFIX

For dUTP and ligation stranded libraries, we supplied the --forward-prob flag values of 0 and 1, respectively.

## SHANNON and BINPACKER assemlbies

Both naive and annotation-grouped rounds of abundance estimation were performed by calling RSEM directly, as described above. Neither method groups contigs into putative gene clusters, so naive estimation produces only contig-level expression estimates (although RSEM does produce isoforms.results and genes.results files, with the latter only referencing one isoform).

## Comparison of Map-to-reference and Map-to-Transcriptome estimates

