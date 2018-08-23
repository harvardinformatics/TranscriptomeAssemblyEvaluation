# Trinity

We build Trinity assemblies with version 2.2.1, with strandedness flags (--RF and --FR for dUTP and ligation-stranded protocols, respectively) where appropriate. A generic example command line for a dUTP library is is as follows.

    /path/to/Trinity --seqType fq --max_memory 320G  --min_kmer_cov 1 --SS_lib_type RF --grid_node_max_memory 5G --left trimmed_corrected_1.fq  --right trimmed_corrected_2.fq --output trinity-out --CPU 32 --grid_conf grid.conf

We parallelized operations, where possible, using the --grid_conf argument, with compute farm settings for SLURM set in the grid.conf file available in this directory.
