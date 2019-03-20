# TranscriptomeAssemblyEvaluation
Workflows and scripts used to perform analyses for **Error, noise and bias in de novo transcriptome assemblies**, by Adam H. Freedman, M. Clamp, and T.B. Sackton posted on bioRxiv on 20 March 2019. 

## Python version and required modules
Python scripts are written to work with python 2.7. Several scripts include calls to Biopython, specifically by calling "from Bio import SeqIO", as well as numpy methods. Therefore, replication of these analyses will require Biopython, and numpy, with the latter easily obtainable via an Anaconda python distribution.  

In some cases, external packages such as bedtools and tabix are called. Please review scripts before execution to determine what must be in the $PATH variable for them to work properly. The scripts in the genotyping directory are organized as a module, within which some cross-importing is done. For all of the genotyping-related scripts to work properly, Make sure the directory has been added to the $PYTHONPATH variable.
