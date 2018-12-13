from sets import Set
import sys
from subprocess import Popen,PIPE
import os
mapin = open('/n/holylfs/LABS/informatics/adamf/refs/mus/GRCm38/gene_intervals/genes.txt','r')
gtf = sys.argv[1] 
genes = Set()
for line in mapin:
    genes.add(line.strip())

sample_bed = sys.argv[2]
fout= open('%s_test3_gene_uniquenucleotide_counts.txt' % sample_bed[:-3],'w')
fout.write('gene\tcount\n')

for gene in genes:
    #grepcmd ="grep %s %s |awk '$3!=\"transcript\"{print $0}'|awk '$3!=\"gene\"{print $0}' |awk '{print $1\"\\t\"$4\"\\t\"$5}' > %s.bed" % (gene,gtf,gene)
    grepcmd ="grep %s %s |awk '$3==\"exon\"{print $0}' |awk '{print $1\"\\t\"$4\"\\t\"$5}' > %s.bed" % (gene,gtf,gene)
    #print grepcmd
    proc = Popen(grepcmd,shell=True,stdout=PIPE,stderr=PIPE)
    grepout,greperr =  proc.communicate()
    #print 'greperr',greperr
    mergecmd = 'mergeBed -i %s.bed > merged.bed' % (gene) 
    proc2 = Popen(mergecmd,shell=True,stdout=PIPE,stderr=PIPE)
    mergeout,mergerr =  proc2.communicate()
    #print 'mergeerr',mergerr
    intersectcmd = 'intersectBed -a merged.bed -b %s > intersect.bed' % sample_bed
    proc3 = Popen(intersectcmd,shell=True,stdout=PIPE,stderr=PIPE)
    intersectout,intersecterr =  proc3.communicate()
    #print 'intererr',intersecterr
    bases = 0
    fopen = open('intersect.bed','r')
    for line in fopen:
        linelist = line.strip().split()
        interval = int(linelist[2]) - int(linelist[1])
        bases += interval
    fout.write('%s\t%s\n' % (gene,bases))
    os.system('rm intersect.bed merged.bed %s.bed' % gene)
fout.close()
