

import sys
from collections import defaultdict
bedin=open(sys.argv[1],'r')
contigs=open(sys.argv[2],'r')
contig_list=[]
for line in contigs:
    contig_list.append(line.strip())
    
fout=open(sys.argv[3],'w')
prefix=sys.argv[4]

fdict=defaultdict(int)

for line in bedin:
    bedlist=line.strip().split()
    fdict[bedlist[3]]+=int(bedlist[-1])
    
    
fout.write('contig\t%s\n' % prefix)
for contig in contig_list:
    if contig in fdict:
        fout.write('%s\t%s\n' % (contig,fdict[contig]))
    else:
         fout.write('%s\t%s\n' % (contig,0)) 
         
         
fout.close()  
