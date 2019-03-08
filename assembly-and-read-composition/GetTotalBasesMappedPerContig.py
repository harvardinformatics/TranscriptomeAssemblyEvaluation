import sys 
from collections import defaultdict
bedin=open(sys.argv[1],'r')
bedout=open(sys.argv[2],'w')
suffix=sys.argv[3]
beddict=defaultdict(int)
for line in bedin:
    linelist=line.strip().split()
    beddict[linelist[3]]+=int(linelist[2])-int(linelist[1])
bedout.write('contig\t%s\n'% suffix)    
for key in beddict.keys():
    bedout.write('%s\t%s\n' % (key,beddict[key]))

bedout.close()
