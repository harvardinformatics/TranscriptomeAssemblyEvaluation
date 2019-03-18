import sys
from collections import defaultdict
covdepthbed=open(sys.argv[1],'r')

matchdict=defaultdict(int)
depth_per_ts_dict=defaultdict(int)
for line in covdepthbed:
    linelist=line.strip().split()
    matchdict[linelist[0]]+=1
    if int(linelist[4])>1:
        depth_per_ts_dict[linelist[0]]+=1
        
        
        
fout=open(sys.argv[1][:-4]+'depth_by_contig.csv','w')
fout.write('EnsTs\tBasesMatch\tContigOverlaps\n')
for key in matchdict.keys():
    fout.write('%s,%s,%s\n' % (key,matchdict[key],depth_per_ts_dict[key]))
    
fout.close()
