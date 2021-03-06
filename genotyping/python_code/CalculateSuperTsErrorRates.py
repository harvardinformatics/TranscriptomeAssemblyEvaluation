#snpid   genomicpositions        supertspositions        mapref_ref      maprefalleles   superts_ref     supertsalleles  maprefcov       mapref_ad       superts_ad      mean_supercov   median_supercov SToverlaps

import sys
from sets import Set
from collections import defaultdict
from Bio import SeqIO
from numpy import median

fopen=open(sys.argv[1],'r')
sid = sys.argv[2]
superts_fasta =  open(sys.argv[3],'r')

error_dict = defaultdict(int)
st_length_dict = {}

for record in SeqIO.parse(superts_fasta,'fasta'):
    st_length_dict[record.id] = len(record.seq)

fields = fopen.readline().strip().split('\t')
fout = open('SuperTsErrors_%s.tsv' % sid,'w')
fout.write('supertsid\tnum_gt_errors\tperr\n')
for line in fopen:
    linelist = line.strip().split('\t')
    linedict =  dict(zip(fields,linelist))
    #print linedict
    ref_alleles = Set(linedict['maprefalleles'].split(';'))
    st_alleles = Set(linedict['supertsalleles'].split(';'))

    if ';' not in linedict['supertspositions'] and linedict['supertspositions'] != 'NA' and 'NA' not in linedict['supertsalleles']:
        if ref_alleles != st_alleles:
            st_error_pos =  int(linedict['supertspositions'].split(':')[1])
            error_dict[linedict['supertspositions'].split(':')[0]]+=1
for st in error_dict:
    fout.write('%s\t%s\t%s\n' % (st,error_dict[st],error_dict[st]/float(st_length_dict[st])))


fout.close()    
    


