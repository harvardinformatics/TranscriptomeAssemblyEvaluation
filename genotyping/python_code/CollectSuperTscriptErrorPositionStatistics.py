import sys
from Bio import SeqIO
from numpy import array,median
from numpy.random import uniform
from scipy.stats import ranksums
from sets import Set

def NormalizeDistance(a,b,xmin,xmax,val):
    normed = a + ((val-xmin)*(b-a))/float(xmax-xmin)
    return normed

def SimulateAndCompareRel2Mids(observed):
    print len(observed)
    simdata = uniform(-1,1,len(observed))
    abssim = abs(simdata)
    absobs=abs(array(observed))
    abs_statistic,abs_p = ranksums(abssim,absobs)
    statistic,p=ranksums(simdata,observed)
    return {'abs_norm':{'median_observed' : median(absobs), 'median_simulated' : median(abssim), 'wilcoxon' : abs_statistic, 'p': abs_p},'norm' : {'median_observed' : median(observed),'median_simulated' : median(simdata),'wilcoxon' : statistic, 'p' : p}} 
    

fopen=open(sys.argv[1],'r')
sid = sys.argv[2]
superts_fasta =  open(sys.argv[3],'r')
st_length_dict = {}

for record in SeqIO.parse(superts_fasta,'fasta'):
    st_length_dict[record.id] = len(record.seq)

rel2mids = []

for record in SeqIO.parse(superts_fasta,'fasta'):
    st_length_dict[record.id] = len(record.seq)


fields = fopen.readline().strip().split('\t')
statsout = open('position_ranksumtest_%s.tsv' % sid,'w')


fout = open('SuperTsErrors_%s.tsv' % sid,'w')
fout.write('supertsid\td2end\treldist2end\tnormd2mid\n')
for line in fopen:
    linelist = line.strip().split('\t')
    linedict =  dict(zip(fields,linelist))
    ref_alleles = Set(linedict['maprefalleles'].split(';'))
    st_alleles = Set(linedict['supertsalleles'].split(';'))

    if ';' not in linedict['supertspositions'] and 'NA' not in linedict['supertsalleles']:
        if ref_alleles != st_alleles:
            st_length = st_length_dict[linedict['supertspositions'].split(':')[0]]
            st_error_pos =  int(linedict['supertspositions'].split(':')[1])
            dist2end = st_length - st_error_pos
            reldist2end = (st_length - st_error_pos)/float(st_length)
            rel2mid = NormalizeDistance(-1,1,1,st_length,st_error_pos)
            rel2mids.append(rel2mid)

            fout.write('%s\t%s\t%s\t%s\n' % (linedict['supertspositions'].split(':')[0],dist2end,reldist2end,rel2mid))
fout.close()    


test_statistics = SimulateAndCompareRel2Mids(rel2mids)
print test_statistics    

statsout.write('### absolute normalized distances: 0.0 to 1.0 ###\n')
statsout.write('Median observed\t%s\n' % test_statistics['abs_norm']['median_observed'])
statsout.write('Median simulated\t%s\n' % test_statistics['abs_norm']['median_simulated'])
statsout.write('Wilcoxon statistic\t%s\n' % test_statistics['abs_norm']['wilcoxon'])
statsout.write('P-value\t%s\n' % test_statistics['abs_norm']['p'])
statsout.write('\n\n### normalized distances: -1 to 1 ###\n')
statsout.write('Median observed\t%s\n' % test_statistics['norm']['median_observed'])
statsout.write('Median simulated\t%s\n' % test_statistics['norm']['median_simulated'])
statsout.write('Wilcoxon statistic\t%s\n' % test_statistics['norm']['wilcoxon'])
statsout.write('P-value\t%s\n' % test_statistics['norm']['p'])


statsout.close()
