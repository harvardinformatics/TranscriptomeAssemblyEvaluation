import argparse
from os.path import basename
from sets import Set

def biallelic_snv_filter(gtdict,superfilter=False):
    filter = False
    mapref_allele_lengths = Set()
    for allele in gtdict['maprefalleles'].split(';'):
        if allele != 'NA':
            mapref_allele_lengths.add(len(allele))
    if len(mapref_allele_lengths) > 1:
        filter = True    
    if len(gtdict['maprefalleles'].split(';')) > 2:
        filter = True
    
    if superfilter == True:
        superts_allele_lengths = Set()
        for allele in gtdict['supertsalleles'].split(';'):
            if allele != 'NA':
                superts_allele_lengths.add(len(allele))
        if len(superts_allele_lengths) > 1:
            filter = True
        if len(gtdict['supertsalleles'].split(';')) > 2:
            filter = True
    
    return filter      


if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='filters gtype intersection file and creates an SNV, bi-alleleic only table')
    parser.add_argument('-gt','--gtable',dest='genotypes',type=str,help='table of mapref and superts genotype intersections')
    parser.add_argument('-sf',action='store_true',help="flag to filter superts so biallelic SNVs as well")
    opts = parser.parse_args()

    fopen=open(opts.genotypes,'r')
    if opts.sf == True:
        fout=open('mapref_superts_biallelicfiltered_%s' % basename(opts.genotypes),'w')
    else:
        fout=open('mapref_biallelicfiltered_%s' % basename(opts.genotypes),'w')
    fields = fopen.readline().strip().split('\t')
    fout.write('%s\n' % '\t'.join(fields))
    for line in fopen:
        linelist = line.strip().split('\t')
        linedict = dict(zip(fields,linelist))
        if opts.sf == True:
            filter_eval = biallelic_snv_filter(linedict,superfilter=True)
        else:
            filter_eval = biallelic_snv_filter(linedict,superfilter=False)
        
        if filter_eval == False:
            fout.write(line)

    fout.close()
