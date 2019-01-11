import argparse
from os.path import basename


def biallelic_snv_filter(gtdict):
    filter = False
    if len(gtdict['mapref_ref']) > 1 and gtdict['mapref_ref'] != 'NA':
        filter = True
    for allele in gtdict['maprefalleles'].split(';'):
        if len(allele) > 1 and allele != 'NA':
            filter = True
    #if len(gtdict['superts_ref']) > 1 and gtdict['superts_ref'] != 'NA':
        #filter = True
    #for allele in gtdict['supertsalleles'].split(';'):
        #if len(allele) > 1 and allele != 'NA':
            #filter = True
    return filter      


if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='filters gtype intersection file and creates an SNV, bi-alleleic only table')
    parser.add_argument('-gt','--gtable',dest='genotypes',type=str,help='table of mapref and superts genotype intersections')
    opts = parser.parse_args()

    fopen=open(opts.genotypes,'r')
    fout=open('maprefbiallelicfiltered_%s' % basename(opts.genotypes),'w')
    fields = fopen.readline().strip().split('\t')
    fout.write('%s\n' % '\t'.join(fields))
    for line in fopen:
        linelist = line.strip().split('\t')
        linedict = dict(zip(fields,linelist))
        filter_eval = biallelic_snv_filter(linedict)
        if filter_eval == False:
            fout.write(line)

    fout.close()
