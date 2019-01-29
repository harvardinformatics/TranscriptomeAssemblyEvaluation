import argparse
from sets import Set
from IntersectMapRefandSuperTranscriptGenotypes import ReverseComplement,GenotypeLineParse,CollapseGenotype


if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Create gtype bed for visualizing gtypes')
    parser.add_argument('-gbed','--genotype_bed',dest='gtypebed',type=str,help='bed genotypes in genomic coordinates')
    parser.add_argument('-o','--output-bed',dest='outbed',help='output bed file')
    parser.add_argument('-s','--superts',action='store_true',help='boolean flag to indicate superts')
    opts=parser.parse_args()

    gtypes = open(opts.gtypebed,'r')
    fout = open(opts.outbed,'w')

    superts_fields = ['gchrom','gposzero','gpos','gstrand','depth','contigid','cposzero','cpos','id','ref','alt','qual','filter','info','gtformats','gtdata']
    mapref_fields = ['gchrom','gposzero','gpos','ref','alt','qual','filter','info','gtformats','gtdata']
    
    if opts.superts:
        opts.fields = superts_fields
    else:
        opts.fields = mapref_fields
        
    for line in gtypes:
        linedict,gtypedict,alleles,ref_allele = GenotypeLineParse(line,opts.fields)
        if opts.superts == True and linedict['gstrand'] == '-':
            ref_allele = ReverseComplement(ref_allele)
            revcomp_alleles = []
            for allele in alleles:
                revcomp_alleles.append(ReverseComplement(allele))
            alleles = revcomp_alleles
        fout.write('%s\t%s\t%s\t%s|%s\n' % (linedict['gchrom'],linedict['gposzero'],linedict['gpos'],ref_allele,';'.join(alleles)))
        
    fout.close()
