import argparse
from CalculateConcordanceMetrics_modular import FalseNegative,BoolFlagGenotypeSetConcordance
from sets import Set

def ExtractMafIfHet(mr_dict,ploidy):
    allele_list = [int(allele) for allele in mr_dict['gtypedata'].split(':')[0].split('/')]
    alleles = Set(allele_list)
    if len(alleles) == 2:
        hetflag = True
        maf = min(allele_list.count(list(alleles)[0]),allele_list.count(list(alleles)[1]))/float(ploidy)    
    else:
        hetflag = False
        maf = 'NA'

    return hetflag,maf

def EvaluateSuperTscriptAtMapRefHetSite(snp_dict):
    if snp_dict['supertsalleles'] == 'NA':
        stclass = 'FN'
    elif len(snp_dict['supertsalleles'].split(';')) == 1:
        stclass = 'FN'
    elif Set(snp_dict['supertsalleles']) == Set(snp_dict['maprefalleles']):
        stclass='concordant'
    elif len(snp_dict['supertsalleles'].split(';')) > 2 or Set(snp_dict['supertsalleles']) != Set(snp_dict['maprefalleles']):
        stclass = 'error'
    else:
        stclass = 'NA'
    return stclass


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Calculates map-ref het MAF for SuperTranscript FN, and concordant genotype calls')
    parser.add_argument('-gtinter','--genotype-intersection',dest='gtintersect',type=str,help='table of mapref and superts genotype intersections')
    parser.add_argument('-mrbed','--mapref-genotype-bed', dest='mrbed',type=str,help='bed file of map-to-ref genotypes')
    parser.add_argument('-p','--ploidy',dest='ploidy',type=int,help='ploidy of sample pool')
    parser.add_argument('-o','--outfile',dest='outfile',type=str,help='name of output table')
    opts = parser.parse_args() 

    intersect_fields = ['snpid', 'genomicpositions', 'supertspositions', 'mapref_ref', 'maprefalleles', 'superts_ref', 'supertsalleles', 'maprefcov', 'mapref_ad', 'superts_ad','mean_supercov','median_supercov','SToverlaps','SToverlapalleles']
    
    mr_open=open(opts.mrbed,'r')
    mr_fields = ['chrom','start','end','ref','alt','qual','filter','info','gtype_fields','gtypedata']
    mr_het_maf_dict = {}
    for line in mr_open:
        mr_line_dict = dict(zip(mr_fields,line.strip().split('\t')))
        ishet,maf =ExtractMafIfHet(mr_line_dict,opts.ploidy)
        if ishet == True:
            mr_het_maf_dict['%s:%s' % (mr_line_dict['chrom'],mr_line_dict['end'])] = maf

    summary_out = open(opts.outfile,'w')
    summary_out.write('chrom\tpos\tmaf\tSTtype\n')
    intersect_open = open(opts.gtintersect,'r')
    for line in intersect_open:
        line_list = line.strip().split('\t')
        snp_dict =  dict(zip(intersect_fields,line_list))
        if len(snp_dict['maprefalleles'].split(';')) == 2:
            maf = mr_het_maf_dict[snp_dict['genomicpositions']]
            snpclass = EvaluateSuperTscriptAtMapRefHetSite(snp_dict)
            if snpclass == 'NA':
                print snp_dict
            else:
                summary_out.write('%s\t%s\t%s\t%s\n' % (snp_dict['genomicpositions'].split(':')[0],snp_dict['genomicpositions'].split(':')[1],maf,snpclass))
 
    summary_out.close()
