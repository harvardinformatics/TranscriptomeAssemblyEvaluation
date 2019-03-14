import argparse
from sets import Set

def BuildProteinHitSet(blastin):
    """
    creates sets of Trinity genes whose
    transdecoder orfs have blastp hits
    """
    fopen = open(blastin,'r')
    hitset = Set()
    for line in fopen:
        hit = line.strip().split()[0].split(':')[0]
        hitset.add(hit) 
    return hitset    


if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='convert SuperTranscript genotypes with no protein blast hit to NAs, and add flag column')
    parser.add_argument('-gt','--gtable',dest='genotypes',type=str,help='table of mapref and superts genotype intersections')
    parser.add_argument('-bhits','-blast-results',dest='blast',type=str,help='name of blastp(x) output for Trinity contigs')
    opts = parser.parse_args()
    trinity_gene_protein_hits = BuildProteinHitSet(opts.blast)
  
    gtopen = open(opts.genotypes,'r')
    fields = gtopen.readline().strip().split('\t')
    gtout = open('proteinhitfiltered_%s' % opts.genotypes,'w')
    gtout.write('%s\tproteinhit\n' % '\t'.join(fields))  
    

    for line in gtopen:
        linelist = line.strip().split('\t')
        linedict = dict(zip(fields,linelist))
        found = 0
        for st in linedict['SToverlaps'].split(';'):
            if st.split(':')[0] in trinity_gene_protein_hits:
                found+=1
        if found > 0:
            gtout.write('%s\tY\n' % '\t'.join(linelist))
        else:
            linedict['supertspositions'] = 'NA'
            linedict['supertsref'] = 'NA'
            linedict['supertsalleles'] = 'NA'
            linedict['superts_ad'] = 'NA'
            linedict['mean_supercov'] = '0'
            linedict['median_supercov'] = '0'
            linedict['SToverlaps'] = 'NA'
            if linedict['supertsalleles'] == 'NA' and linedict['maprefalleles'] == 'NA':
                pass
            else:
                vals_update = []
                for field in fields:
                    vals_update.append(linedict[field])
                gtout.write('%s\tN\n' % '\t'.join(vals_update))

    gtout.close()
