import argparse
from collections import defaultdict
from sets import Set

def kallisto_to_dict(abundancefile):
    fopen = open(abundancefile,'r')
    kal_dict = {}
    colnames = fopen.readline().strip().split()
    for line in fopen:
        linelist = line.strip().split()
        linedict = dict(zip(colnames,linelist))
        linedict['eff_length'] = float(linedict['eff_length'])
        linedict['length'] = int(linedict['length'])
        linedict['tpm'] = float(linedict['tpm']) 
        linedict['est_counts'] = float(linedict['est_counts'])
        kal_dict[linedict['target_id']] = linedict
        kal_dict[linedict['target_id']].pop('target_id')
    return kal_dict

def build_besthit_tscript_dict(hittable):
    hit_dict = {}
    fopen = open(hittable,'r')
    colnames = fopen.readline().strip().split(',')
    for line in fopen:
        linelist = line.strip().split(',')
        linedict = dict(zip(colnames,linelist))
        hit_dict[linedict['EnsTsId']] = linedict['BestHitId']
    return hit_dict
   
def build_gene_to_contigs_map(maptable):
    fopen = open(maptable,'r')
    map_dict = defaultdict(list)
    for line in fopen:
        gene,contig = line.strip().split()
        map_dict[gene].append(contig)
    return map_dict


def GeneTpmsFromContigSet(tsmap,kaldict,totalcount):
    gene_dict = {}
    for gene in tsmap:
        flag = 0 
        for ts in tsmap[gene]:
            if ts in kaldict:
                flag+=1
        if flag > 0:
            gene_dict[gene] = {}
            tpm = 0
            estcount = 0
            for ts in tsmap[gene]:
                tpm += kaldict[ts]['tpm']
                estcount += kaldict[ts]['est_counts']
            gene_dict[gene]['tpm'] = tpm
            gene_dict[gene]['estcount'] = estcount
    
    return gene_dict

def GetTotalCount(kaldict):
    totalcount = 0
    for contig in kaldict:
        totalcount += kaldict[contig]['est_counts']
    return totalcount        

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="append ensembl gene level expression data to best hit to ts tables")
    parser.add_argument('-r','--mapref-quant',dest='mref',type=str,help='kallisto abundance.tsv file for mapref')
    parser.add_argument('-d','--denovof-quant',dest='denovo',type=str,help='kallisto abundance.tsv file for denovo transcriptome assembly')
    parser.add_argument('-mrmap','--mapref_ts_gene_map',dest='mrmap',type=str,help='file mapping ref tscripts to genes')
    parser.add_argument('-demap','--denovo_ts_gene_map',dest='demap',type=str,help='file mapping contigs to genes')
    parser.add_argument('-p','--outfile-prefix',dest='prefix',type=str,help='prefix for outtables')
    opts = parser.parse_args()

    mr_kaldict = kallisto_to_dict(opts.mref)
    mr_totalcount = GetTotalCount(mr_kaldict)
    for tscript in mr_kaldict:
        ts_lstpm = LengthScaleTpmFromContig(mr_kaldict[tscript]['eff_length'],mr_kaldict[tscript]['tpm'],mr_totalcount)

    mr_map = build_gene_to_contigs_map(opts.mrmap)
    mr_genelevel = GeneTpmsFromContigSet(mr_map,mr_kaldict,mr_totalcount) 
 
    denovo_kaldict = kallisto_to_dict(opts.denovo)
    denovo_totalcount = GetTotalCount(denovo_kaldict)
    for contig in denovo_kaldict:
    
    denovo_map = build_gene_to_contigs_map(opts.demap)
    denovo_genelevel = GeneTpmsFromContigSet(denovo_map,denovo_kaldict,denovo_totalcount)

    all_genes = Set()
    for gene in mr_genelevel:
        all_genes.add(gene)
    for gene in denovo_genelevel:
        if 'ENS' in gene:
            all_genes.add(gene)
 
    genetable = open('%s_MRvsDenovo_Expression.csv' % opts.prefix,'w')
    genetable.write('geneid,mr_count,mr_tpm,denovo_count,denovo_tpm\n')
    for gene in all_genes:
        if gene in mr_genelevel and gene in denovo_genelevel:
            genetable.write('%s,%s,%s,%s,%s\n' % (gene,mr_genelevel[gene]['estcount'],mr_genelevel[gene]['tpm'],denovo_genelevel[gene]['estcount'],denovo_genelevel[gene]['tpm']))
        elif gene in mr_genelevel and gene not in denovo_genelevel:
            genetable.write('%s,%s,%sNA,NA\n' % (gene,mr_genelevel[gene]['estcount'],mr_genelevel[gene]['tpm']))
        elif gene not in mr_genelevel and gene in denovo_genelevel:
            genetable.write('%s,NA,NA,%s,%s\n' % (gene,denovo_genelevel[gene]['estcount'],denovo_genelevel[gene]['tpm']))
        else:
           print "not properly cagegorized"

    genetable.close()
