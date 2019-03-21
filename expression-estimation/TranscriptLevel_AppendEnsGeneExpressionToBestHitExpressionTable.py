# /n/holylfs/LABS/informatics/adamf/refs/mus/GRCm38/complete_EnsTsEnsGMap.txt
import argparse
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="append ensembl gene level expression data to best hit to ts tables")
    parser.add_argument('-m','--gene_ts_map',dest='gtmap',type=str,help='ensembl mapping of gene to ts')
    parser.add_argument('-e','--best_hit_expr_table',dest='exprtable',type=str,help='best hit vs ensembl ts expression csv')
    parser.add_argument('-gex','--mapref_gene_expression',dest='genexpr',type=str,help='ensembl gene level expression data')
    opts = parser.parse_args()
    
    ts_gene_map = open(opts.gtmap,'r')
    ts_gene_dict = {}
    for line in ts_gene_map:
        linelist = line.strip().split()
        ts_gene_dict[linelist[0]] = linelist[1]

    genedata = open(opts.genexpr,'r')
    gene_expr = {}
    genekeys = genedata.readline().strip().split('\t')
    for line in genedata:
        linelist = line.strip().split('\t')
        linedict = dict(zip(genekeys,linelist))
        gene_expr[linedict['gene_id']] = {'expected_count':linedict['expected_count'],'TPM':linedict['TPM']}
         
        
    fopen = open(opts.exprtable,'r')
    fout = open('wEnsemblGene_%s' % opts.exprtable,'w')
    header = fopen.readline().strip()
    header = '%s,EnsGene,refgene_expected_count,refgene_TPM\n' % header
    fout.write(header) 
    for line in fopen:
        linelist = line.strip().split(',')
        gene = ts_gene_dict[linelist[0]]
        fout.write('%s,%s,%s,%s\n' % (','.join(linelist),gene,gene_expr[gene]['expected_count'],gene_expr[gene]['TPM']))
        
    fout.close()          
