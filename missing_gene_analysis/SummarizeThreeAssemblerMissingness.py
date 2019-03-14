import argparse

def BuildExpressionDict(matrixin,tpmthreshold):
    fopen = open(matrixin,'r')
    fopen.readline()
    expr_dict = {}
    for line in fopen:
        gene,tpm = line.strip().split()
        if float(tpm) >= tpmthreshold:
            expr_dict[gene] = tpm
    return expr_dict        
    
def ParseMissingnessTable(intable,pconsdict = None):
    if pconsdict == None:
        pconsdict = {}
    code_dict = {'refexpressed':'P','missing':'A'}
    tabledict = {}
    fopen = open(intable,'r')
    fopen.readline()
    for line in fopen:
        geneid,pcons_mean,status,sample = line.strip().split()
        tabledict[geneid] = code_dict[status]
        pconsdict[geneid] = pcons_mean
    return tabledict,pconsdict    
        

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='pipeline for merging missingness by phastcons data at the gene level')
    parser.add_argument('-tr','--trinity',dest='tr',type=str,help='trinity ts coverage data')
    parser.add_argument('-sh','--shannon',dest='sh',type=str,help='shannon ts coverage data') 
    parser.add_argument('-bp','--binpacker',dest='bp',type=str,help='binpacker ts coverage data')
    parser.add_argument('-generefx','--mapref_gene_expression',dest='generefx',type=str,help='expr matrix for map ref')
    parser.add_argument('-fout','--outtable',dest='fout',type=str,help='outfile name')
    parser.add_argument('-sample','--sample_name',dest='sample',type=str,help='sample name')
    parser.add_argument('-tpm','--mintpm',dest='tpm',type=float,help='min tpm threshold')
    opts = parser.parse_args() 


    mapref_expr = BuildExpressionDict(opts.generefx,opts.tpm)

    trinitydict,pconsdict = ParseMissingnessTable(opts.tr)
    shannondict,pconsdict = ParseMissingnessTable(opts.sh,pconsdict=pconsdict)
    bpdict,pconsdict = ParseMissingnessTable(opts.bp,pconsdict=pconsdict)
  
    outtable = open(opts.fout,'w')
    outtable.write('ens_id\ttrinity\tshannon\tbinpacker\tTPM\tphastcons\tsample\n')
    for gene in mapref_expr:
        outtable.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gene,trinitydict[gene],shannondict[gene],bpdict[gene],mapref_expr[gene],pconsdict[gene],opts.sample))

    outtable.close()
