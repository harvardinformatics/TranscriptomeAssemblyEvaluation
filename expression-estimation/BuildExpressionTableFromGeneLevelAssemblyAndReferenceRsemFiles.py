import argparse
from sets import Set

def make_dict_from_rsem_out(rsemoutfile):
    fopen=open(rsemoutfile,'r')
    keys=fopen.readline().split()
    fdict={}
    for line in fopen:
        values=line.strip().split()
        key_val_dict=dict(zip(keys,values))
        fdict[key_val_dict['gene_id']]=key_val_dict
        del(fdict[key_val_dict['gene_id']]['gene_id'])
     
    return fdict  
    

def query_expression_dict(key,dict):
    try:
        querydict=dict[key]
    except:
        querydict={'length':'NA','effective_length':'NA','expected_count':'0.00','TPM':'0.00','FPKM':'0.00'} 
        
    queryvalues=[querydict['length'],querydict['effective_length'],querydict['expected_count'],querydict['TPM'],querydict['FPKM']]
    return queryvalues

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="options for creating table of de novo assenbly isoform expression estimates from 'best hits' csv and map to ref expression estimates")
    parser.add_argument('-asex','--denovo_expression_matrix',dest='denovoexpr',type=str,help='denovo expression matrix')
    parser.add_argument('-refex','--reference_expressIon_matrix',dest='refexpr',type=str,help='reference expression matrix')
    parser.add_argument('-o','--outfile',dest='outfile',type=str,help="name of output file")
    opts = parser.parse_args()
    
    
    refdict=make_dict_from_rsem_out(opts.refexpr)
    
    denovodict=make_dict_from_rsem_out(opts.denovoexpr)   
    gene_keys=set(refdict.keys()).union(set(denovodict.keys()))
    fout=open(opts.outfile,'w')
    fout.write('geneid,ref_len,ref_efflen,ref_expcount,refTPM,refFPKM,denovo_len,denovo_efflen,denovo_expcount,denovoTPM,denovoFPKM\n')
    
    for key in gene_keys:
        refquery=query_expression_dict(key,refdict)
      
        denovoquery=query_expression_dict(key,denovodict)
      
        fout.write('%s,%s,%s\n' % (key,','.join(refquery),','.join(denovoquery)))


             
fout.close()             
                 
