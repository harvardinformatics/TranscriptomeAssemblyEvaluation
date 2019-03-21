import argparse

def parse_rsem_output(infile):
        fopen = open(infile,'r')
        keys = fopen.readline().strip().split()
        expr_dict = {}
        for line in fopen:
            linelist=line.strip().split()
            expr_dict[linelist[0]] =  dict(zip(keys,linelist))
        return expr_dict 

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description="options for creating table of de novo assenbly isoform expression estimates from 'best hits' csv and map to ref expression estimates")
    parser.add_argument('-bmap','--refts_besthit_map',dest='bestmap',type=str,help='best hit per ref transcript csv file, with expression data')
    parser.add_argument('-refex','--reference_isoform_expression',dest='refexpr',type=str,help='reference expression rsem isoform output')
    parser.add_argument('-denovex','--denovo_isoform_expression',dest='denovoexpr',type=str,help='de novo expression rsem isoform output')
    parser.add_argument('-o','--outfile',dest='outfile',type=str,help="name of output file")
    opts = parser.parse_args()
    
    ref_expr_dict = parse_rsem_output(opts.refexpr)
    denovo_expr_dict = parse_rsem_output(opts.denovoexpr)
    ts_to_hit_dict = {}
    hitmap = open(opts.bestmap,'r')
    for line in hitmap:
        linelist=line.strip().split(',')
        ts_to_hit_dict[linelist[0]] =  linelist[1]

    fout=open(opts.outfile,'w')
    fout.write('EnsTsId,BestHitId,EnsExpCount,EnsTPM,EnsEffLen,BestExpCount,BestTPM,BestEffLen\n')
    for refts in ref_expr_dict.keys():
        ens_exp_count = ref_expr_dict[refts]['expected_count']
        ens_tpm = ref_expr_dict[refts]['TPM']
        ens_eff_len =  ref_expr_dict[refts]['effective_length']

        if ts_to_hit_dict[refts] == 'NA':
            best_hit_id = 'NA'
            best_exp_count = 0
            best_tpm = 0
            best_eff_len = 0
        else:
            best_hit_id = ts_to_hit_dict[refts]
            best_exp_count = denovo_expr_dict[best_hit_id]['expected_count']
            best_tpm = denovo_expr_dict[best_hit_id]['TPM']
            best_eff_len = denovo_expr_dict[best_hit_id]['effective_length']
        if float(ens_exp_count) > 0 or float(best_exp_count) > 0:
            fout.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (refts,best_hit_id,ens_exp_count,ens_tpm,ens_eff_len,best_exp_count,best_tpm,best_eff_len))

fout.close()             




