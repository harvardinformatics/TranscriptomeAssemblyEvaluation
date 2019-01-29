import argparse

def CountVariantsFromFilteredVcf(filt_variants_handle,ploidy):
    
    het_count = 0
    alt_homo_count = 0
    keys='CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	GENOTYPE'.split('\t')               
    ploidy_string = '/'.join(['1']*ploidy)
    bedout=open('het_variants.bed','w')
    for line in filt_variants_handle:
        if line[0] != '#' and "PATCH" not in line:
            gtype_dict = dict(zip(keys,line.strip().split('\t')))            
            if gtype_dict['FILTER'] == 'PASS' and len(gtype_dict['ALT']) == 1 and len(gtype_dict['REF']) == 1:
                if ploidy_string in gtype_dict['GENOTYPE']:
                    alt_homo_count+=1
                else:
                    het_count+=1
                    bedout.write('%s\t%s\t%s\n' % (gtype_dict['CHROM'],int(gtype_dict['POS'])-1,gtype_dict['POS'],))
    bedout.close()

    return het_count,alt_homo_count            

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='pipeline for generating sundry transcriptome assembly coverage statistics')
    parser.add_argument('-vvcf','--filt-variants-vcf',dest='vvcf',type=str,help='concatenated psl file bits from job array, w/o headers')
    parser.add_argument('-o','--outfile',dest='outfile',type=str,help='name of outfile to which counts are printed')
    parser.add_argument('-p','--ploidy',dest='ploidy',type=int,default=2,help='# of chromosomes in sample, i.e. 2 x # pooled individuals')
    opts = parser.parse_args()
    
    variants_vcf=open(opts.vvcf,'r')
    het_count,homo_alt_count = CountVariantsFromFilteredVcf(variants_vcf,opts.ploidy)
    fout = open(opts.outfile,'w')
    fout.write('number homo alt sites: %s\n' % homo_alt_count)
    fout.write('number het sites: %s\n' % het_count)
    fout.close()
    
