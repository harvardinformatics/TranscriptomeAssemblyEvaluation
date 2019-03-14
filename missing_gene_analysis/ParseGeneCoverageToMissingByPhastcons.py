import argparse

def BuildGenePhastConstDict(pconstable):
    pcon_dict = {}
    popen=open(pconstable,'r')
    for line in popen:
        gene,score = line.strip().split()
        pcon_dict[gene] = max(0,float(score))
    return pcon_dict

def ParsePresentMissingGenes(covfile,phastconsdict,outhandle,tpm,sample):
    covhandle = open(covfile,'r')
    covhandle.readline()
    background = open('%s.background' % covfile,'w')
    missing = open('%s.missing' % covfile,'w')
    for line in covhandle:
        geneid,refTPM,weighted_coverage = line.strip().split(',')
        refTPM = float(refTPM)
        weighted_coverage = float(weighted_coverage)
        if refTPM >= tpm:
            if geneid in phastconsdict:
                pcons = phastconsdict[geneid]
            else:
                pcons = 0    
            if weighted_coverage == 0:
                missing.write('%s\n' % geneid)
                background.write('%s\n' % geneid)
                flag = 'missing'
            else:
                background.write('%s\n' % geneid)
                flag = 'refexpressed'
            outhandle.write('%s\t%s\t%s\t%s\n' % (geneid,pcons,flag,sample))
    
    background.close()
    missing.close()
                                   
if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='pipeline for obtaining gene-level phastcons scores')    
    parser.add_argument('-gcov','--gene_coverage_table',dest='genecov',type=str,help='table of gene coverage generated with OmnibusAssemblyCoverageAnalysis.py')
    parser.add_argument('-genepcon','--phastcons_by_gene',dest='gene_pcons',type=str,help='mean phastcons by gene table')
    parser.add_argument('-o','--outtable',dest='outtable',type=str,help='outfile with phastcons for missing and present genes in an assembly')
    parser.add_argument('-sample','--sample_label',dest='sample',type=str,help='sample id for column header in out table')
    parser.add_argument('-tpmmin','--min_mapref_tpm',dest='tpm',type=float,help='min mapref TPM to consider gene expressed')
    opts = parser.parse_args() 
    
    outfile = open(opts.outtable,'w')
    outfile.write('geneid\tpcons_mean\tstatus\tsample\n')
    phastcons_by_gene = BuildGenePhastConstDict(opts.gene_pcons)
    ParsePresentMissingGenes(opts.genecov,phastcons_by_gene,outfile,opts.tpm,opts.sample)
    outfile.close()
