import argparse
"""
/n/holylfs/LABS/informatics/adamf/DeNovoTranscriptomeEvaluation/contig_overlap_analysis/python_code/mus_single_isoform_genes.txt
"""

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Extract entries from BLAT psl file for single isoform genes')
    parser.add_argument('-g','--single_isoform_gene_list',dest='genes',type=str,help='gene/tscript tsv table')
    parser.add_argument('-psl','--psl_infile',dest='psl',type=str,help='concatenated psl file bits from blat all-by-all job array, w/o headers')
    opts = parser.parse_args()    
    genets=open(opts.genes,'r') 
    tslist=[]
    for line in genets:
        linelist=line.strip().split()
        tslist.append(linelist[1])

    pslbed=open(opts.psl,'r')
    filtbed=open('singleisoform_%s' % opts.psl,'w')

    for line in pslbed:
        linelist=line.strip().split()
        if linelist[0] in tslist:
            filtbed.write(line)

    filtbed.close()

