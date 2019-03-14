import os
from subprocess import Popen,PIPE
import argparse

def DependencyPathTest(dependency):
    cmd = 'which %s' % dependency
    testout,testerr=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE).communicate()
    return testout
    

def BuildGeneBed(geneid,featurebed):
    cmd = 'grep %s %s > %s.bed' % (geneid,featurebed,geneid)
    makebed = Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    makeout,makeerr = makebed.communicate()
    if makebed.returncode==0:
        sortcmd = 'sortBed -i %s.bed > sorted_%s.bed' % (geneid,geneid)
        sortbed = Popen(sortcmd,shell=True,stderr=PIPE,stdout=PIPE)
        sortout,sorterr=sortbed.communicate()
        if sortbed.returncode==0:
            mergecmd = 'mergeBed -i sorted_%s.bed > merge_sorted_%s.bed' % (geneid,geneid)
            mergebed = Popen(mergecmd,shell=True,stderr=PIPE,stdout=PIPE)
            mergeout,mergeerr=mergebed.communicate()
            if mergebed.returncode!=0:
                print "MergeERROR"
                return mergeerr
        else:
            print "SortERROR"
            return sorterr        

    else:
        print 'Make Gene.bed Error'
        return makeerr

def TabixRunner(genebed,pcons):
    gene_pcons = 0
    num_bases = 0
    gene_intervals = open(genebed,'r')
    for interval in gene_intervals:
        chrom,start,end = interval.strip().split()
        tabix_cmd = 'tabix %s %s:%s-%s' % (pcons,chrom,start,end)
        tabix_run =  Popen(tabix_cmd,shell=True,stderr=PIPE,stdout=PIPE)
        tabixout,tabixerr = tabix_run.communicate()
        if tabix_run.returncode != 0:
            print 'TABIX ERROR'
            print 'tabix cmd', tabix_cmd

        else:
            interval_blocks = tabixout.split('\n')[:-1]
            if len(interval_blocks)>0:
                for block in interval_blocks:
                    chrom,start,end,score = block.strip().split('\t')
                    num_bases+=int(end)-int(start)
                    gene_pcons+=float(score)*(int(end)-int(start))
            #else:
                #print 'empty interval',interval
    if num_bases>0:
        meanscore = gene_pcons/float(num_bases)
    else:
        #print 'no scores for gene',genebed
        meanscore = '-9'
    return meanscore

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='pipeline for obtaining gene-level phastcons scores')
    parser.add_argument('-annots','--gene_annot_bed',dest='featurebed',type=str,help='bed file of gene id and genomic subelement')
    parser.add_argument('-genes','--ensembl_gene_list',dest='genes',type=str,help='list of ensembl genes to iterate over')
    parser.add_argument('-pcon','--phastcons_bedgraph',dest='pcons',type=str,help='name of tabix indexed phastcons bedgraph file')
    parser.add_argument('-o','--gene_cons_out',dest='outfile',type=str,help='name of file describing mean phastcons per gene')
    opts = parser.parse_args() 

    # verify bedtools and tabix in path
    dependencies = ['bedtools','tabix']
    for dependency in dependencies:
        path_eval = DependencyPathTest(dependency)
        if path_eval == '':
            return Exception("Exception: %s not in path" % dependency) 

    # load gene list
    genes=open(opts.genes,'r')
    gene_list = []
    for gene in genes:
        gene_list.append(gene.strip())

    # open outfile #
    fout = open(opts.outfile,'w')
    # get average phastcons scores per gene 
    counter=0
    for gene in gene_list:
        counter+=1
        if counter%1000 == 0:
            print 'processing gene %s' % gene
        makebederr = BuildGeneBed(gene,opts.featurebed)
        mergeid='merge_sorted_%s.bed' % gene
        mean_gene_phastcons = TabixRunner(mergeid, opts.pcons)
        fout.write('%s\t%s\n' % (gene,mean_gene_phastcons))
        os.system('rm *%s.bed' % gene)
    fout.close()


