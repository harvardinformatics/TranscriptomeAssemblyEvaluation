"""
Created by adam h. freedman (adamfreedman@fas.harvard.edu)
This script takes a psl file from a BLAT search of de novo
transcriptome assembly nucleotide sequences against a set
of reference CDS or cDNA sequences, e.g. from Ensembl,
and generates csv files that summarize coverage of reference
sequences for best hits and summed coverage, as a function
of reference transcript expression (TPM).
""" 

import argparse
from numpy import mean,median
from interval import interval, inf, imath
from collections import defaultdict

psl_keys=['matches','misMatches','repMatches','nCount','qNumInsert','qBaseInsert','tNumInsert','tBaseInsert','strand','qName','qSize','qStart','qEnd','tName','tSize','tStart','tEnd','blockCount','blockSizes','qStarts','tStarts']

def build_ref_ts_gene_map(mapfile):
    """
    takes a table mapping reference transcript
    to reference gene id and returns a dictionary
    with ts as key, gene as value
    """
    ref_map={}
    a=open(mapfile,'r')
    for aline in a:
        alist=aline.strip().split()
        ref_dict[alist[0]]=alist[1]
    return ref_map

def psl_addstats(psl_dict):
    """
    takes a psl file from a BLAT or similar search and
    generates alignment 'hitscore', query and target coverage
    of the alignment using matches and mismatches, percent identity (pid)
    """
    query_coverage=(float(psl_dict['matches'])+float(psl_dict['repMatches']))/float(psl_dict['qSize'])
    target_coverage=(float(psl_dict['matches'])+float(psl_dict['repMatches']))/float(psl_dict['tSize'])
    pid=float(psl_dict['matches'])/(float(psl_dict['matches'])+float(psl_dict['repMatches'])+float(psl_dict['misMatches']))
    psl_dict['query_coverage']=query_coverage
    psl_dict['target_coverage']=target_coverage
    psl_dict['pid']=pid
        
    return psl_dict

def build_expression_dict(exprmatrix,expcol=1,header=True):
    expr_dict={}
    fopen=open(exprmatrix,'r')
    if header==True:
        fopen.readline()
    for line in fopen:
        linelist=line.strip().split()
        expr_dict[linelist[0]]=float(linelist[expcol])
    return expr_dict

def build_alignment_intervals(psl_dict):
    """
    build list of intervals, where the interval
    start and end are the bounding coordinates
    in an interval, then creates a merged interval
    object, adds the interval list and merged interval
    to psl_dict and returns updated psl_dict
    """
    block_intervals=[]
    for i in range(int(psl_dict['blockCount'])):
        interval_start=int(psl_dict['tStarts'].split(',')[i])
        interval_end=int(psl_dict['tStarts'].split(',')[i]) + int(psl_dict['blockSizes'].split(',')[i])
        new_interval=interval[interval_start,interval_end]
        block_intervals.append(new_interval)

    merged_interval=interval()
    for newinterval in block_intervals:
        merged_interval = merged_interval | newinterval
    
    psl_dict['interval_list']=block_intervals
    psl_dict['merged_interval']=merged_interval
    psl_dict['coverage_from_interval']=calculate_target_coverage(merged_interval,psl_dict['tSize'])
   
    return psl_dict

def calculate_target_coverage(intervals,targetlength):
    bases=0
    for newinterval in intervals:
        bases=bases+(newinterval[1]-newinterval[0])
    return bases/float(targetlength)

def build_besthit_per_target(target_to_queries):
    """
    this function takes a default dict with targets as keys
    and lists of blat hit dictionaries as values, and extracts
    the best hit per target. filtering on pid and mincov is already
    done before passing target_to_queries, such that starting by adding
    the first entry in any list does not lead to initiation of best per target
    dict with hits that should have been filtered
    """
    besthit_per_target={}
    for target in target_to_queries:
        besthit_per_target[target]=target_to_queries[target][0]
        if len(target_to_queries[target])>1:
            for i in range(len(target_to_queries[target])-1):
                if target_to_queries[target][i+1]['coverage_from_interval']>besthit_per_target[target]['coverage_from_interval']: 
                    besthit_per_target[target]=target_to_queries[target][i+1]
    return besthit_per_target    
    

def calculate_summed_coverage(target_to_queries):
    summed_coverage_dict={}
    for targetkey in target_to_queries:
        target_length=int(target_to_queries[targetkey][0]['tSize'])
        queries=[]
        totalinterval=interval()        
        for querydict in target_to_queries[targetkey]:
            queries.append(querydict['qName'])    
            totalinterval=totalinterval | querydict['merged_interval']
        sumcov=calculate_target_coverage(totalinterval,target_length)
        queries_string=';'.join(queries)
        summed_coverage_dict[targetkey]={'coverage':sumcov,'queries':queries_string,'numqueries':len(queries)}
    return summed_coverage_dict

def build_ts_gene_map(mapfilehandle):
    """
    assumes a tab separated file, first
    column == tsid, 2nd column == geneid
    """
    ts_to_gene_dict = {}
    for line in mapfilehandle:
        tsid,geneid =  line.strip().split('\t')
        ts_to_gene_dict[tsid] = geneid
    return ts_to_gene_dict


def calc_tscript_tpm_weights(ts_expr_matrix,ts_gene_map):
    expr_group_by_gene = {}
    for tscript in ts_gene_map.keys():
        if tscript in ts_expr_matrix:
            if ts_gene_map[tscript] in expr_group_by_gene:
                expr_group_by_gene[ts_gene_map[tscript]][tscript] = float(ts_expr_matrix[tscript])
            else:
                expr_group_by_gene[ts_gene_map[tscript]] = {}
                expr_group_by_gene[ts_gene_map[tscript]][tscript] = float(ts_expr_matrix[tscript])
    for gene in expr_group_by_gene:
        tpm_sum = 0
        for tscript in expr_group_by_gene[gene]:
            tpm_sum += expr_group_by_gene[gene][tscript]
        if tpm_sum > 0:
            for tscript in expr_group_by_gene[gene]:
                expr_group_by_gene[gene][tscript] = expr_group_by_gene[gene][tscript]/float(tpm_sum)    
    
    return expr_group_by_gene        

def calculate_weighted_coverage(weight_dict,sumcovdict):
    gene_cov_dict = {}
    for gene in weight_dict:
        gene_coverage = 0
        for ts in weight_dict[gene]:
            if ts in sumcovdict:
                gene_coverage += weight_dict[gene][ts] * sumcovdict[ts]['coverage']    
        gene_cov_dict[gene] = gene_coverage 
    return gene_cov_dict

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='pipeline for generating sundry transcriptome assembly coverage statistics')
    parser.add_argument('-psl','--psl_infile',dest='hits',type=str,help='concatenated psl file bits from job array, w/o headers')
    parser.add_argument('-minc','--min_query_coverage',dest='mincov',default=0.0,type=float,help='proportion of the query covered by matches in the alignment')
    parser.add_argument('-pid','--percent_identity',dest='pid',default=0.95,type=float,help='min identity in the alignment filtering threshold')
    parser.add_argument('-isorefx','--reference_isoform_matrix',dest='isorefx',type=str,help='expression matrix for isoform map-to-reference')
    parser.add_argument('-generefx','--reference_gene_matrix',dest='generefx',type=str,help='expression matrix for gene map-to-reference')
    parser.add_argument('-o','--outfile_prefix',dest='outprefix',type=str,help='prefix for output files summarizing coverage')
    parser.add_argument('-gtmap','--gene_transcript_map',dest='gtmap',type=str,help='gene transcript map')
    opts = parser.parse_args()

    psl_open = open(opts.hits,'r')
    maphandle = open(opts.gtmap,'r')

    gene_ts_dict = build_ts_gene_map(maphandle)
    
    # create dict for storing coverage intervals at the gene level
    gene_interval_dict = defaultdict()

    # build reference expression dictionaries
    target_expression = build_expression_dict(opts.isorefx)
    gene_expression_dict = build_expression_dict(opts.generefx)

    # get tpm weights for ref tscripts
    weight_dict = calc_tscript_tpm_weights(target_expression,gene_ts_dict)
    # create dictionary with targets as keys 
    target_to_queries=defaultdict(list)

    for hit in psl_open:
        values=hit.strip().split()
        hit_dict=dict(zip(psl_keys,values))
        hit_dict=psl_addstats(hit_dict)
        if hit_dict['pid']>=opts.pid:
            hit_dict['tName']=hit_dict['tName'].split('.')[0] # this strips the suffix off of Ensembl CDS entries causing mismatch with transcript ids
            hit_dict=build_alignment_intervals(hit_dict)
            if hit_dict['query_coverage']>=opts.mincov and hit_dict['pid']>=opts.pid:
                target_to_queries[hit_dict['tName']].append(hit_dict)
                # update gene-level coverage intervals
                genetarget = gene_ts_dict[hit_dict['tName']]
                if genetarget in gene_interval_dict:
                    gene_interval_dict[genetarget] = gene_interval_dict[genetarget] | hit_dict['merged_interval']
                else:
                    gene_interval_dict[genetarget] = hit_dict['merged_interval']
    besthit_per_target = build_besthit_per_target(target_to_queries)   
    bestout=open(opts.outprefix+'besthit_coverage.csv','w')
    bestout.write('target,query,coverage,TPM\n') 
    bestmis=open(opts.outprefix+'_best.missing','w')
    for targetkey in target_expression:
        if target_expression[targetkey] > 0:
            if targetkey in besthit_per_target:
                bestout.write('%s,%s,%s,%s\n' % (targetkey,besthit_per_target[targetkey]['qName'],besthit_per_target[targetkey]['coverage_from_interval'],target_expression[targetkey]))
            else:
                bestout.write('%s,%s,%s,%s\n' % (targetkey,'NA',0,target_expression[targetkey]))
                bestmis.write('%s\n' % targetkey)
    
    bestout.close()    
    bestmis.close()
    
    summed_coverage_per_target = calculate_summed_coverage(target_to_queries)
    sumout=open(opts.outprefix+'_summedcoverage.csv','w')
    sumout.write('target,numqueries,coverage,TPM,querylist\n')
    summis=open(opts.outprefix+'summed.missing','w')
    for targetkey in target_expression:
        if target_expression[targetkey] > 0:
            if targetkey in summed_coverage_per_target:
                sumout.write('%s,%s,%s,%s,%s\n' % (targetkey,summed_coverage_per_target[targetkey]['numqueries'],summed_coverage_per_target[targetkey]['coverage'],target_expression[targetkey],summed_coverage_per_target[targetkey]['queries']))
            else:
                sumout.write('%s,%s,%s,%s,%s\n' % (targetkey,0,0,target_expression[targetkey],'NA'))
                summis.write('%s\n' % targetkey)

    sumout.close()
    summis.close()
   
    weighted_gene_coverage = calculate_weighted_coverage(weight_dict,summed_coverage_per_target)
    gmissing = open('missing_genes.txt','w') 
    gene_coverage_out = open('%s_genecoverage.csv' % opts.outprefix,'w')
    gene_coverage_out.write('geneid,refTPM,weighted_coverage\n')
    for gene in gene_expression_dict.keys():
        if gene_expression_dict[gene] > 0:
            if gene in gene_interval_dict:    
                wt_cov = weighted_gene_coverage[gene]
            else:
                gmissing.write('%s\n' % gene)
                wt_cov = 0
            gene_coverage_out.write('%s,%s,%s\n' % (gene,gene_expression_dict[gene],wt_cov))

    gene_coverage_out.close()    
    gmissing.close()
