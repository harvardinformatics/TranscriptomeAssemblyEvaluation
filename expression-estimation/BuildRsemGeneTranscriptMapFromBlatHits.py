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
        ref_map[alist[0]]=alist[1]
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

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='pipeline for generating sundry transcriptome assembly coverage statistics')
    parser.add_argument('-psl','--psl_infile',dest='hits',type=str,help='concatenated psl file bits from job array, w/o headers')
    parser.add_argument('-minc','--min_query_coverage',dest='mincov',default=0.0,type=float,help='proportion of the query covered by matches in the alignment')
    parser.add_argument('-pid','--percent_identity',dest='pid',default=0.95,type=float,help='min identity in the alignment filtering threshold')
    parser.add_argument('-refmap','--reference-gene-ts-map',type=str,dest='refmap',help='mappings of reference tscripts to genes')
    parser.add_argument('-c','--assembly-contig-file',type=str,dest='contigs',help='list of de novo assembly contig ids')
    parser.add_argument('-mapout','--gene_transcript_map_out',dest='outmap',type=str,help='output gene transcript map')
    opts = parser.parse_args()

    contigs_handle = open(opts.contigs,'r')
    contig_list = []
    for line in contigs_handle:
        contig_list.append(line.strip())

    refmap = build_ref_ts_gene_map(opts.refmap)
    
    queries_to_targets = {}
    psl_open = open(opts.hits,'r')

    for hit in psl_open:
        values=hit.strip().split()
        hit_dict=dict(zip(psl_keys,values))
        hit_dict=psl_addstats(hit_dict)
        if hit_dict['pid']>=opts.pid and hit_dict['query_coverage']>=opts.mincov:
            # this strips the suffix off of Ensembl CDS entries #
            # causing mismatch with transcript ids #
            hit_dict['tName']=hit_dict['tName'].split('.')[0] 
            hit_dict=build_alignment_intervals(hit_dict)
            if hit_dict['qName'] in queries_to_targets:
                if hit_dict['target_coverage'] > queries_to_targets[hit_dict['qName']]['target_coverage']:
                    queries_to_targets[hit_dict['qName']]['tName'] = hit_dict['tName']
                    queries_to_targets[hit_dict['qName']]['target_coverage'] = hit_dict['target_coverage']
            else:
                queries_to_targets[hit_dict['qName']] = {}
                queries_to_targets[hit_dict['qName']]['tName'] = hit_dict['tName']
                queries_to_targets[hit_dict['qName']]['target_coverage'] = hit_dict['target_coverage']

    mapout = open(opts.outmap,'w')
    for contig in contig_list:
        if contig in queries_to_targets:
            mapout.write('%s\t%s\n' % (refmap[queries_to_targets[contig]['tName']],contig))
        else:
            mapout.write('%s\t%s\n' % ('_'.join(contig.split('_')[:-1]),contig))

    mapout.close()
    
                    
