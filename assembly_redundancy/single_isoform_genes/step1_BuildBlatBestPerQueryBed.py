import argparse

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

def BuildBestHitPerQueryDict(psl_handle):
    psl_keys=['matches','misMatches','repMatches','nCount','qNumInsert','qBaseInsert','tNumInsert','tBaseInsert','strand','qName','qSize','qStart','qEnd','tName','tSize','tStart','tEnd','blockCount','blockSizes','qStarts','tStarts']
    query_to_target_dict={}
    for line in psl_handle:
        psl_dict = dict(zip(psl_keys,line.strip().split()))
        psl_dict = psl_addstats(psl_dict)
       
        if psl_dict['query_coverage'] >= opts.mincov and psl_dict['pid'] >= opts.pid:
            if psl_dict['qName'] not in query_to_target_dict:
                query_to_target_dict[psl_dict['qName']] = psl_dict 
            else:
                if psl_dict['target_coverage'] > query_to_target_dict[psl_dict['qName']]['target_coverage']:
                    query_to_target_dict[psl_dict['qName']] = psl_dict
           
    return query_to_target_dict


def MakeBlockIntervals(psldict):
    block_starts = [int(i) for i in psldict['tStarts'].split(',')[:-1]]
    block_sizes = [int(i) for i in psldict['blockSizes'].split(',')[:-1]]
    block_intervals=[]
    for i in range(len(block_starts)):
        block_intervals.append([block_starts[i],block_starts[i]+block_sizes[i]])
    return block_intervals

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='pipeline for generating sundry transcriptome assembly coverage statistics')
    parser.add_argument('-psl','--psl_infile',dest='hits',type=str,help='concatenated psl file bits from job array, w/o headers')
    parser.add_argument('-o','--outfile_prefix',dest='outprefix',type=str,help='prefix for output file')
    parser.add_argument('-minc','--min_query_coverage',dest='mincov',default=0,type=float,help='proportion of the query covered by matches in the alignment') 
    parser.add_argument('-pid','--percent_identity',dest='pid',default=0.95,type=float,help='min % identity in the alignment filtering threshold')
    opts = parser.parse_args() 

    psl_keys=['matches','misMatches','repMatches','nCount','qNumInsert','qBaseInsert','tNumInsert','tBaseInsert','strand','qName','qSize','qStart','qEnd','tName','tSize','tStart','tEnd','blockCount','blockSizes','qStarts','tStarts']
    
    pslin = open(opts.hits,'r')
    besthit_per_query = BuildBestHitPerQueryDict(pslin)
    
    bedout=open('%s_besthitperquery.bed' % opts.outprefix,'w')
    for key in besthit_per_query:
        intervals = MakeBlockIntervals(besthit_per_query[key])
        for interval in intervals:
            bedout.write('%s\t%s\t%s\t%s\n' % (besthit_per_query[key]['tName'],interval[0],interval[1],key))
        
    bedout.close()


