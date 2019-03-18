import argparse

psl_keys=['matches','misMatches','repMatches','nCount','qNumInsert','qBaseInsert','tNumInsert','tBaseInsert','strand','qName','qSize','qStart','qEnd','tName','tSize','tStart','tEnd','blockCount','blockSizes','qStarts','tStarts']


def CreateAlignmentSets(blat_handle,minlength,alignment_proportion):
    matches_dict={}
    linecounter = 0
    for line in blat_handle:
        linecounter+=1
        if linecounter % 10000 == 0:
            print 'processing line %s' % linecounter
        linelist = line.strip().split("\t")
        line_dict = dict(zip(psl_keys,linelist))
        if line_dict['qName'] != line_dict['tName']:
            if int(line_dict['qSize']) >= minlength and int(line_dict['tSize']) >= minlength:
                alignment_length = 0
                for block in line_dict['blockSizes'].split(',')[:-1]:
                    alignment_length+=int(block)
                if alignment_length/float(line_dict['qSize']) >= alignment_proportion or alignment_length/float(line_dict['tSize']) >= alignment_proportion:
                    pid=str(float(line_dict['matches'])/(float(line_dict['matches'])+float(line_dict['repMatches'])+float(line_dict['misMatches'])))
                    matches_dict[(line_dict['qName'],line_dict['tName'])] = [pid,line_dict['qSize'],line_dict['tSize'],alignment_length]

    return matches_dict     


if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='script characterizing frequency of contig sequence similarlity in tscriptome assembly')
    parser.add_argument('-psl','--psl_infile',dest='hits',type=str,help='concatenated psl file bits from blat all-by-all job array, w/o headers')
    parser.add_argument('-o','--outfile',dest='outfile',type=str,help='prefix for output files summarizing alignment overlap stats')
    parser.add_argument('-p','--alignment_proportion',dest='alignment_prop',type=float,help='min proportion of contig constituted by alignment')
    parser.add_argument('-minsize','--minimum_contig_size',dest='minsize',type=int,help='min length of the aligned contigs')
    parser.add_argument('-T','--Trinity',dest='trinity',action='store_true',help='are contigs in Trinity data structures?')
    opts = parser.parse_args() 

    psl_in = open(opts.hits,'r')
    pid_dict = CreateAlignmentSets(psl_in,opts.minsize,opts.alignment_prop)
    fout=open(opts.outfile,'w')
    if opts.trinity == True:
        fout.write('query\ttarget\tpid\tqSize\ttSize\talign_length\ttrinity_comp\n')
        for alignment in pid_dict.keys():
            if '_'.join(alignment[0].split('_')[:-1]) == '_'.join(alignment[1].split('_')[:-1]):
                fout.write('%s\t%s\t%s\t%s\n' % (alignment[0],alignment[1],'\t'.join([str(i) for i in pid_dict[alignment]]),'samecomp'))
            else:         
                fout.write('%s\t%s\t%s\t%s\n' % (alignment[0],alignment[1],'\t'.join([str(i) for i in pid_dict[alignment]]),'diffcomp'))
    else:
        fout.write('query\ttarget\tpid\tqSize\ttSize\talign_length\n')
        for alignment in pid_dict.keys():
            fout.write('%s\t%s\t%s\n' % (alignment[0],alignment[1],'\t'.join([str(i) for i in pid_dict[alignment]])))

    fout.close()
