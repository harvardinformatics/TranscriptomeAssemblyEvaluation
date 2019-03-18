import argparse
from collections import defaultdict
from numpy import mean,std

def CalculateCompositionalVariance(contigs,intron_dict,cds_dict):
    intron_bp = []
    cds_bp = []
    for contig in contigs:
        intron_bp.append(int(intron_dict[contig]))
        cds_bp.append(int(cds_dict[contig]))
    
    intron_std = std(intron_bp)
    cds_std = std(cds_bp)

    #intron_diffs = []
    #cds_diffs = []
    #for count in intron_bp:
        #intron_diffs.append(abs(count-mean(intron_bp)))
    #for count in cds_bp:
        #cds_diffs.append(abs(count-mean(cds_bp)))

    #return mean(intron_diffs),mean(cds_diffs) 
    return intron_std,cds_std

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='calc intron/cds comp stats for contigs overlapping single isoform genes')
    parser.add_argument('-bpqbed','--single-isoform-best-per-query-bed',dest='bpq',type=str,help='single isoform best per query bed')
    parser.add_argument('-c','--comp-table',dest='comp',type=str,help='composition table')
    parser.add_argument('-out','--output-table',dest='fout',type=str,help='name output table')
    parser.add_argument('-overlaps','--contig-overlap-table',dest='overlaps',type=str,help='table of contig overlap bases')
    parser.add_argument('-minoverlap','--minimum-overlap-fraction',dest='pover',type=float,help='minimum proportion mapped bases overlapping among contigs')
    opts = parser.parse_args()
    
    overlap_open=open(opts.overlaps,'r')
    overfields = overlap_open.readline().strip().split(',')
    
    filtered_isoforms = []
    for line in overlap_open:
        linelist = line.strip().split(',')
        over_dict = dict(zip(overfields,linelist))
        if int(over_dict['ContigOverlaps'])/float(over_dict['BasesMatch']) >= opts.pover:
            filtered_isoforms.append(over_dict['EnsTs'])


    iso_to_contigs_dict = defaultdict(list)
    bpq_open=open(opts.bpq,'r')
    for line in bpq_open:
        linelist = line.strip().split()
        iso_to_contigs_dict[linelist[0]].append(linelist[-1]) 

    intron_dict = {}           
    cds_dict = {}
    
    comp_open = open(opts.comp,'r')
    fields = comp_open.readline().split(',')
    for line in comp_open:
        linelist = line.strip().split(',')
        line_dict = dict(zip(fields,linelist))
        intron_dict[line_dict['read']] = line_dict['intron']
        cds_dict[line_dict['read']] = line_dict['cds'] 

    
    fout = open(opts.fout,'w')
    fout.write('isoform,introndiff,cdsdiff\n')
    for isoform in iso_to_contigs_dict:
        if isoform in filtered_isoforms:
            if len(iso_to_contigs_dict[isoform]) > 1:
                try:
                    introndiff,cdsdiff =  CalculateCompositionalVariance(iso_to_contigs_dict[isoform],intron_dict,cds_dict)
                    fout.write('%s,%s,%s\n' % (isoform,introndiff,cdsdiff))
                except:
                    fout.write('%s,NA,NA\n' % isoform)

    fout.close()
    
    
