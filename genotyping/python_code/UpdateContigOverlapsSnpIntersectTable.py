import sys
from sets import Set
from collections import defaultdict
intersect_open=open(sys.argv[1],'r')
from Bio import SeqIO
intersect_header= intersect_open.readline().strip()
intersect_fields = intersect_header.split('\t')
new_header='%s\tSToverlaps\tSToverlapalleles\n' % intersect_header
superts_cov_open = open(sys.argv[2],'r')
##/n/holylfs/LABS/informatics/adamf/DeNovoTranscriptomeEvaluation/genotyping/supertranscript_trinity_clusters/MDC/sites_count/exon_filter/mdc_superts_exons_coverage.bed
rev_comp_dict = {'A':'T','C':'G','G':'C','T':'A'}
st_fasta_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[3], "fasta")) # superts fasta

intersect_update=open('wcontigoverlaps_allelesfromoverlaps_%s' % sys.argv[1],'w')
intersect_update.write(new_header)

cov_fields = ['supertsid','ststart','stend','chrom','chromstart','chromend','strand','stiddup','ststartdup','stenddup','st_cov']

genome_to_superts = defaultdict(list)

for line in superts_cov_open:
    linelist = line.strip().split('\t')
    linedict = dict(zip(cov_fields,linelist))
    pos = '%s:%s' % (linedict['chrom'],linedict['chromend'])
    genome_to_superts[pos].append('%s:%s:%s' % (linedict['supertsid'],linedict['stend'],linedict['strand']))
    
for line in intersect_open:
    linelist =  line.strip().split('\t')
    linedict = dict(zip(intersect_fields,linelist))
    if len(genome_to_superts[linedict['genomicpositions']]) > 0:
        alleles_from_st_refs = Set()
        for st_pos in genome_to_superts[linedict['genomicpositions']]:
            st,stcoord,strand =  st_pos.split(':')
            if strand == '+':
                alleles_from_st_refs.add(st_fasta_dict[st].seq[int(stcoord)-1])
            else:
                alleles_from_st_refs.add(rev_comp_dict[st_fasta_dict[st].seq[int(stcoord)-1]]) 
        intersect_update.write('%s\t%s\t%s\n' % (line.strip(),';'.join(genome_to_superts[linedict['genomicpositions']]),';'.join(alleles_from_st_refs)))
    else:
        intersect_update.write('%s\tNA\tNA\n' % line.strip())
intersect_update.close()
        
