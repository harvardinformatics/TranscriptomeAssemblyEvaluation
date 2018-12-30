from sets import Set
from collections import defaultdict
import sys
min200ts = open('/n/holylfs/LABS/informatics/adamf/refs/mus/GRCm38/rsem_ensembl_indices/mus_GRCm38.83_minlen200_tsids.txt','r')

pcod_ts = open('/n/holylfs/LABS/informatics/adamf/refs/mus/GRCm38/mus_proteincoding_ts.ids','r')
pcod_ts_list = []
for line in pcod_ts:
   pcod_ts_list.append(line.strip())
   
pcod_gene = open('/n/holylfs/LABS/informatics/adamf/refs/mus/GRCm38/mus_proteincoding_gene.ids','r')   
pcod_genes = []
for line in pcod_gene:
    pcod_genes.append(line.strip())

min200ts_list = []
for line in min200ts:
    min200ts_list.append(line.strip())
    
pcod_ts_set = Set(pcod_ts_list)
min200_ts_set = Set(min200ts_list)

passfilter_ts = pcod_ts_set.intersection(min200_ts_set)

#####################
gene_ts_map = open('/n/holylfs/LABS/informatics/adamf/refs/mus/GRCm38/gene_ts_map_2018.11.28.tsv','r')

gtmap_dict = defaultdict(list)
for line in gene_ts_map:
    ts,gene = line.strip().split()
    gtmap_dict[gene].append(ts)
    
pf_genes = []
pfts_nogene = 0
for gene in gtmap_dict.keys():
    if len(Set(gtmap_dict[gene]).intersection(passfilter_ts)) > 0:
        if gene not in pcod_genes:
            #print 'wtf!',gene,Set(gtmap_dict[gene]).intersection(passfilter_ts)
            pfts_nogene +=1
        else:
            pf_genes.append(gene)

####################    
    
iso_matrix = open(sys.argv[1],'r')
iso_header = iso_matrix.readline()
expr_iso_matrix=open('expressed_%s' % sys.argv[1],'w')
filt_iso_matrix = open('pcod_expressed_minlen200_%s' % sys.argv[1],'w')
filt_iso_matrix.write(iso_header)
expr_iso_matrix.write(iso_header)

gene_matrix = open(sys.argv[2],'r')
expr_gene_matrix=open('expressed_%s' % sys.argv[2],'w')
filt_gene_matrix = open('pcod_expressed_minlen200_%s' % sys.argv[2],'w')
gene_header = gene_matrix.readline()
filt_gene_matrix.write(gene_header)
expr_gene_matrix.write(gene_header)

for line in iso_matrix:
    iso,tpm = line.strip().split()
    if float(tpm) > 0:
        expr_iso_matrix.write(line)
        if iso in passfilter_ts:
            filt_iso_matrix.write(line)    

for line in gene_matrix:
    gene,tpm = line.strip().split()
    if float(tpm) > 0:
        expr_gene_matrix.write(line)
        if gene in pf_genes:
            filt_gene_matrix.write(line)
        
filt_iso_matrix.close()
filt_gene_matrix.close() 
expr_iso_matrix.close()
expr_gene_matrix.close()

print 'pf ts but no pcod gene count = %s' % pfts_nogene
