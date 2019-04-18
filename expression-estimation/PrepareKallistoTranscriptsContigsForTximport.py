import os
import sys
from sets import Set
from collections import defaultdict
fields = ['target_id','length','eff_length','est_counts','tpm']
genes = Set()
mapref_quant = open(sys.argv[1],'r')
mapref_quant_dict = {}
header  = mapref_quant.readline()
mapref_fields = header.strip().split()
for line in mapref_quant:
    linelist = line.strip().split()
    linedict =  dict(zip(mapref_fields,linelist))
    mapref_quant_dict[linedict['target_id']] = linedict

denovo_quant = open(sys.argv[2],'r')
denovo_quant_dict = {}
denovo_fields = denovo_quant.readline().strip().split()
for line in denovo_quant:
    linelist = line.strip().split()
    linedict =  dict(zip(denovo_fields,linelist))
    denovo_quant_dict[linedict['target_id']] = linedict 


mapref_map = open(sys.argv[3],'r')
mapref_map_dict = defaultdict(list)
for line in mapref_map:
    gene,ts = line.strip().split()
    mapref_map_dict[gene].append(ts)
    genes.add(gene)

denovo_map = open(sys.argv[4],'r')
denovo_map_dict = defaultdict(list)
for line in denovo_map:
    gene,ts =  line.strip().split()
    if 'ENS' in gene:
        denovo_map_dict[gene].append(ts)
        genes.add(gene)


mapref_rev = open('overlap_revised_%s' % os.path.basename(sys.argv[1]),'w')
mapref_rev.write(header)
denovo_rev = open('overlap_revised_%s' % os.path.basename(sys.argv[2]),'w')
denovo_rev.write(header)

tximportmap=open('tximport_map','w')
tximportmap.write('TXNAME\tGENEID\n')

for gene in genes:
    if gene in mapref_map_dict and gene in denovo_map_dict:
        for ts in mapref_map_dict[gene]:
            ts_data = [mapref_quant_dict[ts][field] for field in fields]
            mapref_rev.write('%s\n' % '\t'.join(ts_data))
            tximportmap.write('%s\t%s\n' % (ts,gene)) 
            denovo_rev.write('%s\t%s\t0\t0\t0\n' % (ts,mapref_quant_dict[ts]['length']))
        for ts in denovo_map_dict[gene]:
            ts_data = [denovo_quant_dict[ts][field] for field in fields]
            denovo_rev.write('%s\n' % '\t'.join(ts_data))
            tximportmap.write('%s\t%s\n' % (ts,gene))
            mapref_rev.write('%s\t%s\t0\t0\t0\n' % (ts,denovo_quant_dict[ts]['length']))

tximportmap.close()
denovo_rev.close()
mapref_rev.close()

## write r script ##
rscript = open('tximport.Rscript','w')
rscript.write('setwd("%s")\n' % os.getcwd())
rscript.write('library(tximport)\n')
rscript.write('library(readr)\n')
rscript.write('myfiles = c("overlap_revised_%s","overlap_revised_%s")\n' % (os.path.basename(sys.argv[1]),os.path.basename(sys.argv[2])))
rscript.write('tx2gene = read.table("tximport_map",header=TRUE)\n')
# counts
rscript.write('txi.kallisto.counts<-tximport(myfiles,type="kallisto",txIn = TRUE,tx2gene = tx2gene)\n')
rscript.write('txi.counts.dataframe<-as.data.frame(txi.kallisto.counts$counts)\n')
rscript.write('txi.counts.dataframe$geneid<-row.names(txi.counts.dataframe)\n')
rscript.write('colnames(txi.counts.dataframe)<-c("mapref_count","denovo_count","gene_id")\n')
rscript.write('txi.length.dataframe<-as.data.frame(txi.kallisto.counts$length)\n')
rscript.write('txi.length.dataframe$geneid<-row.names(txi.length.dataframe)\n')
rscript.write('colnames(txi.length.dataframe)<-c("mapref_efflen","denovo_efflen")\n')
rscript.write('colnames(txi.length.dataframe)<-c("mapref_efflen","denovo_efflen","gene_id")\n')
rscript.write('countmerged<-merge(txi.counts.dataframe,txi.length.dataframe,by=c("gene_id"))\n')
# tpm
rscript.write('txi.tpm.dataframe<-as.data.frame(txi.kallisto.counts$abundance)\n')
rscript.write('txi.tpm.dataframe$geneid<-row.names(txi.tpm.dataframe)\n')
rscript.write('colnames(txi.tpm.dataframe)<-c("mapref_tpm","denovo_tpm","gene_id")\n')
rscript.write('countmerged<-merge(countmerged,txi.tpm.dataframe,by=c("gene_id"))\n')
rscript.write('countmerged<-subset(countmerged,countmerged$mapref_count > 0 | countmerged$denovo_count > 0)\n')
rscript.write('write.csv(countmerged,"%s_kallistoMRvsDeNovo.csv",row.names=FALSE,quote=FALSE)\n' % sys.argv[5])
rscript.close() 
