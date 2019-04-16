import os
import sys
from sets import Set
#fields  = ['transcript_id','gene_id','length','effective_length','expected_count','TPM','FPKM','IsoPct']
iso_mapref=open(sys.argv[1],'r')
header =  iso_mapref.readline()
fields = header.strip().split('\t')
mapref_dict={}
isoforms  = Set()
isoform_lengths = {}
for line in iso_mapref:
    linelist =  line.strip().split('\t')
    linedict = dict(zip(fields,linelist))
    isoforms.add(linelist[0])
    mapref_dict[linelist[0]] = linedict

denovo_dict = {}
iso_denovo = open(sys.argv[2],'r')
iso_denovo.readline()
for line in iso_denovo:
    linelist = line.strip().split('\t')
    linedict = dict(zip(fields,linelist))
    if 'ENS' in linedict['gene_id']:
        isoforms.add(linedict['transcript_id'])
        denovo_dict[linelist[0]] = linedict

mr_rev = open('overlap_revised_%s' % sys.argv[1],'w')
mr_rev.write(header)

denovo_rev = open('overlap_revised_%s' % sys.argv[2],'w')
denovo_rev.write(header)

tximportmap=open('tximport_map','w')
tximportmap.write('TXNAME\tGENEID\n')
for isoform in isoforms:
    if isoform in mapref_dict:
        tximportmap.write('%s\t%s\n' % (isoform,mapref_dict[isoform]['gene_id']))
    else:
        tximportmap.write('%s\t%s\n' % (isoform,denovo_dict[isoform]['gene_id']))

tximportmap.close()

rscript = open('tximport.Rscript','w')
rscript.write('setwd("%s")\n' % os.getcwd())
rscript.write('library(tximport)\n')
rscript.write('myfiles = c("overlap_revised_%s","overlap_revised_%s")\n' % (sys.argv[1],sys.argv[2]))
rscript.write('tx2gene = read.table("tximport_map",header=TRUE)\n')
rscript.write('txi.rsem<-tximport(myfiles,type="rsem",txIn = TRUE,tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")\n')
rscript.write('txi.dataframe<-as.data.frame(txi.rsem$counts)\n')
rscript.write('txi.dataframe$geneid<-row.names(txi.rsem$counts)\n')
rscript.write('colnames(txi.dataframe)<-c("mapref","denovo","gene_id")\n')
rscript.write('txi.dataframe<-subset(txi.dataframe,txi.dataframe$mapref > 0 | txi.dataframe$denovo > 0)\n') 
rscript.write('write.table(txi.dataframe,"%s_lstpmMRvsDeNovo.tsv",row.names=FALSE,quote=FALSE)\n' % sys.argv[3])
rscript.close() 


for isoform in isoforms:
    if isoform in mapref_dict:
        mrline =[mapref_dict[isoform][field] for field in fields]
        mr_rev.write('%s\n' % '\t'.join(mrline))
        denovo_rev.write('%s\t0\t0\t0\t0\t0\n' % '\t'.join(mrline[0:3]))
    else:
        denovoline =[denovo_dict[isoform][field] for field in fields]
        mr_rev.write('%s\t0\t0\t0\t0\t0\n' % '\t'.join(denovoline[0:3]))
        denovo_rev.write('%s\n' % '\t'.join(denovoline))
        
denovo_rev.close()
mr_rev.close()
