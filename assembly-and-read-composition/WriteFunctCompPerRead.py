from collections import defaultdict
utr=open('utr_contigtallies.txt','r')
sstop=open('startstop_contigtallies.txt','r')
ncintron=open('noncoding_introns_contigtallies.txt','r')
nc=open('noncoding_contigtallies.txt','r')
intron=open('intron_contigtallies.txt','r')
intrainter=open('intra_intergenic_contigtallies.txt','r')
interext=open('intergenic_external_contigtallies.txt','r')
cds=open('cds_contigtallies.txt','r')

mappedbases=open('mappedbasecount_per_contig.txt','r')

compdict=dict()

def parse_function_file(file,compdict):
    header=file.readline()
    funcclass=header.strip().split()[1]
    print funcclass
    for line in file:
        linelist=line.strip().split()
        if linelist[0].split('/')[0] not in compdict:
            compdict[linelist[0].split('/')[0]]=defaultdict(int)
            compdict[linelist[0].split('/')[0]][funcclass]+=int(linelist[1])
        else:
             compdict[linelist[0].split('/')[0]][funcclass]+=int(linelist[1])

intrainter.readline()
for line in intrainter:
    linelist=line.strip().split()
    if linelist[0].split('/')[0] not in compdict:
        compdict[linelist[0].split('/')[0]]=defaultdict(int)
        compdict[linelist[0].split('/')[0]]['intergenic']+=int(linelist[1])
    else:
        compdict[linelist[0].split('/')[0]]['intergenic']+=int(linelist[1])
interext.readline()
for line in interext:
    linelist=line.strip().split()
    if linelist[0].split('/')[0] not in compdict:
        compdict[linelist[0].split('/')[0]]=defaultdict(int)
        compdict[linelist[0].split('/')[0]]['intergenic']+=int(linelist[1])  
    else:
        compdict[linelist[0].split('/')[0]]['intergenic']+=int(linelist[1])  


for handle in [utr,sstop,ncintron,nc,intron,cds]:
    parse_function_file(handle,compdict)


mappedbases.readline()
for line in mappedbases:
    linelist=line.strip().split()
    compdict[linelist[0].split('/')[0]]['mappedbases']+=int(linelist[1]) 

fout=open('basecomp_per_read.csv','w')
fout.write('read,mappedbases,cds,utr,sstop,intron,noncoding,noncoding_intron,intergenic\n')

for key in compdict.keys():
    fout.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (key,compdict[key]['mappedbases'],compdict[key]['cds'],compdict[key]['utr'],compdict[key]['startstop'],compdict[key]['intron'],compdict[key]['noncoding'],compdict[key]['noncoding_introns'],compdict[key]['intergenic']))

fout.close()

