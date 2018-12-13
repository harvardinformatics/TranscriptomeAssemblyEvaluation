import sys
isoform_list = []
expr_iso_handle = open(sys.argv[1],'r')
for line in expr_iso_handle:
    isoform_list.append(line.strip())
#print isoform_list 

gtfin = open(sys.argv[2],'r')
gtfout = open('expressed_%s' % sys.argv[2],'w')
for line in gtfin:
    if line[0] == '#':
        gtfout.write(line)
        
    elif 'ENSMUST' in line:
        ts=line.strip().split('\t')[8].split('transcript_id')[1].split(';')[0].replace('"','')
        ts = ts.replace(' ','')
        if ts in isoform_list:
            gtfout.write(line)
            #print line
gtfout.close()
