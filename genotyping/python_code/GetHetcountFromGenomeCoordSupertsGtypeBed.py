import sys
hetcount = 0
fopen=open(sys.argv[1],'r')
for line in fopen:
    if 'PASS' in line:
        linelist = line.strip().split()
        refallele=linelist[9]
        altalleles=linelist[10]
        if len(refallele) ==1 and len(altalleles)==1:
            print refallele,altalleles
            alleles=linelist[15].split(':')[0].split('/')
            #print alleles
            unique = []
            for allele in alleles:
                if allele not in unique:
                    unique.append(allele)
            if len(unique)==2:
                hetcount+=1
                
print 'hetcount is ...%s' % hetcount
