import argparse
from sets import Set

def evaluate_heterozygosity(genotypestring):
    hetset=Set(['0','1'])
    alleles = Set()
    for allele in genotypestring.split(':')[0].split('/'):
        alleles.add(allele)
    if alleles == hetset:
        return 1
    else:
        return 0

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Collapse genotypes to diploid and get counts of heterozygous sites from exon mapref gtype bed obtained from gatk vcf')
    parser.add_argument('-b','--exon-gtype-bed',dest='bedin',type=str,help='input bed file of exon region genotypes for map-to-ref')
    parser.add_argument('-p','--passfilter',action='store_true',help='only consider genotypes with gatk PASS flag')
    opts = parser.parse_args()

    hetcount = 0
    bed_open = open(opts.bedin,'r')
    for line in bed_open:
        linelist = line.strip().split()
        if len(linelist[4]) == 1 and len(linelist[3]) == 1:
            if opts.passfilter == True and 'PASS' in line:
                hetcount += evaluate_heterozygosity(linelist[9])
            elif opts.passfilter == False:
                print 'warning'
                hetcount += evaluate_heterozygosity(linelist[9])
            else:
                pass 

    print 'recorded het sites: %s\n' % hetcount   
