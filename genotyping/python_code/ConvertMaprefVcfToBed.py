import argparse
from ConvertSuperTranscriptDataToGenomicCoordinates import IntersectWithExons

def RefVcfToBed(vcf,patchfilter=True):
    vcfin = open(vcf,'r')
    bedout = open('%sbed' % vcf[:-3],'w')
    if patchfilter == True:
        for line in vcfin:
            if line[0] != '#':
                if 'CHR' not in line:
                    linelist = line.strip().split('\t')
                    bedout.write('%s\t%s\t%s\t%s\n' % (linelist[0],int(linelist[1])-1,linelist[1],'\t'.join(linelist[3:])))
    else:
        for line in vcfin:
            if line[0] != '#': 
                linelist = line.strip().split('\t')
                bedout.write('%s\t%s\t%s\t%s\n' % (linelist[0],int(linelist[1])-1,linelist[1],'\t'.join(linelist[3:])))
    bedout.close()

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='Converts map-to-ref vcf file from RNA-seq data to exonic genotypes in bed format')
    parser.add_argument('-rvcf','--ref_vcf',dest='refvcf',type=str,help='vcf of map-to-reference genotypes')
    parser.add_argument('-e','--exon_bed',dest='exons',type=str,help='bed file of genomic exons')
    opts = parser.parse_args()

    RefVcfToBed(opts.refvcf)
    IntersectWithExons('%sbed' % opts.refvcf[:-3],opts.exons)
