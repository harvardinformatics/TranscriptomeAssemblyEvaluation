import argparse
"""
QC metric definitions:
1. Recall == proportion of all called mapref genotypes that are called correctly with supertranscripts
2. Recall_hets == proportion of all mapref het genotypes that are called correctly
3. FP == frequency of sites without a mapref genotype call that have superts genotypes
4. FP_hets == frequency of sites without a mapref genotype or with a homozygous alt call that have het superts genotypes
5. FN == frequency of sites with a mapref genotype call that have no superts genotype call
6. FN_hets == frequency of mapref het calls that don't have a corresponding superts het call
7. error_rate = frequency of sites that have conflicting mapref and superts calls, including
   no call in mapref but called in superts, or called in mapref but no call in superts. 1-error_rate = global_specificity
8. specificity == proportion of supertranscripts genotypes that are observed correctly relative to mapref
9. snv_to_indel_error == frequency of mapref SNV genotypes that are called as indel polymorphisms in superts
10. allele_included == frequency of cases where an erroneous superts genotype call has, as one of its alleles, an allele in the mapref genotype call
11. multi_perror = frequency of genotyping errors where there is a mapref genotype but multiple supertranscripts aligned to the position
"""

def CalculateMetrics(filehandle):
    fields = ['snpid', 'genomicpositions', 'supertspositions', 'mapref_ref', 'maprefalleles', 'superts_ref', 'supertsalleles', 'maprefcov', 'mapref_ad', 'superts_ad']
    metric_dict={
        'recall': [0,0],
        'recall_hets': [0,0],
        'fp': [0,0],
        'fp_hets': [0,0],
        'fn': [0,0],
        'fn_hets': [0,0],
        'error': [0,0],
        'specificity': [0,0],
        'snv2indel': [0,0],
        'allele_incl': [0,0],
        'multi_perror': [0,0]
    }


    for line in filehandle:
        linedict = dict(zip(fields,line.strip().split('\t')
        genomic_positions = len(linedict['genomicpositions'].split(';')
        superts_positions = len(linedict['supertspositions'].split(';') 
        maprefalleles = linedict['maprefalleles'].split(';')
        maprefalleles.sort()
        supertsalleles = linedict['supertsalleles'].split(';')
        supertsalleles.sort()

        if maprefalleles == supertsalleles:
            metric_dict['recall'][0]+=1 ; metric_dict['recall'][1]+=1
            metric_dict['specificity'[[0]+=1 ; metric_dict['specificity'][1]+=1

            if len(maprefalleles) > 1:
                metric_dict['recall_hets'][0]+=1 ; metric_dict['recall_hets'][1]+=1
            elif len(maprefalleles) == 1 and maprefalleles[0] == linedict['mapref_ref']:
                metric_dict['fp_hets'][1]+=1    
            




if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='Converts map-to-ref vcf file from RNA-seq data to exonic genotypes in bed format')
    parser.add_argument('-gt','--gtable',dest='genotypes',type=str,help='table of mapref and superts genotype intersections')
    opts = parser.parse_args()
