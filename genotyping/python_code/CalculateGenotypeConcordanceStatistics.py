import argparse
"""
QC metric definitions:
1. Recall (sensitivity, true positive rate) == proportion of all called mapref genotypes that are called correctly with supertranscripts
2. Recall_hets == proportion of all mapref het genotypes that are called correctly
3. FP == frequency of sites that have superts genotypes but no mapref genotype
4. FP_hets == frequency of sites with a het superts genotype that are without a mapref genotype or with a homozygous alt mapref call
5. FN == frequency of sites with a mapref genotype call that have no superts genotype call
6. FN_hets == frequency of mapref het calls that don't have a corresponding superts het call
7. error_rate = frequency of sites that have conflicting mapref and superts calls, including
   no call in mapref but called in superts, or called in mapref but no call in superts. 
8. specificity == proportion of supertranscripts genotypes that are observed correctly relative to mapref, includes superts 'NA'
9. snv_to_indel_error == frequency of called mapref SNV genotypes that are called as indel polymorphisms in superts
10. allele_included == frequency of cases where an erroneous superts genotype call has, as one of its alleles, an allele in the mapref genotype call
11. multi_perror == frequency of genotyping errors where there is a mapref genotype but multiple supertranscripts aligned to the position
1st and second element of each list are for the numerator and denominator used to calculate each statistic
12. diffhets == proportion of cases where map ref calls one het, and superts calls another
"""

def CalculateMetrics(filehandle):
    fields = ['snpid', 'genomicpositions', 'supertspositions', 'mapref_ref', 'maprefalleles', 'superts_ref', 'supertsalleles', 'maprefcov', 'mapref_ad', 'superts_ad']
    """
    field defs:
        snpid == unique numbered id; prefix is 'genomic' if not found only in map-to-ref, i.e no superts genptype
        genomicpositions == semi-colon separated list of all genomic positions in the snp cluster
        supertspositions == semi-colon separated list or all superts positions in the snp cluster
        mapref_ref == map-to-reference reference allele
        maprefallele == semi-colon separated list of all observed alleles
        superts_ref == supertranscript reference allele, reverse-complemented if superts maps to - strand of genome
        supertsalleles == semi-colon separated list of observed superts alleles, rev comp as above if on - strand
        maprefcov == coverage depth in reference genome, whether or not a mapref genotype is called
        mapref_ad == depth by allele for mapref gatk genotypes, semi-colon separated
        superts_ad == depth by allele for superts gatk genotypes, semi-colon separated
    """
 
    metric_dict={
        'recall': [0,0],
        'recall_hets': [0,0],
        'fp': [0,0],
        'fp_hets': [0,0],
        'fn': [0,0],
        'fn_hets': [0,0],
        'error_rate': [0,0],
        'specificity': [0,0],
        'snv2indel': [0,0],
        'allele_incl': [0,0],
        'multi_perror': [0,0],
        'diff_hets': [0,0],
    }


    for line in filehandle:
        linedict = dict(zip(fields,line.strip().split('\t')))
        genomic_positions = len(linedict['genomicpositions'].split(';'))
        superts_positions = len(linedict['supertspositions'].split(';')) 
        maprefalleles = linedict['maprefalleles'].split(';')
        maprefalleles.sort()
        supertsalleles = linedict['supertsalleles'].split(';')
        supertsalleles.sort()

        #### concordant ####
        if maprefalleles == supertsalleles:
            metric_dict['recall'][0]+=1 ; metric_dict['recall'][1]+=1
            metric_dict['specificity'][0]+=1 ; metric_dict['specificity'][1]+=1
            metric_dict['error_rate'][1]+=1
            metric_dict['fp'][1]+=1
            metric_dict['fn'][1]+=1

            if len(maprefalleles) > 1:
                metric_dict['recall_hets'][0]+=1 ; metric_dict['recall_hets'][1]+=1
                metric_dict['fp_hets'][1]+=1
                metric_dict['fn_hets'][1]+=1
                metric_dict['diff_hets'][1]+=1

            mapref_allele_lengths = []
            
            for allele in maprefalleles:
                mapref_allele_lengths.append(len(allele))
            mapref_allele_lengths.append(len(linedict['mapref_ref']))
            if max(mapref_allele_lengths) == 1:
                metric_dict['snv2indel'][1]+=1

        #### missing from map-to-reference #### 
        elif maprefalleles == ['NA']:
            metric_dict['specificity'][1]+=1                 
            metric_dict['fp'][0]+=1 ; metric_dict['fp'][1]+=1
            metric_dict['error_rate'][0]+=1 ; metric_dict['error_rate'][1]+=1
            if superts_positions > 1:
                metric_dict['multi_perror'][0]+=1 ; metric_dict['multi_perror'][1]+=1
            else:
                metric_dict['multi_perror'][1]+=1
            if len(supertsalleles) > 1:
                metric_dict['fp_hets'][0]+=1 ; metric_dict['fp_hets'][1]+=1
        #### missing from superts ####
        elif supertsalleles == ['NA']:
            metric_dict['fn'][0]+=1 ; metric_dict['fn'][1]+=1
            metric_dict['specificity'][1]+=1
            metric_dict['error_rate'][0]+=1 ; metric_dict['error_rate'][1]+=1
            metric_dict['recall'][1]+=1
            if len(maprefalleles) > 1:
                metric_dict['fn_hets'][0]+=1 ; metric_dict['fn_hets'][1]+=1
                metric_dict['recall_hets'][1]+=1

        #### discordance not due to missing genotypes ####
        else:
            metric_dict['error_rate'][0]+=1 ; metric_dict['error_rate'][1]+=1
            metric_dict['specificity'][1]+=1            
            metric_dict['recall'][1]+=1
            if len(maprefalleles) > 1 and len(supertsalleles) > 1:
                metric_dict['diff_hets'][0]+=1 ; metric_dict['diff_hets'][1]+=1
            elif len(maprefalleles) > 1 and len(supertsalleles) == 1:
                metric_dict['fn_hets'][0]+=1 ; metric_dict['fn_hets'][1]+=1
                metric_dict['recall_hets'][1]+=1
            elif len(maprefalleles) == 1 and len(supertsalleles) > 1:
                metric_dict['fp_hets'][0]+=1 ; metric_dict['fp_hets'][1]+=1
         
            included_count = 0
            for allele in supertsalleles:
                if allele in maprefalleles:
                    included_count+=1
            if included_count > 0:
                metric_dict['allele_incl'][0]+=1 ; metric_dict['allele_incl'][1]+=1
            
            if superts_positions > 1:
                metric_dict['multi_perror'][0]+=1 ; metric_dict['multi_perror'][1]+=1
            else:
                metric_dict['multi_perror'][1]+=1

            mapref_allele_lengths = []
            superts_allele_lengths = []

            for allele in maprefalleles:
                mapref_allele_lengths.append(len(allele))
            mapref_allele_lengths.append(len(linedict['mapref_ref']))
            
            for allele in supertsalleles:
                superts_allele_lengths.append(len(allele))
            superts_allele_lengths.append(len(linedict['superts_ref']))
            
            if max(mapref_allele_lengths) == 1 and max(superts_allele_lengths) > 1:
                metric_dict['snv2indel'][0]+=1 ; metric_dict['snv2indel'][1]+=1
    return metric_dict

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='Converts map-to-ref vcf file from RNA-seq data to exonic genotypes in bed format')
    parser.add_argument('-gt','--gtable',dest='genotypes',type=str,help='table of mapref and superts genotype intersections')
    opts = parser.parse_args()


    tablein = open(opts.genotypes,'r')
    tablein.readline()
    tableout = open('concordance_metrics_%s' % opts.genotypes,'w')
    tableout.write('metric\tmetric_count\tn_observations\n')
    concordance_metrics = CalculateMetrics(tablein)
    for metric in concordance_metrics:
        tableout.write('%s\t%s\t%s\n' % (metric,concordance_metrics[metric][0],concordance_metrics[metric][1]))

    tableout.close() 
