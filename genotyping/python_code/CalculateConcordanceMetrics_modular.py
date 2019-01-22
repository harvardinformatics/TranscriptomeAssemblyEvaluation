import argparse
from sets import Set

fields = ['snpid', 'genomicpositions', 'supertspositions', 'mapref_ref', 'maprefalleles', 'superts_ref', 'supertsalleles', 'maprefcov', 'mapref_ad', 'superts_ad','mean_supercov','median_supercov']

def BoolFlagGenotypeSetConcordance(snp_dict):
    """
    assesses whether allele sets are the same
    for mapref and SuperTranscript based genotyping methods
    for each snp cluster id
    """
    mapref_allele_set = Set(snp_dict['maprefalleles'].split(';'))
    superts_allele_set = Set(snp_dict['supertsalleles'].split(';'))
    if mapref_allele_set == superts_allele_set:
       return True
    else:
        return False

def Recall(snp_dict):
    """
    determine if a called map-ref genotype
    is recovered by SuperTranscript approach
    """
    if snp_dict['maprefalleles'] != 'NA':
        mapref_allele_set = Set(snp_dict['maprefalleles'].split(';'))
        superts_allele_set = Set(snp_dict['supertsalleles'].split(';'))
        if mapref_allele_set == superts_allele_set:
            return True
        else:
            return False

    else:
        return 'NA'

def HetRecall(snp_dict):
    """
    determine if a called map-ref het genotype
    is recovered by SuperTranscript approach
    """
    if snp_dict['maprefalleles'] != 'NA':
        mapref_allele_set = Set(snp_dict['maprefalleles'].split(';'))
        superts_allele_set = Set(snp_dict['supertsalleles'].split(';'))
        if len(mapref_allele_set) == 2 and mapref_allele_set == superts_allele_set:
            return True
        else:
            return False

    else:
        return 'NA'

def FalsePositive(snp_dict):
    if snp_dict['maprefalleles'] == 'NA' and snp_dict['supertsalleles'] != 'NA':
        return True
    else:
        return False

def FalsePositiveHet(snp_dict):
    """
    cases where no genotype, or a fixed alternative site
    in map-to-ref is called as a bi-allelic het via
    SuperTranscript method
    """
    if snp_dict['maprefalleles'] == 'NA' or len(snp_dict['maprefalleles']) == 1:
        if snp_dict['supertsalleles'] != 'NA' and len(snp_dict['supertsalleles'].split(';')) == 2:
            return True
        else:
            return False
    else:
        return 'NA'

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='Converts map-to-ref vcf file from RNA-seq data to exonic genotypes in bed format')
    parser.add_argument('-gt','--gtable',dest='genotypes',type=str,help='table of mapref and superts genotype intersections')
    opts = parser.parse_args()

    qc_dict = {'fp_het':{'fp_het':0,'counted':0},'fp':{'fp':0,'counted':0},'het_recall':{'pass':0,'counted':0},'recall':{'pass':0,'counted':0},'concordance':{'pass':0,'counted':0}}

    fopen = open(opts.genotypes,'r')
    header_fields = fopen.readline().strip().split()
    if header_fields != fields:
       raise ValueError('table fields do not match required input fields')

    for line in fopen:
        line_list = line.strip().split('\t')
        snp_dict =  dict(zip(fields,line_list)) 
        ### concordance ###
        concordant = BoolFlagGenotypeSetConcordance(snp_dict)
        if concordant == True:
            qc_dict['concordance']['pass'] +=1
            qc_dict['concordance']['counted']+=1
        else:
            qc_dict['concordance']['counted']+=1
        
        ### recall ###
        recall = Recall(snp_dict)
        if recall == True:
            qc_dict['recall']['pass']+=1
            qc_dict['recall']['counted']+=1
        elif recall == False:
            qc_dict['recall']['counted']+=1
        else:
            pass

        ### het recall ###
        hetrecall = HetRecall(snp_dict)
        if hetrecall == True:
            qc_dict['het_recall']['pass']+=1
            qc_dict['het_recall']['counted']+=1
        elif hetrecall == False:
            qc_dict['het_recall']['counted']+=1
        else:
            pass

        ### false positive ###
        fp = FalsePositive(snp_dict)
        if fp == True:
            qc_dict['fp']['fp']+=1
            qc_dict['fp']['counted']+=1
        else:
            qc_dict['fp']['counted']+=1

        ### fp hets ###
        fphet = FalsePositiveHet(snp_dict)
        if fphet == True:
            qc_dict['fp_het']['fp_het']+=1
            qc_dict['fp_het']['counted']+=1
        elif fphet == False:
            qc_dict['fp_het']['counted']+=1
        else:
            pass 



    print qc_dict
