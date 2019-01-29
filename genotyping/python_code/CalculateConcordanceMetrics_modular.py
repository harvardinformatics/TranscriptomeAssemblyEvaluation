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

def Specificity(snp_dict):
    """
    proportion of called SuperTranscript
    genotypes that match map-to-ref
    """
    if snp_dict['supertsalleles'] != 'NA':
        mapref_allele_set = Set(snp_dict['maprefalleles'].split(';'))
        superts_allele_set = Set(snp_dict['supertsalleles'].split(';'))
        if mapref_allele_set == superts_allele_set:
            return True
        else:
            return False
    else:
        return 'NA'

def MrefSnvSuperTsIndel(snp_dict):
    """
    frequency of sites where SuperTs gtype error
    due to a map-to-ref SNV called as an 
    indel with SuperTranscript
    """
    if snp_dict['supertsalleles'] != 'NA' and snp_dict['maprefalleles'] != 'NA':
        mapref_allele_lengths = Set([len(allele) for allele in Set(snp_dict['maprefalleles'].split(';'))])
        superts_allele_lengths = Set([len(allele) for allele in Set(snp_dict['supertsalleles'].split(';'))])
        
        if Set(snp_dict['maprefalleles'].split(';')) != Set(snp_dict['supertsalleles'].split(';')):
            if len(mapref_allele_lengths) == 1 and len(superts_allele_lengths) > 1:
                return True
            else:
                return False
    else:
        return 'NA'
        

def RefAlleleIncludedInSuperTsError(snp_dict):
    """
    when gtype called for both map-to-ref
    and SuperTranscripts, when SuperTs gtype error,
    whether a map-to-ref allele is included in 
    the erroneous genotype
    """
    if snp_dict['supertsalleles'] != 'NA' and snp_dict['maprefalleles'] != 'NA':
        mapref_allele_set = Set(snp_dict['maprefalleles'].split(';'))
        superts_allele_set = Set(snp_dict['supertsalleles'].split(';'))
        if mapref_allele_set != superts_allele_set:
            if len(mapref_allele_set.intersection(superts_allele_set)) > 0:
                return True
            else:
                return False
    else:
        return 'NA'    

def FalsePositiveIndel(snp_dict):
    """
    frequency of false positives that
    are called as indels
    """
    if snp_dict['maprefalleles'] == 'NA' or len(snp_dict['maprefalleles'].split(';')) == 1:
        if snp_dict['supertsalleles'] != 'NA':
            superts_allele_lengths = Set([len(allele) for allele in Set(snp_dict['supertsalleles'].split(';'))])
            if len(superts_allele_lengths) > 1:
                return True
            else:
                return False
        else:
            return 'NA'

    else:
        return 'NA'

    mapref_allele_lengths = Set([len(allele) for allele in Set(snp_dict['maprefalleles'].split(';'))])
    superts_allele_lengths = Set([len(allele) for allele in Set(snp_dict['supertsalleles'].split(';'))])

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
    elif snp_dict['maprefalleles'] != 'NA' and len(snp_dict['maprefalleles'].split(';')) == 1 and len(snp_dict['supertsalleles'].split(';')) > 1:
        return True
    else:
        return False

def FalsePositiveHet(snp_dict):
    """
    cases where no genotype, or a fixed alternative site
    in map-to-ref is called as a bi-allelic het via
    SuperTranscript method
    """
    if snp_dict['maprefalleles'] == 'NA' or len(snp_dict['maprefalleles'].split(';')) == 1:
        if snp_dict['supertsalleles'] != 'NA' and len(snp_dict['supertsalleles'].split(';')) == 2:
            return True
        else:
            return False
    else:
        return 'NA'

def FalsePositiveNatoHet(snp_dict):
    if snp_dict['maprefalleles'] == 'NA':
        if snp_dict['supertsalleles'] != 'NA' and len(snp_dict['supertsalleles'].split(';')) == 2:
            return True
        else:
            return False
    else:
        return 'NA'

def FalsePositiveNaToFixedAlt(snp_dict):
    if snp_dict['maprefalleles'] == 'NA':
        if snp_dict['supertsalleles'] != 'NA' and len(snp_dict['supertsalleles'].split(';')) == 1:
            return True
        else:
            return False
    else:
        return 'NA'


def FalseNegative(snp_dict):
    """
    sites where there is a map-to-ref genotype
    call and no genotype call in SuperTranscripts
    """
    if snp_dict['maprefalleles'] != 'NA':
        if snp_dict['supertsalleles'] == 'NA':
            return True
        else:
            return False
    else:
         return 'NA'

def MultiSuperError(snp_dict):
    """
    proportion of discordant genotypes
    where multiple SuperTranscripts align
    to a genomic location
    """
    mapref_allele_set = Set(snp_dict['maprefalleles'].split(';'))
    superts_allele_set = Set(snp_dict['supertsalleles'].split(';'))
    if mapref_allele_set != superts_allele_set:
        if len(snp_dict['supertspositions'].split(';')) > 1:
            return True
        else:
            return False
    else:
        return 'NA'

def FalseNegativeHet(snp_dict):
    """
    cases where a map-ref het is called but
    in SuperTranscript neither a het call
    nor a fixed alternative allele gtype
    """
    if len(snp_dict['maprefalleles'].split(';')) == 2:
        if len(snp_dict['supertsalleles'].split(';')) == 1: # note this includes NAs
            return True
        elif len(snp_dict['supertsalleles'].split(';')) > 2: # considers multi-allelic fns
            return True
        else:
            return False
    else:
        return 'NA'
        

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='Converts map-to-ref vcf file from RNA-seq data to exonic genotypes in bed format')
    parser.add_argument('-gt','--gtable',dest='genotypes',type=str,help='table of mapref and superts genotype intersections')
    parser.add_argument('-sid','--sample_id', dest='sampleid',type=str,help='sample label for table writing')
    opts = parser.parse_args()

    qc_dict = {'fp_natofixedalt':{'fp_natofixedalt':0,'counted':0},'fp_natohet':{'fp_natohet':0,'counted':0},'err_multisuper':{'err_multisuper':0,'counted':0},'fp_indel':{'fp_indel':0,'counted':0},'err_snv2indel':{'err_snv2indel':0,'counted':0},'err_mrefincl':{'err_mrefincl':0,'counted':0},'specificity':{'specificity':0,'counted':0},'fn_het':{'fn_het':0,'counted':0},'fn':{'fn':0,'counted':0},'fp_het':{'fp_het':0,'counted':0},'fp':{'fp':0,'counted':0},'het_recall':{'het_recall':0,'counted':0},'recall':{'recall':0,'counted':0},'concordance':{'concordance':0,'counted':0}}
    print 'num metrics', len(qc_dict.keys())
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
            qc_dict['concordance']['concordance'] +=1
            qc_dict['concordance']['counted']+=1
        else:
            qc_dict['concordance']['counted']+=1
        
        #### superts error with mapref allele included ###
        mrefincl = RefAlleleIncludedInSuperTsError(snp_dict)
        if mrefincl == True:
            qc_dict['err_mrefincl']['err_mrefincl']+=1
            qc_dict['err_mrefincl']['counted']+=1
        elif mrefincl ==False:
            qc_dict['err_mrefincl']['counted']+=1
        else:
            pass

        ### fraction super gtype errors call snvs as indels ###
        snv2indel = MrefSnvSuperTsIndel(snp_dict)
        if snv2indel == True:
            qc_dict['err_snv2indel']['err_snv2indel']+=1
            qc_dict['err_snv2indel']['counted']+=1
        elif snv2indel == False:
            qc_dict['err_snv2indel']['counted']+=1
        else:
            pass

        ### specificity ###
        specificity = Specificity(snp_dict)
        if specificity == True:
            qc_dict['specificity']['specificity'] +=1
            qc_dict['specificity']['counted']+=1
        elif specificity == False:
            qc_dict['specificity']['counted']+=1
        else:
            pass

        ### recall ###
        recall = Recall(snp_dict)
        if recall == True:
            qc_dict['recall']['recall']+=1
            qc_dict['recall']['counted']+=1
        elif recall == False:
            qc_dict['recall']['counted']+=1
        else:
            pass

        ### het recall ###
        hetrecall = HetRecall(snp_dict)
        if hetrecall == True:
            qc_dict['het_recall']['het_recall']+=1
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

        ### false positives as indels ###
        fpindel = FalsePositiveIndel(snp_dict)
        if fpindel == True:
            qc_dict['fp_indel']['fp_indel']+=1
            qc_dict['fp_indel']['counted']+=1
        elif fpindel == False:
            qc_dict['fp_indel']['counted']+=1
        else:
            pass

        ### fp hets ###
        fphet = FalsePositiveHet(snp_dict)
        if fphet == True:
            qc_dict['fp_het']['fp_het']+=1
            qc_dict['fp_het']['counted']+=1
        elif fphet == False:
            qc_dict['fp_het']['counted']+=1
        else:
            pass 
        
        ### fp na to het ###
        fpnatohet = FalsePositiveNatoHet(snp_dict)
        if fpnatohet == True:
            qc_dict['fp_natohet']['fp_natohet']+=1
            qc_dict['fp_natohet']['counted']+=1
        elif fpnatohet == False:
            qc_dict['fp_natohet']['counted']+=1
        else:
            pass

        ### fp na to fixed alt ###
        fpnatofixedalt = FalsePositiveNaToFixedAlt(snp_dict)
        if fpnatofixedalt == True:
            qc_dict['fp_natofixedalt']['fp_natofixedalt']+=1
            qc_dict['fp_natofixedalt']['counted']+=1
        elif fpnatofixedalt == False:
            qc_dict['fp_natofixedalt']['counted']+=1
        else:
            pass

        ### false negative ###
        fn = FalseNegative(snp_dict)
        if fn == True:
            qc_dict['fn']['fn']+=1
            qc_dict['fn']['counted']+=1
        elif fn == False:
            qc_dict['fn']['counted']+=1
        else:
            pass

        ### false negative het ###
        fnhet = FalseNegativeHet(snp_dict)
        if fnhet == True:
            qc_dict['fn_het']['fn_het']+=1
            qc_dict['fn_het']['counted']+=1
        elif fnhet == False:
            qc_dict['fn_het']['counted']+=1
        else:
            pass

        ### multi superts alignment error ###
        errmultsuper = MultiSuperError(snp_dict)
        if errmultsuper == True:
            qc_dict['err_multisuper']['err_multisuper']+=1
            qc_dict['err_multisuper']['counted']+=1
        elif errmultsuper == False:
            qc_dict['err_multisuper']['counted']+=1
        else:
            pass
    print qc_dict

    fout= open('%s.concordancemetrics.tsv' % opts.genotypes[:-4],'w')
    
    qc_keys =  qc_dict.keys()
    qc_keys.sort()
    header='%s\tSample\n' % '\t'.join(qc_keys)
    fout.write(header)
    qc_values = [str(qc_dict[metric][metric]/float(qc_dict[metric]['counted'])) for metric in qc_keys]
    fout.write('%s\t%s\n' % ('\t'.join(qc_values),opts.sampleid))
    fout.close()
