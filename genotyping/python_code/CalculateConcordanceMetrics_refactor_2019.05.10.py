import argparse
from sets import Set

fields = ['snpid', 'genomicpositions', 'supertspositions', 'mapref_ref', 'maprefalleles', 'superts_ref', 'supertsalleles', 'maprefcov', 'mapref_ad', 'superts_ad','mean_supercov','median_supercov','SToverlaps','SToverlapalleles','proteinhit']

##### Updating Allele Sets #####

def ReplaceNaWithReferenceAllele(snp_dict):
    ## update MR to ref allele if no gtype call ##
    if snp_dict['maprefalleles'] == 'NA':
        snp_dict['maprefalleles'] = snp_dict['mapref_ref']
    
    ## update ST to ST ref if an ST maps to the locus ##
    if snp_dict['superts_ref'] == 'NA' and snp_dict['SToverlapalleles'] != 'NA':
        snp_dict['supertsalleles'] = snp_dict['SToverlapalleles']
   
    return snp_dict


##### Defining Metrics #####

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

def Precision(snp_dict):
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

def IndelError(snp_dict):
    """
    frequency of sites where ST gtype
    error due to ST indel polymorphism
    """  
    if snp_dict['supertsalleles'] != 'NA':
        mapref_allele_set = Set(snp_dict['maprefalleles'].split(';'))
        superts_allele_set = Set(snp_dict['supertsalleles'].split(';'))
        mapref_allele_lengths = Set([len(allele) for allele in Set(snp_dict['maprefalleles'].split(';'))])
        superts_allele_lengths = Set([len(allele) for allele in Set(snp_dict['supertsalleles'].split(';'))])
 
        if mapref_allele_set != superts_allele_set:
            if len(mapref_allele_lengths) == 1 and len(superts_allele_lengths) > 1:
                return True
            else:
                return False
    else:
        return 'NA'
        

def MrefSnvSuperTsIndel(snp_dict): # OK
    """
    frequency of sites where SuperTs gtype error
    due to a map-to-ref SNV called as an 
    indel with SuperTranscript;
    these polymorphic sites would normally
    get filtered out bec an indel
    """

    if snp_dict['supertsalleles'] != 'NA':
        mapref_allele_lengths = Set([len(allele) for allele in Set(snp_dict['maprefalleles'].split(';'))])
        superts_allele_lengths = Set([len(allele) for allele in Set(snp_dict['supertsalleles'].split(';'))])
        
        if Set(snp_dict['maprefalleles'].split(';')) != Set(snp_dict['supertsalleles'].split(';')):
            if len(mapref_allele_lengths) == 1 and len(superts_allele_lengths) > 1:
                return True
            else:
                return False
    else:
        return 'NA'


def RefAlleleIncludedInSuperTsError(snp_dict): #OK
    """
    when gtype called for both map-to-ref
    and SuperTranscripts, when SuperTs gtype error,
    whether a map-to-ref allele is included in 
    the erroneous genotype
    """
    if snp_dict['supertsalleles'] != 'NA':
        mapref_allele_set = Set(snp_dict['maprefalleles'].split(';'))
        superts_allele_set = Set(snp_dict['supertsalleles'].split(';'))
        if mapref_allele_set != superts_allele_set:
            if len(mapref_allele_set.intersection(superts_allele_set)) > 0:
                return True
            else:
                return False
    else:
        return 'NA'    


def FalsePositiveIndel(snp_dict): #OK 
    """
    frequency of false positives that
    are called as indels
    """
    if len(Set(snp_dict['maprefalleles'].split(';'))) == 1:
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


def Recall(snp_dict): 
    """
    determine if a called map-ref genotype, including
    fixed-alternative, is recovered by 
    the SuperTranscript approach
    """
    if snp_dict['maprefalleles'] != snp_dict['mapref_ref']:
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
    mapref_allele_set = Set(snp_dict['maprefalleles'].split(';'))
    superts_allele_set = Set(snp_dict['supertsalleles'].split(';'))
    if len(mapref_allele_set) == 2 and mapref_allele_set == superts_allele_set:
        return True
    else:
        return False

    else:
        return 'NA'


def FalsePositive(snp_dict):
    """
    cases where no genotype, or a fixed alternative site
    in map-to-ref is called as a bi-allelic het via
    SuperTranscript method
    """

    if snp_dict['supertsalleles'] =! 'NA':
        if len(Set(snp_dict['maprefalleles'].split(';'))) == 1 and len(Set(snp_dict['supertsalleles'].split(';'))) >= 2:
            return True
        elif Set(snp_dict['maprefalleles'].split(';'))) == Set(snp_dict['supertsalleles'].split(';'))):
            return False
        else:
            return 'NA'
    else:
        return 'NA'


def FalseNegative(snp_dict):
    """
    sites:
        -- no overlapping ST allele
        -- invariant ST and het MR
    """
    if snp_dict['maprefalleles'] != snp_dict['mapref_ref']:
        refbase = False
    else:
        refbase = True
    
    if snp_dict['supertsalleles'] == 'NA':
        FN = True
    elif len(Set(snp_dict['supertsalleles'].split(';'))) == 1 and len(Set(snp_dict['maprefalleles'].split(';'))) == 2:
        FN = True
    elif Set(snp_dict['supertsalleles'].split(';'))) == Set(snp_dict['maprefalleles'].split(';'))):
        FN = False
    else:
        FN = 'NA'
    return FN,refbase

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='Converts map-to-ref vcf file from RNA-seq data to exonic genotypes in bed format')
    parser.add_argument('-gt','--gtable',dest='genotypes',type=str,help='table of mapref and superts genotype intersections')
    parser.add_argument('-sid','--sample_id', dest='sampleid',type=str,help='sample label for table writing')
    opts = parser.parse_args()
    
    table_corrected = open('corrected_%s' % opts.genotypes,'w')
    qc_dict = {'fp_indel':{'fp_indel': 0,'counted': 0},
               'err_snv2indel' : {'err_snv2indel': 0,'counted': 0},
               'err_mrefincl' : {'err_mrefincl': 0,'counted': 0},
               'precision' : {'precision': 0,'counted': 0},
               'fn' : {'fn' : 0,'counted' : 0},
               'fn_refbase_incl': {'fn_refbase_incl' : 0 , 'counted' : 0},
               'fp' : {'fp' : 0,'counted' : 0},
               'recall' : {'recall' : 0,'counted' : 0},
               'het_recall' : {'het_recall' : 0,'counted' : 0}
               'concordance' : {'concordance' : 0,'counted' : 0},
               'fn_noalignment' : {'fn_noalignment' : 0,'counted' : 0}
               }

    fopen = open(opts.genotypes,'r')
    header_fields = fopen.readline().strip().split()
    table_corrected.write('%s\n' % '\t'.join(header_fields))

    for line in fopen:
        line_list = line.strip().split('\t')
        snp_dict =  ReplaceNaWithReferenceAllele(dict(zip(fields,line_list)))
        table_corrected.write('%s\n' % '\t'.join([snp_dict[field] for field in header_fields]))

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

        ### precision ###
        precision = Precision(snp_dict)
        if precision == True:
            qc_dict['precision']['precision'] +=1
            qc_dict['precision']['counted']+=1
        elif precision == False:
            qc_dict['precision']['counted']+=1
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

        ### false positive ##
        #
        fp = FalsePositive(snp_dict)
        if fp == True:
            qc_dict['fp']['fp']+=1
            qc_dict['fp']['counted']+=1
        elif fp = False:
            qc_dict['fp']['counted']+=1
        else:
            pass

        ### false positives as indels ###
        fpindel = FalsePositiveIndel(snp_dict)
        if fpindel == True:
            qc_dict['fp_indel']['fp_indel']+=1
            qc_dict['fp_indel']['counted']+=1
        elif fpindel == False:
            qc_dict['fp_indel']['counted']+=1
        else:
            pass


        ### false negative ###
        fn = FalseNegative(snp_dict)
        if fn == True,False: # MR != refbase w/ no call
            qc_dict['fn']['fn']+=1
            qc_dict['fn']['counted']+=1
            qc_dict['fn_refbase_incl']['fn']+=1
            qc_dict['fn_refbase_incl']['counted']+=1
        elif fn == True,True:
            qc_dict['fn_refbase_incl']['fn']+=1
            qc_dict['fn_refbase_incl']['counted']+=1 
        elif fn == False,False:
            qc_dict['fn']['counted']+=1
            qc_dict['fn_refbase_incl']['counted']+=1
        elif fn == False,True:
            qc_dict['fn_refbase_incl']['counted']+=1
              

        if snp_dict['SToverlaps'] == 'NA':
                qc_dict['fn_noalignment']['fn_noalignment']+=1
                qc_dict['fn_noalignment']['counted']+=1
            else:
                qc_dict['fn_noalignment']['counted']+=1

    fout= open('%s.concordancemetrics.tsv' % opts.genotypes[:-4],'w')
    
    qc_keys =  qc_dict.keys()
    qc_keys.sort()
    header='%s\tSample\n' % '\t'.join(qc_keys)
    fout.write(header)
    qc_values = [str(qc_dict[metric][metric]/float(qc_dict[metric]['counted'])) for metric in qc_keys]
    fout.write('%s\t%s\n' % ('\t'.join(qc_values),opts.sampleid))
    fout.close()
