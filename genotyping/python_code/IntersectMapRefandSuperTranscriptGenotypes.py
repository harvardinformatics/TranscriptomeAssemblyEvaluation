import argparse
from collections import defaultdict
from sets import Set

def CollapseGenotype(ref_allele,alt_allele_string,genotype_dictionary):
        IntToNucleotideMap = {'0':ref_allele}
        alt_allele_list=alt_allele_string.split(',')
        for i in range(len(alt_allele_list)):
            IntToNucleotideMap[str(i+1)] = alt_allele_list[i]
        genotype_integers = Set(genotype_dictionary['GT'].split('/'))
        gtypes  = []
        for integer in genotype_integers:
            gtypes.append(IntToNucleotideMap[integer])
        return gtypes
        
def GenotypeLineParse(line,fields):
        line_dict = dict(zip(fields,line.strip().split()))
        gtype_dict = dict(zip(line_dict['gtformats'].split(':'),line_dict['gtdata'].split(':')))
        alleles = CollapseGenotype(line_dict['ref'],line_dict['alt'],gtype_dict)
        return line_dict,gtype_dict,alleles

def CrossMapContigsToGenomes(maprefbed,supertsbed,mreffields,superfields):
        ref_to_superts = defaultdict(list)
        superts_to_ref = defaultdict(list)

        mrin = open(maprefbed,'r')    
        mref_gtype_dict = {}
        superin =  open(supertsbed,'r')

        superts_alleles_depth_dict = {}

        for line in superin:
            linedict,gtypedict,alleles = GenotypeLineParse(line,superfields)
            ref_to_superts['%s:%s' % (linedict['gchrom'],linedict['gpos'])].append('%s:%s' % (linedict['contigid'],linedict['cpos']))
            superts_to_ref['%s:%s' % (linedict['contigid'],linedict['cpos'])].append('%s:%s' % (linedict['gchrom'],linedict['gpos']))
            superts_alleles_depth_dict['%s:%s' % (linedict['gchrom'],linedict['gpos'])] = {'alleles' : alleles ,'depth' : linedict['depth']}

        for line in mrin:
            linedict,gtypedict,alleles = GenotypeLineParse(line,mreffields)
            chrompos = '%s:%s' % (linedict['gchrom'],linedict['gpos'])
            mref_gtype_dict[chrompos] = {'alleles' : alleles}
            if chrompos not in ref_to_superts:
                ref_to_superts[chrompos] = 'NA'
    
        return ref_to_superts,superts_to_ref,superts_alleles_depth_dict,mref_gtype_dict

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='Summarize map-to-ref and SuperTranscript genotype intersection')
    parser.add_argument('-m','--map-ref-gtypes-bed',dest='mapref',type=str,help='bed of map-to-ref genotypes')
    parser.add_argument('-s','--supertranscript-gtypes-bed',dest='superts',type=str,help='bed of supertranscript genotypes')
    opts = parser.parse_args()

    superts_fields = ['gchrom','gposzero','gpos','gstrand','depth','contigid','cposzero','cpos','id','ref','alt','qual','filter','info','gtformats','gtdata']
    mapref_fields = ['gchrom','gposzero','gpos','ref','alt','qual','filter','info','gtformats','gtdata']

    ref_to_superts,superts_to_ref,superts_alleles_depth_dict,mref_gtype_dict = CrossMapContigsToGenomes(opts.mapref,opts.superts,mapref_fields,superts_fields)

    print 'ref to superts',len(ref_to_superts),ref_to_superts[ref_to_superts.keys()[0]]
    print 'superts to ref',len(superts_to_ref),superts_to_ref[superts_to_ref.keys()[0]]
    print 'superts alleles depth dict',len(superts_alleles_depth_dict),superts_alleles_depth_dict[superts_alleles_depth_dict.keys()[0]]
    print 'mref gtype dict',len(mref_gtype_dict),mref_gtype_dict[mref_gtype_dict.keys()[0]]


