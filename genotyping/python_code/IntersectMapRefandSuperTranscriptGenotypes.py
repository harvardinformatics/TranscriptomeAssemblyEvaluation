from Bio import SeqIO
import argparse
from collections import defaultdict
from sets import Set

def CollapseGenotype(ref_allele,alt_allele_string,genotype_dictionary):
        IntToNucleotideMap = {'0':ref_allele}
        alt_allele_list=alt_allele_string.split(',')
        for i in range(len(alt_allele_list)):
            IntToNucleotideMap[str(i+1)] = alt_allele_list[i]
        genotype_integers = Set(genotype_dictionary['GT'].split('/'))
        alleles  = []
        for integer in genotype_integers:
            alleles.append(IntToNucleotideMap[integer])
        return IntToNucleotideMap['0'],alleles
        
def GenotypeLineParse(line,fields):
        line_dict = dict(zip(fields,line.strip().split()))
        gtype_dict = dict(zip(line_dict['gtformats'].split(':'),line_dict['gtdata'].split(':')))
        refallele,alleles =CollapseGenotype(line_dict['ref'],line_dict['alt'],gtype_dict)
        return line_dict,gtype_dict,alleles,refallele

def BuildCovDict(maprefcovbed):
    mapref_covdict = {}
    fopen=open(maprefcovbed,'r')
    for line in fopen:
        linelist=line.strip().split()
        mapref_covdict['%s:%s' % (linelist[0],int(linelist[1])+int(linelist[3]))]=linelist[4]
    return mapref_covdict

def MergeClusters(cluster_dict):
    popped_ids = []
    ids = cluster_dict.keys()
    ids.sort()
    for i in range(len(ids)):
        for j in range(len(ids)):
            if i != j and ids[i] not in popped_ids and ids[j] not in popped_ids:
                found = 0
                for position in cluster_dict[ids[j]]:
                    if position in cluster_dict[ids[i]]:
                        found +=1
                if found > 0:
                    cluster_dict[ids[i]] = cluster_dict[ids[i]].union(cluster_dict[ids[j]])
                    cluster_dict.pop(ids[j])
                    popped_ids.append(ids[j])

    return cluster_dict
        

def CrossMapContigsToGenomes(maprefbed,supertsbed,mreffields,superfields,maprefcovdict):
    mrin = open(maprefbed,'r') 
    superin =  open(supertsbed,'r')
   
    mref_gtype_dict = {}
    false_negative_dict = {}
    sequence_clusters = {}
    sequence_clusters = {}
        
    ref_to_superts = defaultdict(list)
    
    cluster_counter = 0
    superts_alleles_depth_dict = defaultdict(list) # need to do this bec ind genome positions may have more than 1 gtype
        
    for line in superin:
        linedict,gtypedict,alleles,ref_allele = GenotypeLineParse(line,superfields)
        superts_alleles_depth_dict['%s:%s' % (linedict['contigid'],linedict['cpos'])].append({'refallele': ref_allele,'alleles' : alleles ,'depth' : linedict['depth']})
        ref_to_superts['%s:%s' % (linedict['gchrom'],linedict['gpos'])].append('%s:%s' % (linedict['contigid'],linedict['cpos'])) 

        positions = ['%s:%s' % (linedict['gchrom'],linedict['gpos']),'%s:%s' % (linedict['contigid'],linedict['cpos'])]
        found = 0
        for clusterkey in sequence_clusters:
            for position in positions:
                if position in sequence_clusters[clusterkey]:
                    sequence_clusters[clusterkey] = sequence_clusters[clusterkey].union(positions) 
                    found += 1

        if found == 0:
            cluster_counter +=1
            sequence_clusters['snp%s' % cluster_counter]=Set(positions)
    false_negative_dict = {}
    for line in mrin:
        linedict,gtypedict,alleles,ref_allele = GenotypeLineParse(line,mreffields)
        chrompos = '%s:%s' % (linedict['gchrom'],linedict['gpos'])
        mref_gtype_dict[chrompos] = {'refallele':ref_allele,'alleles' : alleles, 'depth':maprefcovdict[chrompos]}
        if chrompos not in ref_to_superts:
            false_negative_dict[chrompos] = {'refallele':ref_allele,'alleles' : alleles, 'depth':maprefcovdict[chrompos]}

    sequence_clusters_merged=MergeClusters(sequence_clusters)

    return sequence_clusters_merged,superts_alleles_depth_dict,mref_gtype_dict,false_negative_dict

def BuildSnpTableFromSnpClusters(snp_id_dict,refalleles,supertsalleles,genome_fasta_dict,mref_covdict,false_negative_dict,outfile='test.tsv'):
    """
    refalleles == dict w/ refsnp pos as key, dict as val, with alleles and depth key-val pairs
    supertsalleles == same as for refalleles but each genomic position has as value a list of dictionaries denoting alleles and depths
    """
    table_out = open(outfile,'w')
    table_out.write('snpid\tgenomicpositions\tsupertspositions\trefnucleotides\tmaprefalleles\tsupertsalleles\tmaprefcov\n')
    for snpid in snp_id_dict:
        mapref_alleles_set = Set()
        ref_alleles_set=Set()
        mapref_coverage = []
        superts_alleles_set = Set()
        genomic_positions = Set()
        superts_positions = Set()
        for genotype_position in snp_id_dict[snpid]:
            if 'TRINITY' in genotype_position:
                superts_positions.add(genotype_position)
                for allele_dict in supertsalleles[genotype_position]:
                    for allele in allele_dict['alleles']:
                        superts_alleles_set.add(allele)
            else:
                genomic_positions.add(genotype_position) # this deals with multi-mapping of supertranscripts to genome allowing for > 1 genomic position
                if genotype_position in refalleles:
                    for allele in refalleles[genotype_position]['alleles']:
                        mapref_alleles_set.add(allele)
                        ref_alleles_set.add(refalleles[genotype_position]['refallele'])
                    mapref_coverage.append(mref_covdict[genotype_position])
                else:
                    mapref_alleles_set.add('NA')
                    chrom = genotype_position.split(':')[0]
                    pos = int(genotype_position.split(':')[1])
                    ref_alleles_set.add(str(genome_fasta_dict[chrom][pos-1]))
                    mapref_coverage.append(mref_covdict[genotype_position])
        table_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (snpid,';'.join(genomic_positions),';'.join(superts_positions),';'.join(ref_alleles_set),';'.join(mapref_alleles_set),';'.join(superts_alleles_set),';'.join(mapref_coverage)))
    for pos in false_negative_dict:
        table_out.write('genomic%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (pos,pos,'NA',false_negative_dict[pos]['refallele'],';'.join(false_negative_dict[pos]['alleles']),'NA',false_negative_dict[pos]['depth']))
    table_out.close()
            
        


        
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Summarize map-to-ref and SuperTranscript genotype intersection')
    parser.add_argument('-m','--map-ref-gtypes-bed',dest='mapref',type=str,help='bed of map-to-ref genotypes')
    parser.add_argument('-s','--supertranscript-gtypes-bed',dest='superts',type=str,help='bed of supertranscript genotypes')
    parser.add_argument('-mrc','--mapref-coverage-bed',dest='maprefcov',type=str,help='mapref raw coverageBed -d bed')
    parser.add_argument('-o','--tsv_outfile',dest='tableout',type=str,help='output of alleles and positions by superts and mapref')
    parser.add_argument('-gf','--genome-fasta',dest='genome',type=str,help='reference genome sequence')
    opts = parser.parse_args()
    
    superts_fields = ['gchrom','gposzero','gpos','gstrand','depth','contigid','cposzero','cpos','id','ref','alt','qual','filter','info','gtformats','gtdata']
    mapref_fields = ['gchrom','gposzero','gpos','ref','alt','qual','filter','info','gtformats','gtdata']

    genome_dict=SeqIO.to_dict(SeqIO.parse(opts.genome, "fasta"))
    mapref_cov_dict=BuildCovDict(opts.maprefcov)
    snp_clusters,superts_alleles_depth_dict,mref_gtype_dict,false_negative_dict = CrossMapContigsToGenomes(opts.mapref,opts.superts,mapref_fields,superts_fields,mapref_cov_dict)
   
    BuildSnpTableFromSnpClusters(snp_clusters,mref_gtype_dict,superts_alleles_depth_dict,genome_dict,mapref_cov_dict,false_negative_dict,outfile=opts.tableout)
