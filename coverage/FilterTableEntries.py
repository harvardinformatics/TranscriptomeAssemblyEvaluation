import argparse
#'/n/holylfs/LABS/informatics/adamf/DeNovoTranscriptomeEvaluation/reference_transcript_coverage/mus_proteincoding_ts.ids','r')

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='generic filtering on protein coding database hits')
    parser.add_argument('-input','--unfiltered_table',dest='tablein',type=str,help='table to filter on transcript hits')
    parser.add_argument('-sep','--field_separator',dest='sep',type=str,help='separator definiition,using codes: TAB, SPACE, COMMA')
    parser.add_argument('-k','--keep-list',dest='keepers',type=str,help='file of entries that are protein coding')
    parser.add_argument('-c','--target_column',dest='column',type=int,help = 'nonzero indexed column where target in dbase found')
    parser.add_argument('-header','--header_present',action='store_true',help='boolean flag to indicate if header')
    parser.add_argument('-pfix','--output_prefix',dest='prefix',type=str,help='prefix for outfile')
    opts = parser.parse_args()
    
    num_filtered_out = 0
    keepers_list = []
    keeper_handle = open(opts.keepers,'r')
    for line in keeper_handle:
        keeper_list.append(line.strip())

    sep_dict = {'TAB':'\t','SPACE': ' ','COMMA':','}
    sep = sep_dict[opts.sep]

    tablein = open(opts.tablein, 'r')
    tableout = open('%s_%s' % (opts.prefix,opts.tablein),'w')
    if opts.header_present:
        header = tablein.readline()
        tableout.write(header)
    for line in tablein:
        linelist = line.strip().split(sep)
        if linelist[opts.column - 1] in pcod_list:
            tableout.write(line)
        else:
            num_filtered_out+=1
    tableout.close()
    print('%s queries removed from table' % num_filtered_out)
    

