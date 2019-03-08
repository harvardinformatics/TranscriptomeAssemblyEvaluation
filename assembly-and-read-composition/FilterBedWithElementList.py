"""
creator: adam h freedman, afreedman405 at gmail dot com
date: 2016.09.26

This script filters or extracts rows from a bed file based upon strings in particular columns.
E.g. to extract only gene level entries, for the purpose of getting gene boundaries
"""

import argparse

def search_bed_line(linein,column,elements):
    bedlist=linein.strip().split()
    if bedlist[column-1] in elements:
        return True
    else:
        return False

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description="options for filtering bed file")
    parser.add_argument('-i','--bedin',dest='bedfile',type=str,help='input bed file')
    parser.add_argument('-c','--column',dest='filt_column',type=int,help='column, will be converted to python zero-based')
    parser.add_argument('-e','--elements',action='append',help='element used to filter, can be invoked multiple times to create list')
    parser.add_argument('-x','--exclude',action='store_true')
    opts = parser.parse_args()

    print opts
    counter=0
    bedin=open(opts.bedfile,'r')

    if opts.exclude==True:
        bedout=open('excluded_%s_%s' % ('.'.join(opts.elements),opts.bedfile),'w')
        for line in bedin:
            counter+=1
            if counter%10000==0:
                print 'processing bed line...',counter
            eval_bedline=search_bed_line(line,opts.filt_column,opts.elements)    
            if eval_bedline!=True:
                bedout.write(line)
    else:
        bedout=open('included_%s_%s' % ('.'.join(opts.elements),opts.bedfile),'w')
        for line in bedin:
            counter+=1
            if counter%10000==0:
                print 'processing bed line...',counter
            eval_bedline=search_bed_line(line,opts.filt_column,opts.elements)
            if eval_bedline==True:
                bedout.write(line)

    bedout.close()
