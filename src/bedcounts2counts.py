#!/usr/bin/env python
desc="""Report coverage at given position of each transcript.
Strand specific.

USAGE:
for f in TOTAL_?1.bam; do
  echo `date` $f;
  bedtools coverage -d -s -abam $f -b mRNA.bed | bedcounts2counts.py > $f.wholereads.tab;
done;
date
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Dublin, 16/07/2012
"""

import argparse, os, sys
from datetime import datetime

def bedcounts2counts( handle,out,verbose ):
    """Report coverage at each position of transcript.
    """
    t2counts = {}
    for line in sys.stdin:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        #TETp9_9.1_fwd_ref	0	95	TETp9	.	+	1	1951
        chrName,s,e,transcript,score,strand,pos,count = line.split("\t")
        s,e,pos = int(s),int(e),int(pos)
        #add new transcript
        if transcript not in t2counts:
            t2counts[transcript] = { "strand": strand, "starts": { s: [] } } #[ 0 for i in range(e-s) ]
        #add new transcript fragment
        if s not in t2counts[transcript]["starts"]:
            t2counts[transcript]["starts"][s] = []
        #add base count
        t2counts[transcript]["starts"][s].append( count )
        
    #add last transcript
    for transcript in sorted( t2counts.keys() ):
        counts=[]
        reverse=False
        if t2counts[transcript]["strand"]=="-":
            reverse=True
        for s,region_counts in sorted( t2counts[transcript]["starts"].iteritems(),reverse=reverse ):
            #sys.stderr.write( " %s %s %s %s\n" % ( transcript,t2counts[transcript]["strand"],s,region_counts[:10] ) )
            counts += region_counts
        out.write( "%s\t%s\t%s\n" % ( transcript,len(counts),";".join(counts) ) )

def main():
    usage  = "bedtools coverage -d -s -abam $f -b mRNA.bed | %(prog)s [options] S1.bam V1.bam S1.+.bed V1.+.bed"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="input",   default=sys.stdin,  type=file,
                        help="input              [stdin]")
    parser.add_argument("-o", dest="output",  default=sys.stdout, type=argparse.FileType('w'),
                        help="input              [stdout]")
   
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    bedcounts2counts( o.input,o.output,o.verbose )
    
if __name__=='__main__': 
    t0 = datetime.now()
    main()
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
  
