#!/usr/bin/env python
desc="""Report pars scores for each transcript position.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 21/11/2012
"""

import argparse, gzip, os, sys
from math     import log
import numpy  as np
from datetime import datetime
from counts2pars import load_pars, _check_load, normalise_counts

def load_pearson( handle ):
    """Return dictionary with pearson coefficients
    and stdev for each gene.
    """
    gene2pearson = {}
    for l in handle:
        gene,coeff,stdev = l.split()
        coeff,stdev = float(coeff),float(stdev)
        gene2pearson[gene] = ( coeff,stdev )
    return gene2pearson

def process_gene( output,gene2pearson,s1name2counts,v1name2counts,s1Count,v1Count,loadTh,verbose ):
    """
    """
    i=k=0
    coeffs = [] #,stdevs = [],[]
    for gene in s1name2counts:
        s1 = s1name2counts[gene]
        v1 = v1name2counts[gene]
        i += 1
        #check if both samples passes load criteria
        if not _check_load( s1,v1,loadTh ):
            continue
        k+=1
        if gene not in gene2pearson:
            if verbose:
                sys.stderr.write( " Warning: %s not found in gene2pearson!\n" % gene )
            continue
        #add coeff and stdev
        coeff,stdev = gene2pearson[gene]
        coeffs.append( coeff )
        #stdevs.append( stdev )

    #write output
    line = "%s\t%s\t%s\t%s\n" % ( loadTh,k,np.mean(coeffs),np.std(coeffs) )
    output.write( line )
    
def counts2pearson( s1fn,v1fn,pearson,ctrl,loadThs,output,verbose ):
    """Calculate log2(v1/s1) for every base from s1/v1.
    Skip transcripts having load below loadTh in any of the samples.
    """
    #first load counts
    if verbose:
        sys.stderr.write("Loading pearson coefficients...\n")
    gene2pearson = load_pearson( pearson ); print len(gene2pearson),gene2pearson.keys()[:10]
    #return

    #then load counts
    if verbose:
        sys.stderr.write("Loading counts...\n")
    s1name2counts,s1Count = load_pars( s1fn,1 )
    v1name2counts,v1Count = load_pars( v1fn,1 )
    
    #normalise if control provided
    if ctrl:
        if verbose:
            sys.stderr.write("Normalising by control alignments...\n")
        c0name2counts,c0Count = load_pars( ctrl.name,1 )
        s1name2counts,v1name2counts = normalise_counts( s1name2counts,v1name2counts,c0name2counts,s1Count,v1Count,c0Count,verbose )

    #then process transcripts
    if verbose:
        sys.stderr.write("Selecting passing gene... \n" )
    output.write( "#laod\tgenes\tmean\tstdev\n" )
    for loadTh in loadThs:
        process_gene( output,gene2pearson,s1name2counts,v1name2counts,s1Count,v1Count,loadTh,verbose )
        
def main():

    usage  = "%(prog)s [options] -i S1.counts V1.counts -j pearson.txt"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="files",   nargs=2, type=file,
                        help="input files        [%(default)s]")
    parser.add_argument("-j", dest="pearson", type=file,
                        help="correlation file   [%(default)s]")
    parser.add_argument("-o", dest="output",  default=sys.stdout, type=argparse.FileType('w'),
                        help="output             [stdout]")
    parser.add_argument("-l", dest="loads",    nargs="+", default=[1.0], type=float,
                        help="load thresholds    [%(default)s]")
    parser.add_argument("-c", dest="ctrl",    type=file,
                        help="normalize with control counts [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #process two tabs files
    fn1,fn2 = [ f.name for f in o.files ]
    counts2pearson( fn1,fn2,o.pearson,o.ctrl,o.loads,o.output,o.verbose )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
