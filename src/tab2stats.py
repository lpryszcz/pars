#!/usr/bin/env python
desc="""Calculate stats from tab file
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 26/02/2012
"""

import argparse, gzip, os, sys
from datetime import datetime

import numpy             as np
from scipy.stats import stats
from tab2correlation import load_tab

def get_zero_tail_stats(gene2scores, load=1.0):
    """Return mean, stdev, min and max 0-tail for given load. """
    zeros = []
    for gene, scores in gene2scores.iteritems():
        #load filter
        if load and  1.0 * sum(scores) / len(scores) < load:
            continue
        #count tailing zeros
        i, k = -1, 0
        while not scores[i]:
            k += 1
            i -= 1
        #store count
        zeros.append(k)
    return np.mean(zeros), np.std(zeros), min(zeros), max(zeros)

def tab2stats(tab, out, controls, load, verbose):
    """Load tab files and report correlation between them"""
    #load = minCount = 0
    tabdata = load_tab(tab, 0, 0, controls) 

    #get length of 0-tail for load=1
    mean, stdev, minc, maxc = get_zero_tail_stats(tabdata, load)
    out.write("%s\t%s\t%s\t%s\t%s\n" % (tab.name, mean, stdev, minc, maxc))
    

def main():

    usage  = "%(prog)s [options] -i file1.tab file2.tab"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", "--verbose",      default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", "--input",        nargs="+", type=file,
                        help="input files")
    parser.add_argument("-o", "--out",          default=sys.stdout, 
                        help="output stream      [stdout]" )
    parser.add_argument("-l", "--load",         default=1.0, type=float,
                        help="min number of cuts per bp of transcript [%(default)s]")
    parser.add_argument("-c", "--controls",     default="", 
                        help="compare only controls loaded from file")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #load controls
    controls = []
    if o.controls:
       controls = set(l.strip() for l in open(controlsFn))
       if o.verbose:
           sys.stderr.write(" %s control names loaded.\n" % len(controls))

    #process tab files
    if o.verbose:
        sys.stderr.write("Processing %s tab files...\n" % len(o.input))
        
    header = "fname\tmean\tstdev\tmin\tmax"
    o.out.write("# 0-tail for LOAD>=%s\n%s\n"%(o.load, header))        
    for tab in o.input:
        tab2stats(tab, o.out, controls, o.load, o.verbose)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
