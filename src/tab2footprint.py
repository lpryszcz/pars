#!/usr/bin/env python
desc="""Generate footprints plots from multiple sources.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 26/02/2012
"""

import argparse, gzip, os, sys
from datetime import datetime
import math
import numpy as np
from scipy.stats import stats

import matplotlib.pyplot as plt

def load_tab(handle, gene):
    """Return tab-separated file as list"""
    #deal with fname instead of fileobj
    if type(handle) is str:
        handle = open(handle)
    #open as gzip if needed
    if handle.name.endswith('.gz'):
        handle = gzip.open(handle.name)
    #load
    data = []
    for l in handle:
        ldata = l.split('\t')
        if len(ldata) != 2:
            info = "Warning: 2 fields expected, but %s found! %s\n"
            sys.stderr.write(info % (len(ldata), str(ldata)))
            continue
        #filter controls
        geneid = ldata[0]
        if gene!=geneid:
            continue
        #unload tab file
        cdata = [float(x) for x in ldata[1].strip().strip(';').split(';')]
        #store enrichment
        data.append([x/np.mean(cdata) for x in cdata])
    return data
    
def tab2footprint(files, out, gene, verbose):
    """ """
    #get tabdata
    tabdata = { f.name: load_tab(f, gene) for f in files }

    #get correlations
    fn_scores = [] 
    
    #get plots
    fig, axlist = plt.subplots(len(files), 1, sharex='col') # sharey='row',
    
    for i in range(len(files)):
        fn = files[i].name
        #get plot
        ax = axlist[i]
        #set title
        ax.set_title("%s" % os.path.basename(fn).split('.')[0])
        #set log scale
        #ax.set_yscale('log')
        #ax.set_ylim(0.5, 1.0)
        ax.set_ylabel("Reads enrichment")
        #add scores
        for scores in tabdata[fn]:
            fn_scores.append((fn, scores))
            #plot rho and number of genes
            p1, = ax.plot(xrange(1, len(scores)+1), scores) #, "b", label="Rho")
                          
    #set xlabel for last plot
    ax.set_xlabel("Position")

    #define array to store scores
    scores  = np.zeros((len(fn_scores), len(fn_scores)))
    #get correlations
    for i in range(len(fn_scores)):
        #set one for self comparison for scores, pval left 0
        scores[i][i] = 1
        for j in range(i+1, len(fn_scores)):
            fn1, scores1 = fn_scores[i]
            fn2, scores2 = fn_scores[j]
            if scores1 and scores2:
                rho, pval = stats.spearmanr(scores1, scores2)
            else:
                rho, pval = 0, 1
            #store rho and number of common genes
            scores[i][j] = rho
            scores[j][i] = pval
                
    #save header
    header = "\t".join(fn for fn, scores in fn_scores)
    #save array for correlation coefficient
    title = "Spearman correlation coefficient vs P value"
    np.savetxt(out, scores,  delimiter="\t", header="%s\n%s"%(title, header), fmt='%.3f')
    #return
    #plot
    if out == sys.stdout:
        plt.show()
    else:
        plt.savefig(out, dpi=200)#, orientation="landscape")
    
def main():

    usage  = "%(prog)s [options] -i file1.tab file2.tab"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
    
    parser.add_argument("-v", "--verbose",      default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", "--input",        nargs="+", type=file,
                        help="input files")
    parser.add_argument("-o", "--out",          default=sys.stdout, 
                        help="output stream      [stdout]" )
    parser.add_argument("-g", "--gene",         default="TETp9p9.1", 
                        help="report info for gene")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #process two tabs files
    tab2footprint(o.input, o.out, o.gene, o.verbose)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
    
   