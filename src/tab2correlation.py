#!/usr/bin/env python
desc="""Calculate PARS corelation over number of .tab files
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 20/02/2012
"""

import argparse, gzip, os, sys
from datetime import datetime

import numpy             as np
from scipy.stats import stats

def load_tab(handle, load, minCount, controls):
    """Return tab-separated file as dict"""
    #deal with fname instead of fileobj
    if type(handle) is str:
        handle = open(handle)
    #open as gzip if needed
    if handle.name.endswith('.gz'):
        handle = gzip.open(handle.name)
    #load
    data = {}
    for l in handle:
        ldata = l.split('\t')
        if len(ldata) != 2:
            info = "Warning: 2 fields expected, but %s found! %s\n"
            sys.stderr.write(info % (len(ldata), str(ldata)))
            continue
        #filter controls
        geneid = ldata[0]
        if controls and geneid not in controls:
            continue
        #unload tab file
        cdata = [int(x) for x in ldata[1].strip().strip(';').split(';')]
        #skip if too small load
        if load and 1.0 * sum(cdata) / len(cdata) < load:
            continue
        #skip if too few cut points
        if minCount and len(filter(lambda x: x!=0, cdata)) < minCount:
            continue
        #store
        data[geneid] = cdata
    return data

def concatenate_scores(cgenes, gene2scores1, gene2scores2, ignore_zeros):
    """Concatenate scores from common genes"""
    scores1, scores2 = [], []
    for gene in cgenes:
        for c1, c2 in zip(gene2scores1[gene], gene2scores2[gene]):
            if ignore_zeros and not c1 and not c2:
                continue
            scores1.append(c1)
            scores2.append(c2)
    return scores1, scores2
    
def tabs2correlation(files, out, load, minCount, controlsFn, ignore_zeros, verbose):
    """Load tab files and report correlation between them"""
    #load controls
    controls = []
    if controlsFn:
       controls = set(l.strip() for l in open(controlsFn))
       if verbose:
           sys.stderr.write(" %s control names loaded.\n" % len(controls))
           
    #load tab files
    if verbose:
        sys.stderr.write("Loading %s tab files...\n" % len(files))
    tabdata = { f.name: load_tab(f, load, minCount, controls) for f in files }

    #define array to store scores
    scores  = np.zeros((len(files), len(files)))
    pvalues = np.zeros((len(files), len(files)))
    #calcute correlation
    if verbose:
        sys.stderr.write("Calculating Spearman correlation...\n")        
    for i in range(len(files)):
        #set one for self comparison for scores, pval left 0
        scores[i][i] = 1
        for j in range(i+1, len(files)):
            #get fnames
            fn1, fn2 = files[i].name, files[j].name
            #get common genes
            cgenes = set(tabdata[fn1].keys()).intersection(tabdata[fn2])
            if cgenes:
                #get scores
                scores1, scores2 = concatenate_scores(cgenes, tabdata[fn1], tabdata[fn2], ignore_zeros)
                #get correlation
                rho, pval = stats.spearmanr(scores1, scores2) #stats.pearsonr
            else:
                scores1 = []
                rho, pval = 0, 1
            if verbose:
                info = " %s - %s: %s common positions in %s genes with %.3f %s\n"
                sys.stderr.write(info%(fn1, fn2, len(scores1), len(cgenes), rho, pval))
            #store rho and number of common genes
            scores[i][j] = rho
            scores[j][i] = len(cgenes)
            #store p-value and number of positions
            pvalues[i][j] = pval
            pvalues[j][i] = len(scores1)
    #save header
    header = "\t".join(f.name for f in files)
    #save array for correlation coefficient
    title = "Spearman correlation coefficient vs number of genes"
    np.savetxt(out, scores,  delimiter="\t", header="%s\n%s"%(title, header), fmt='%.3f')
    # and for P-values
    title = "\nSpearman correlation P-value vs number of compared positions"
    np.savetxt(out, pvalues, delimiter="\t", header="%s\n%s"%(title, header))#, fmt='%.1e')

def main():

    usage  = "%(prog)s [options] -i file1.tab file2.tab"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", "--verbose",      default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", "--input",        nargs="+", type=file,
                        help="input files")
    parser.add_argument("-o", "--out",          default=sys.stdout, 
                        help="output stream      [stdout]" )
    parser.add_argument("-c", "--controls",     default="", 
                        help="compare only controls loaded from file")
    parser.add_argument("-s", "--ignore_zeros", default=False, action="store_true",
                        help="ignore zero values [%(default)s]")
    parser.add_argument("-l", "--load",         default=1.0, type=float,
                        help="min number of cuts per bp of transcript [%(default)s]")
    parser.add_argument("-m", "--minCount",     default=1, type=int,
                        help="min positions to compare [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #process two tabs files
    tabs2correlation(o.input, o.out, o.load, o.minCount, o.controls, o.ignore_zeros, o.verbose)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
