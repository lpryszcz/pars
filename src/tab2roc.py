#!/usr/bin/env python
desc="""Calculate ROC curve for two tab files with increasing load. 
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

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt


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

def load_filter(data, load=0):
    """Return tabdata filtered by load cut-off"""
    newdata = {}
    for gene, scores in data.iteritems():
        if load and 1.0 * sum(scores) / len(scores) < load:
            continue
        newdata[gene] = scores
    return newdata
    
def tabs2correlation(files, out, minload, maxload, stepload, controlsFn, \
                     ignore_zeros, verbose):
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
    tabdata = { f.name: load_tab(f, 0, 0, controls) for f in files }

    #start figure
    #http://stackoverflow.com/questions/7733693/matplotlib-overlay-plots-with-different-scales
    n, r = len(files), 2
    ##n! / ((n-r)!*r!)
    ncombinations = math.factorial(n) / (math.factorial(n-r)*math.factorial(r)) 
    #fig, ax = plt.subplots(1, ncombinations, sharey='all')
    fig, ax = plt.subplots(len(files)-1, len(files)-1, sharey='row', sharex='col')

    #define loads
    loads = np.arange(minload, maxload, stepload)
    
    #calcute correlation
    if verbose:
        sys.stderr.write("Calculating Spearman correlation...\n")
    axi = 0
    for i in range(len(files)):
        #set one for self comparison for scores, pval left 0
        for j in range(i+1, len(files)):
            #get fnames
            fn1, fn2 = files[i].name, files[j].name
            if verbose:
                sys.stderr.write(" %s - %s\n" % (fn1, fn2))
            scores, genes = [], []
            for load in loads:
                #filter tabdata
                tabs = [load_filter(tabdata[fn1], load), load_filter(tabdata[fn2], load)]
                #get common genes
                cgenes = set(tabs[0].keys()).intersection(tabs[1])
                if cgenes:
                    #get scores
                    scores1, scores2 = concatenate_scores(cgenes, tabs[0], \
                                                          tabs[1], ignore_zeros)
                    #get correlation
                    rho, pval = stats.spearmanr(scores1, scores2) 
                else:
                    scores1 = scores2 = []
                    rho, pval = 0, 1
                if verbose:
                    info = " LOAD>=%s: %s common positions in %s genes with %.3f %s\n"
                    sys.stderr.write(info%(load, len(scores1), len(cgenes), rho, pval))
                #store rho and number of common genes
                scores.append(rho)
                genes.append(len(cgenes))

            #add subplot with two X-axes
            #http://matplotlib.org/examples/axes_grid/demo_parasite_axes2.html
            host = ax[i][j-1] #ax[axi] #
            par1 = host.twinx()

            #share secondary y axis
            if axi:
                host.get_shared_y_axes().join(par1, ppar1)
            ppar1 = par1

            #set log scale
            par1.set_yscale('log')
            par1.set_ylim(10, 10000)
            host.set_ylim(0.5, 1.0)
            
            #plot rho and number of genes
            p1, = host.plot(loads, scores, "b", label="Rho")
            p2, = par1.plot(loads, genes,  "r", label="Genes")
            
            #set axes labes
            ##x for bottow row
            if i+2 == len(files):
                host.set_xlabel("LOAD")
            ##y-left  for left most column
            #if axi == 0:
            if j==1:
                host.set_ylabel("Spearman correlation coefficient [rho]")
                host.yaxis.label.set_color(p1.get_color())
            ##y-right for right most column
            #if axi+1 == ncombinations:
            if j+1 == len(files):
                par1.set_ylabel("No. of genes")
                par1.yaxis.label.set_color(p2.get_color())

            #host.legend()

            #set title
            s1 = os.path.basename(fn1).split('.')[0]
            s2 = os.path.basename(fn2).split('.')[0]
            host.set_title("%s vs %s" % (s1, s2))

            #next plot axi
            axi += 1
    #plot
    #plt.draw()
    if out == sys.stdout:
        plt.show()
    else:
        plt.savefig(out, dpi=200, orientation="landscape")

def main():

    usage  = "%(prog)s [options] -i file1.tab file2.tab"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
    
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
    parser.add_argument("--minLoad",      default= 0.0, type=float,
                        help="min load [%(default)s]")
    parser.add_argument("--maxLoad",      default=20.1, type=float,
                        help="max load [%(default)s]")
    parser.add_argument("--stepLoad",     default= 0.5, type=float,
                        help="load step [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #process two tabs files
    tabs2correlation(o.input, o.out, o.minLoad, o.maxLoad, o.stepLoad, o.controls, \
                     o.ignore_zeros, o.verbose)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
