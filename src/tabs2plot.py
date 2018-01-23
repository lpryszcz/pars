#!/usr/bin/env python
desc="""Plot Pearson corelation between PARS scores from 2 sources.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Dublin, 4/07/2012
"""

import argparse, gzip, os, sys
from datetime import datetime
import matplotlib.pyplot as plt
import numpy             as np
from scipy.stats.stats import pearsonr
from counts2pars import load_pars

def plotData(x,y,color="red",marker="o",name="Data"):
    "From: http://stackoverflow.com/questions/8154511/drawing-a-correlation-graph-in-matplotlib"
    fit    = np.polyfit(x,y,1)#; print fit
    fit_fn = np.poly1d(fit)#; print fit_fn # fit_fn is now a function which takes in x and returns an estimate for y
    #pearR  = np.corrcoef(x,y)[1,0]
    pearR,pval = pearsonr(x,y)
    plt.plot( x,y,"%s%s" % (color[0],marker) )
    #"--k",color=color,
    plt.plot( x,fit_fn(x),label="y=%s; r=%.3f;\nP=%s" % (str(fit_fn).strip(),pearR,pval) )
    return pearR,pval

def skip_zeros( c1,c2 ):#*argv ):
    """Return tuple of *argv lists with only ith elements
    that are non-zero in all list from *argv."""
    c11,c22 = [],[]
    for a,b in zip(c1,c2):
        if a and b:
            c11.append(a)
            c22.append(b)
    return c11,c22

def scatter_plot( outdir,gene,c1,c2,fn1,fn2,color="blue",marker="." ):
    """Save scatter-plot of c1,c2 values"""
    #plotting stuff
    plt.clf()
    #add data
    pearR,pval = plotData(c1,c2,color,marker,gene)
    #add title and so on
    plt.title( "Correlation of PARS scores for %s" % gene )
    plt.xlabel( 'PARS from %s' % fn1 )
    plt.ylabel( 'PARS from %s' % fn2 )
    plt.legend(loc=2)
    #plt.show()
    #save fig
    outfn = os.path.join( outdir,gene+".png" )
    plt.savefig( outfn )
    return pearR,pval
    
def plot_pars( fn1,fn2,outdir,ignore_zeros,mincount,verbose ):
    """Calculate log2(v1/s1) for every base from s1/v1.
    Skip transcripts having load below loadTh in any of the samples.
    """
    #first load counts
    if verbose:
        sys.stderr.write("Loading counts...\n")
    counts1 = load_pars( fn1 )
    counts2 = load_pars( fn2 )

    #then process transcripts
    if verbose:
        sys.stderr.write("Plotting PARS scores for %s entries...\n" % len(counts2) )
    i=j=pos=0
    coeffs = []
    for gene in sorted(counts2.keys()):
        i+=1
        if gene not in counts1:
            continue
        if verbose:
            sys.stderr.write(" %s %s %s    \r" % (i,j,gene) )        
        c1 = counts1[gene]
        c2 = counts2[gene]
        if len(c1)!=len(c2):
            sys.stderr.write("Error: Different transcript length for %s (%s,%s)\n" % (gene,len(c1),len(c2)) )
            continue
            
        #skip missing data 
        if ignore_zeros:
            c1,c2 = skip_zeros( c1,c2 )
            if len(c1)<mincount:
                continue
        j += 1
        #sys.stderr.write("%s\n%.2f\n%.2f\n" % ( gene,"\t".join(str(x) for x in c1),"\t".join(str(x) for x in c2) ) )
        #save scatterplot and return pearson coefficient and P
        coeff,p = scatter_plot( outdir,gene,c1,c2,fn1,fn2,color="blue",marker="." )
        print "%s\t%s\t%s\t%s" % ( gene,coeff,p,len(c1) )
        coeffs.append( coeff )
        pos += len(c1)
        
    sys.stderr.write("Correlation mean: %.2f [+- %.2f] for %s positions\n\n" % ( np.mean(coeffs),np.std(coeffs),pos ) )
       
def main():

    usage  = "%(prog)s [options] tab1 tab2"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="files",   nargs=2, type=file,
                        help="input files        [%(default)s]")
    parser.add_argument("-o", dest="outdir",  default="pars_pics", 
                        help="output directory   [%(default)s]" )
    parser.add_argument("-s", dest="ignore_zeros", default=False, action="store_true",
                        help="ignore zero values [%(default)s]")
    parser.add_argument("-l", dest="mincount",     default=5, type=int,
                        help="min positions to compare [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    if not os.path.isdir(o.outdir):
        os.makedirs(o.outdir)

    #process two tabs files
    fn1,fn2 = [ f.name for f in o.files ]
    plot_pars( fn1,fn2,o.outdir,o.ignore_zeros,o.mincount,o.verbose )
    
if __name__=='__main__': 
    t0 = datetime.now()
    main()
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )

  
