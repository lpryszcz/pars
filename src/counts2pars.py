#!/usr/bin/env python
desc="""Report pars scores for each transcript position.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Dublin, 2/07/2012
"""

import argparse, gzip, os, sys
from datetime import datetime
from math     import log

def load_pars( fn,getTotCount=0 ):
    """Return dictionary of counts for every transcript.
    In addition, return totCount from given experiment in millions.
    """
    name2counts = {}
    #open file
    if fn.endswith(".gz"):
        handle = gzip.open(fn)
    else:
        handle =     open(fn)
    #load scores/counts
    totCount = 0
    for l in handle:
        l = l.strip()
        l = l.strip(";")
        name,length,counts = l.split("\t")[:3]
        c = [ float(x) for x in counts.split(";") ]
        name2counts[name] = c
        totCount += sum(c)
    if getTotCount:
        return ( name2counts,totCount/10.0e6 )
    return name2counts

def _get_load( s1,v1 ):
    """Return True if s1 and v1 pass load threshold.
    """
    #if not sum(s1) or not sum(v1):
    #    return False
    eff_length = len(s1)
    s1load = sum(s1)/eff_length
    v1load = sum(v1)/eff_length    
    #if s1load>=loadTh and v1load>=loadTh:
    #    return True
    return s1load,v1load

def normalise_counts( s1name2counts,v1name2counts,c0name2counts,s1Count,v1Count,c0Count,verbose ):
    """Normalise by reads from control library.
    """
    s1n2c = {}
    v1n2c = {}
    for gene in s1name2counts:
        #get counts for gene
        s1 = s1name2counts[gene]
        v1 = v1name2counts[gene]
        c0 = c0name2counts[gene]
        
        #normalise each position count by control
        s1n = []
        v1n = []
        for si,vi,ci in zip( s1,v1,c0 ):
            #normalize by total number of reads and by control
            cn = ci/c0Count
            #sn = ( si/s1Count - cn ) * s1Count
            sn = si - ci*s1Count/c0Count
            vn = vi - ci*v1Count/c0Count
            #zero counts < 0 and append to lists
            if sn < 0:
                sn = 0
            s1n.append( sn )    
            if vn < 0:
                vn = 0
            v1n.append( vn )

        #store
        s1n2c[gene] = s1n
        v1n2c[gene] = v1n

    return ( s1n2c,v1n2c )
            
def process_gene( output,s1name2counts,v1name2counts,s1Count,v1Count,loadTh,minReads,readsPerSample,verbose ):
    """
    """
    i=k=0
    for gene in s1name2counts:
        s1 = s1name2counts[gene]
        v1 = v1name2counts[gene]
        i += 1
        if verbose:
            sys.stderr.write(" %s %s %s    \r" % (i,k,gene) )        
        #check if both samples passes load criteria
        s1load,v1load = _get_load( s1,v1 )
        if s1load < loadTh or v1load < loadTh:
            continue
        k+=1
        #calculate pars score for each base
        pars = []
        for si,vi in zip( s1,v1 ):
            if si and vi and si>=readsPerSample and vi>=readsPerSample and si+vi>=minReads:
                #p = "%.2f" % round( log( vi*s1Count / si*v1Count,2 ),2 )
                p = "%.2f" % round( log( vi / si,2 ),2 )
            else:
                p = "0"
            pars.append( p )
        #write output
        line = "%s\t%s\t%s\t%s\t%s\n" % ( gene,len(pars),";".join( pars ),s1load,v1load )
        output.write( line )
    return ( i,k )
    
def counts2pars( s1fn,v1fn,ctrl,loadTh,output,minReads,readsPerSample,verbose ):
    """Calculate log2(v1/s1) for every base from s1/v1.
    Skip transcripts having load below loadTh in any of the samples.
    """
    #first load counts
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
        sys.stderr.write("Calculating PARS scores...\n")
    i,k = process_gene( output,s1name2counts,v1name2counts,s1Count,v1Count,loadTh,minReads,readsPerSample,verbose )
    
    if verbose:
        sys.stderr.write("Processed %s transcripts. %s passed load filter.\n" % (i,k) )        
        
def main():

    usage  = "%(prog)s [options] [options] -i S1.counts V1.counts"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="files",   nargs=2, type=file,
                        help="input files        [%(default)s]")
    parser.add_argument("-o", dest="output",  default=sys.stdout, type=argparse.FileType('w'),
                        help="input              [stdout]")
    parser.add_argument("-l", dest="load",    default=1.0, type=float,
                        help="minimum load       [%(default)s]")
    parser.add_argument("-c", dest="ctrl",    type=file,
                        help="normalize with control counts [%(default)s]")
    parser.add_argument("-t", dest="minReads", default=0, type=int,
                        help="min s reads at position in all samples [%(default)s]")
    parser.add_argument("-s", dest="readsPerSample", default=0, type=int,
                        help="min s reads at position in each sample [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #process two tabs files
    fn1,fn2 = [ f.name for f in o.files ]
    counts2pars( fn1,fn2,o.ctrl,o.load,o.output,o.minReads,o.readsPerSample,o.verbose )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
