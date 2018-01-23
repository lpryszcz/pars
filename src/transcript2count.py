#!/usr/bin/env python
desc="""Count number of 5'-end reads for all transcripts at given load.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Dublin/Barcelona, 4/07/2012

versions:
- 0.3 - unstranded option added (to be done!)
- 0.2 - solved issue of including reverse reads for + transcripts (thanks to Jose) 
- 0.1 
"""

import argparse, os, sys
import pysam
from datetime import datetime
from genome_annotation import load_transcripts_bed
        
def get_reads_5ends( samfile,ref,start,end,mapq,reverse ):
    """Return number of reads starting at each position
    of defined region.
    """
    #define empty counts list
    counts = [ 0 for x in range(end-start) ]
    #process all reads in region
    for sam in samfile.fetch( ref,start,end+1 ):
        #skip if low quality
        if sam.mapq < mapq:
            continue
        ##skip if read comes from opposite strand
        #if the transcript is reverse and read is not reverse
        #if the transcript is forward and the read is reverse
        if sam.is_reverse and not reverse or not sam.is_reverse and reverse: 
            continue
        #get 5' position of read and 1bp upstream
        if reverse:
            pos = sam.aend - 1 
        else:
            pos = sam.pos
        #check if start in transcript range
        tpos = pos-start
        if 0 <= tpos < len(counts):
            #update counts
            counts[tpos] += 1
    return counts

def process_transcripts( bed,bam,mapq,offset,verbose ):
    """
    """
    #first load transcripts
    if verbose:
        sys.stderr.write("Loading BED file...\n")
    transcripts = load_transcripts_bed( bed )#oneoff=True ; print transcripts

    #then parse bam file for every transcript
    if verbose:
        sys.stderr.write("Parsing BAM file for %s BED entries...\n" % len(transcripts) )
    #process transcripts
    i = 0
    samfile = pysam.Samfile( bam,"rb" )
    refs    = set( samfile.references )
    for transcript in transcripts:
        #if transcript not in ("YAL001C","YAR002W","YAR002C-A","YAR003W","YAL003W"): continue
        i += 1
        if verbose:
            sys.stderr.write(" %s %s    \r" % (i,transcript) )
        #ref
        ref    = transcripts[transcript]["chromosome"]
        strand = transcripts[transcript]["strand"]
        #check if ref indeed in bam
        if ref not in refs:
            if verbose:
                sys.stderr.write( " Warning: %s not in BAM file\n" % ( ref, ) )
            continue
        #define strand
        reverse = False
        if strand == "-":
            reverse = True            
        #process all exons (and UTRs)
        counts = []
        exonCount = len(transcripts[transcript]["intervals"])
        for ii in range(exonCount):
            #get start and end
            start,end,score = transcripts[transcript]["intervals"][ii]
            #add 1bp for transcript start if reverse
            if ii   == 0 and reverse:
                start -= offset
            #add 1bp for transcript end if not reverse
            if ii+1 == exonCount and not reverse:
                end   += offset
            #get counts
            counts += get_reads_5ends( samfile,ref,start,end,mapq,reverse )
        #reverse (or not) and add tailing base count
        if reverse:
            counts.reverse()
        #move all by 1bp - as 1bp added to transcript start)
        counts = counts[offset:]
        #output line
        line   = "%s\t%s\t%s\n" % ( transcript,len(counts),";".join(str(x) for x in counts) ) 
        sys.stdout.write( line )

def main():

    usage  = "%(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='0.2')
    parser.add_argument("-a", dest="bam",  type=file,
                        help="bam file (sorted)     [%(default)s]")
    parser.add_argument("-b", dest="bed",  type=file,
                        help="bed file              [%(default)s]")
    parser.add_argument("-q", dest="mapq", default=0, type=int,
                        help="min mapping quality   [%(default)s]")
    parser.add_argument('--offset', dest='offset', default=1,
                        help="reads are counted for base upstream [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #process all transcripts
    process_transcripts( o.bed.name,o.bam.name,o.mapq,o.offset,o.verbose )
        
if __name__=='__main__': 
    t0 = datetime.now()
    main()
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
  
