#!/usr/bin/env python
"""Transpose tab and semi-colon separated data"""

import sys

outcols=[]
maxlen=0
for l in sys.stdin:
  outcols.append([])
  #split by tab
  for row in l.split('\t'):
    #split by semicolon
    for r in row.split(';'):
      outcols[-1].append(r)
      #print r

  if len(outcols[-1])>maxlen:
    maxlen=len(outcols[-1])
    
lines=[ [] for x in range(maxlen) ]    
for i in range(maxlen):
  for j in range(len(outcols)):
    cell=""
    if i < len(outcols[j]):
      cell=outcols[j][i]
    lines[i].append( cell ) 
    
for l in lines:
  print "\t".join(l)
