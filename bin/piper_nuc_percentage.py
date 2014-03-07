#! /usr/bin/env python
import sys
from collections import defaultdict, Counter
ct = defaultdict ( lambda: Counter () )
ext = int (sys.argv[1]) * 2 + 1
for line in sys.stdin:
        seq = str (line.split()[1])
        for i in range (0,ext):
                ct[i][seq[i].upper()]+=1
for n in range (0, ext):
        A = ct[n]['A']
        C = ct[n]['C']
        G = ct[n]['G']
        T = ct[n]['T']
        print str (A) + '\t' + str (C) + '\t' + str (G) + '\t' + str (T)
