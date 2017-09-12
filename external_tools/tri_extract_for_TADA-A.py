#!/public/home/jiangyi/bin/Python3.4.0/bin/python3.4

from __future__ import print_function
import sys
import re
import subprocess
import sys, re, getopt, math

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if sys.argv[1] == '--help':
    print('tri_extract_changed_by_YL.py [fasta_file]')
    print ('fasta_file is produced from bedtools getfasta -tab, from a bed file, the first column has the format like this: chr1:10001-10010')
    sys.exit()

fasta_file = sys.argv[1]


f = open(fasta_file,"r")
m = 1; # an index for each base
while 1:
    line = f.readline()
    if not line:
        break
    line = line.rstrip()
    s = line.split("\t")
    ch,st,en = re.split('[:-]',s[0]) #st is the start position, whichi is 0-based, inherited from the bed file. 
    myseq = s[1].upper()
    for i in range(len(myseq)-2):
        tri = myseq[i:i+3]
        coor = int(st)+2+i # 1-based coordinates.
        print(ch,coor,coor,tri,m,sep="\t")
        m = m+1

f.close()


