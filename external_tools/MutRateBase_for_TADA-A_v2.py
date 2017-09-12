# in V2 version, any trinucleotide that has at least an "N" will be assigned a mutation rate 0.
from __future__ import print_function
import sys
import re
import subprocess
import sys, re, getopt, math

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if sys.argv[1] == '--help':
    print('MutRateBase.py [rateFile] [triFile]')
    print('rateFile has the allele-specific mutation rate based on trinucleotide model from one trinucleotide to another, e.g., fordist_1KG_mutation_rate_table.txt\n')
    print('triFile has the trinulceotide list generate from a bed file, triFile is the output from tri_extract_changed_by_YL.py')
    sys.exit()

# deal with the options
rateFile = sys.argv[1]
triFile = sys.argv[2]


# process
rateDict = dict()
base = ["A", "T", "C", "G"]
for m in base:
    for n in base:
        for i in base:
            rateDict[m+n+i] = []
            
# mutation rate
if re.search(r'\.gz$',rateFile):
    f = gzip.open(rateFile,'rt')
else:
    f = open(rateFile,'r')
while 1:
    l = f.readline().rstrip()
    if not l:
        break
    if l=="":
        continue
    s = l.split(" ")
    rateDict[s[0]].append("\t".join(s[1:]))

    
f.close()

# trinucleotide list
if re.search(r'\.gz$',triFile):
    f = gzip.open(triFile,'rt')
else:
    f = open(triFile,'r')

while 1:
    l = f.readline().rstrip()
    if not l:
        break
    if l=="":
        continue
    s = l.split("\t")
    # if there is N in the trinucleotide, will report mutation rate as 0, should print 3 lines, each line has a mutation rate as 0. 
    if re.search(r'N',s[3],re.I):
        print(s[0],s[1],s[2],s[3],"NNN", 0, s[4], sep="\t")
        print(s[0],s[1],s[2],s[3],"NNN", 0, s[4], sep="\t")
        print(s[0],s[1],s[2],s[3],"NNN", 0, s[4], sep="\t")
        continue
    for i in rateDict[s[3].upper()]:
        # only report the alternative allele
        i = i.split("\t")
        ref = list(s[3])[1]
        alt = list(i[0])[1]
        # the start position is changed to 1-based
        #print(s[0],s[1],s[2],ref,alt,i[2],sep="\t")
        print(s[0],int(s[1])-1,s[2],ref,alt,i[1],sep="\t")
    
f.close()






