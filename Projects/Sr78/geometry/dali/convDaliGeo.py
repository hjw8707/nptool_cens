#!/usr/bin/python3


import sys
import re

zoffset = 0

f= open(sys.argv[1], 'r')
lines = f.readlines()
f.close()

daliEls = [ 'X', 'Y', 'Z', 'a', 'b', 'c', 't', 'type' ]

dali = []
for line in lines:
    line = line.strip()
    words = line.split()
    daliEl = {}
    for idx, word in enumerate(words):
        daliEl[daliEls[idx]] = word
    dali.append(daliEl)

for el in dali:
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("Dali2")
    z = (float(el['Z']) + zoffset) if float(el['Z']) > 0 else (float(el['Z']) - zoffset) 
    print(" Pos = %10s %10s %10.2f cm" % (el['X'], el['Y'], z))
    print(" Ang = %10s %10s %10s deg" % (el['a'], el['b'], el['c']))
    print(" Type = %3s" % el['type'])

print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")    

#for el in dali:
#    x = float(el['X'])
#    y = float(el['Y'])         
#    z = (float(el['Z']) + zoffset) if float(el['Z']) > 0 else (float(el['Z']) - zoffset)     
#    print("%10.2f" % ((x*x + y*y + z*z)**0.5))
