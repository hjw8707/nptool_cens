#!/usr/bin/python3
import sys
import re
import argparse

parser = argparse.ArgumentParser(description='convCsIGeo')

parser.add_argument('-d', '--doffset', type = float, default = 0, help = 'distance offset from the target')
parser.add_argument('-z', '--zoffset', type = float, default = 0, help = 'Z(beam direction) offset')
parser.add_argument('-t', '--target', type = str, default = 'lh2', help = 'lh2 or ch2')
parser.add_argument('geometry', type = str, default = 'csi_geometry.txt', help = 'input geometry text file')

args = parser.parse_args()

f= open(args.geometry, 'r')
lines = f.readlines()
f.close()

daliEls = [ 'X', 'Y', 'Z', 'a', 'b', 'c' ]

dali = []
for line in lines:
    line = line.strip()
    if line[0] == '#': continue
    words = line.split()
    daliEl = {}
    for idx, word in enumerate(words):
        daliEl[daliEls[idx]] = float(word)
    dali.append(daliEl)

for el in dali:
    ##################################################
    # distance offset
    if   float(el['a']) == 90 : el['Y'] -= args.doffset
    elif float(el['b']) == 90 : el['X'] += args.doffset
    elif float(el['b']) == 270: el['X'] -= args.doffset
    ##################################################

    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("CsI")
    print(" Pos= %6.1f %6.1f %6.1f mm" % (el['X'], el['Y'], el['Z']))
    print(" Ang= %6.1f %6.1f %6.1f deg" % (el['a'], el['b'], el['c']))
    print(" Shape= Square")
    print(" Thickness= 50 mm")
    print(" FaceFront= 15 mm")
    print(" FaceBack= 15 mm")
    print(" Scintillator= CsI_Scintillator")
    print(" LeadThickness= 0 mm")

print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")    

#for el in dali:
#    x = float(el['X'])
#    y = float(el['Y'])         
#    z = (float(el['Z']) + zoffset) if float(el['Z']) > 0 else (float(el['Z']) - zoffset)     
#    print("%10.2f" % ((x*x + y*y + z*z)**0.5))
