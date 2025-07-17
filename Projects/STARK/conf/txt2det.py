#!/usr/bin/python3
import sys
import re
import argparse


parser = argparse.ArgumentParser(
    prog = 'txt2det',
    description='STARK configuration table to NPTool detector file')

parser.add_argument('input', help='STARK configuration')
parser.add_argument('-b', '--noBB10', help='no BB10', action='store_true')
parser.add_argument('-q', '--noQQQ5', help='no QQQ5', action='store_true')

args = parser.parse_args()

f= open(args.input, 'r')
lines = f.readlines()
f.close()

cols = [ 'Type', 'Rho', 'Phi', 'Z', 'Rev', 'Flip', 'Beta', 'Group' ]

stark = []
for line in lines:
    line = line.strip()
    if line[0] == '#': continue
    words = line.split()
    starkEl = {}
    for idx, word in enumerate(words):
        starkEl[cols[idx]] = word
    stark.append(starkEl)

for el in stark:
    if args.noBB10 and el['Type'] == 'BB10': continue
    if args.noQQQ5 and el['Type'] == 'QQQ5': continue

    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("STARK")
    print(" Type= %s" % el['Type'])
    print(" Rho= %10s mm" % el['Rho'])
    print(" Phi= %10s deg" % el['Phi'])
    if el['Type'] == 'QQQ5':
        print(" Beta= %10s deg" % el['Beta'])
    print(" Z= %10s mm" % el['Z'])
    print(" Rev= %10s" % el['Rev'])
    print(" Flip= %10s" % el['Flip'])
    print(" Group= %10s" % el['Group'])
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
