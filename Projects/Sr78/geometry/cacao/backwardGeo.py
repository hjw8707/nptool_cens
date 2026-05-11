#!/usr/bin/env python
def PrintBlock(x, y, z, a, b, c):
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("CACAO")
    print(" Pos= %6.1f %6.1f %6.1f mm" % (x, y, z))
    print(" Ang= %6.1f %6.1f %6.1f deg" % (a, b, c))
    print(" Dim= %6.1f %6.1f %6.1f mm" % (15, 15, 50))    
    print(" ShieldThickness= 0 mm")
    
(a, b, c) = (0, 90, 90)
slope = -95/160.
offset = []
for i in range(2):
    offset.append(0)
for i in range(1,4):
    offset.append(slope*i*16)
for i in range(8):
    offset.append(slope*3*16)
    
for iring in range(0,10):
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("%% for ring %d" % iring)
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    xDef = 100
    xOff = offset[iring]
    z = -80.0 - iring * 16

    n_block = 7
    if iring > 1: n_block = 6
    if iring > 3: n_block = 5
    for iblock in range(-n_block+1,n_block):
        y = iblock * 16.
        x = ((xDef+xOff)**2 - y**2)**0.5
        PrintBlock(x,y,z,a,b,-c)
        PrintBlock(-x,y,z,a,b,c)

