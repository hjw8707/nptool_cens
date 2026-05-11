#!/usr/bin/python
def PrintBlock(x, y, z, a, b, c):
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("CsI")
    print(" Pos= %6.1f %6.1f %6.1f mm" % (x, y, z))
    print(" Ang= %6.1f %6.1f %6.1f deg" % (a, b, c))
    print(" Shape= Square")
    print(" Thickness= 50 mm")
    print(" FaceFront= 15 mm")
    print(" FaceBack= 15 mm")
    print(" Scintillator= CsI_Scintillator")
    print(" LeadThickness= 0 mm")

    
(a, b, c) = (0, 60, 0)
slope = 95/160.+ 0.25
offset = []
for i in range(2):
    offset.append(0)
for i in range(1,9):
    offset.append(slope*i*16)
          
for iring in range(0,8):
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("%% for ring %d" % iring)
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    xDef = 100
    xOff = offset[iring]
    z = 80.0 + iring * 16 * 1.2 + xOff/3**0.5

    base = 8
    n_block = base
    if iring > 1: n_block = base + 1
    if iring > 3: n_block = base + 2
    if iring > 4: n_block = base + 3
    if iring > 5: n_block = base + 4
    if iring > 6: n_block = base + 5
    for iblock in range(-n_block+1,n_block):
        y = iblock * 16.
        x = ((xDef+xOff)**2 - (y/1.13)**2)**0.5
        if abs(x) < 30: x = 25*x/abs(x)
        PrintBlock(x,y,z,a,b,c)
        PrintBlock(-x,y,z,a,-b,c)

