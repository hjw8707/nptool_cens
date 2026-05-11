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

    
(a, b, c) = (0, 90, 0)
slope = 95/160.
offset = []
for i in range(2):
    offset.append(0)
for i in range(1,9):
    offset.append(slope*i*16)
          
for iring in range(0,10):
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("%% for ring %d" % iring)
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    xDef = 100
    xOff = offset[iring]
    z = 80.0 + iring * 16

    n_block = 7
    if iring > 2: n_block = 8
    if iring > 4: n_block = 9
    for iblock in range(-n_block+1,n_block):
        y = iblock * 16.
        x = ((xDef+xOff)**2 - y**2)**0.5
        PrintBlock(x,y,z,a,b,c)
        PrintBlock(-x,y,z,a,b,c)

