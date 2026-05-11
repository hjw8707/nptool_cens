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

z_start = -48.
n_zlayer = 8
    
# bottom wall
(a, b, c) = (90, 90, 0)
y = -100.
for i in range(n_zlayer):
    z = z_start + 16.*i
    for iblock in range(-4,5):
        x = iblock * 16.
        PrintBlock(x,y,z,a,b,c)

# left, right wall
(a, b, c) = (0, 90, 0)
x = 100.
for i in range(n_zlayer):
    z = z_start + 16.*i
    for iblock in range(-7,8):
        y = iblock * 16.
        PrintBlock(x,y,z,a,b,c)
        PrintBlock(-x,y,z,a,b,c)        

