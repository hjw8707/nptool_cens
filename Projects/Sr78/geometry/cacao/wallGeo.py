#!/usr/bin/env python
def PrintBlock(x, y, z, a, b, c):
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("CACAO")
    print(" Pos= %6.1f %6.1f %6.1f mm" % (x, y, z))
    print(" Ang= %6.1f %6.1f %6.1f deg" % (a, b, c))
    print(" Dim= %6.1f %6.1f %6.1f mm" % (15, 15, 50))    
    print(" ShieldThickness= 0 mm")

z_start = -48 - 32
n_zlayer = 12
    
###############################################
# bottom wall: 120 blocks
# 
n_xlayers = 10
(a, b, c) = (0, 90, 0)
y = -105.
for i in range(n_zlayer):
    z = z_start + 16.*i
    for iblock in range(n_xlayers): # n layers
        x = (iblock - (n_xlayers-1)/2.) * 16.
        PrintBlock(x,y,z,a,b,c)

################################################
# left, right wall: 168 x 2 blocks = 336 blocks
#
n_ylayers = 14
(a, b, c) = (0, 90, 90)
x = 105.
for i in range(n_zlayer):
    z = z_start + 16.*i
    for iblock in range(n_ylayers): # n layers
        y = (iblock - (n_ylayers-1)/2.) * 16.
        PrintBlock(x,y,z,a,b,-c)
        PrintBlock(-x,y,z,a,b,c)        

#######################
# total = 456 blocks
#######################