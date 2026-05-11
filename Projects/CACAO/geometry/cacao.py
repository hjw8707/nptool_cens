#!/usr/bin/env python
import sys

def PrintBlock(x, y, z, a, b, c):
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("CACAO")
    print(" Pos= %6.1f %6.1f %6.1f mm" % (x, y, z))
    print(" Ang= %6.1f %6.1f %6.1f deg" % (a, b, c))

cacao_xy = 36 # 36 mm for each module
cacao_z_center = 50/2. + 1.5 # from the module bottom to the center of the detector
z_start = -48 - 36
n_zlayer = 7
n_cacao = 0

###############################################
# bottom wall: 30 modules (5 modules per row)
#
def BottomWall(n_xlayers = 5, distance = 105.):
    global n_cacao
    (a, b, c) = (0, -90, 0)
    y = -abs(distance)
    for i in range(n_zlayer):
        z = z_start + cacao_xy*i
        for iblock in range(n_xlayers): # n layers
            x = (iblock - (n_xlayers-1)/2.) * cacao_xy
            PrintBlock(x,y,z,a,b,c)
            n_cacao = n_cacao + 1

################################################
# side wall: 14 x 12 blocks => 7 x 6 modules (2x2) layout
# left, right wall: 42 x 2 blocks = 84 blocks
#
def SideWall(n_ylayers = 7, distance = 130):
    global n_cacao
    (a, b, c) = (0, -90, -90)
    x = abs(distance)
    for i in range(n_zlayer):
        z = z_start + cacao_xy*i
        for iblock in range(n_ylayers): # n layers
            y = (iblock - (n_ylayers-1)/2.) * cacao_xy
            PrintBlock(x,y,z,a,-b,-c)
            PrintBlock(-x,y,z,a,-b,c)
            n_cacao = n_cacao + 2

########################################################
# Top Wall
# top wall: 10 modules (2 modules per row)
def TopWall(n_xlayers = 5, distance = 105):
    global n_cacao
    (a, b, c) = (0, 90, 0)
    y = abs(distance)
    for i in [0, 4, 5, 6]:
        z = z_start + cacao_xy*i
        if i == 1 or i == 4:
            for iblock in [0, n_xlayers-1]: # n layers
                x = (iblock - (n_xlayers-1)/2.) * cacao_xy
                PrintBlock(x,y,z,a,b,c)
                n_cacao = n_cacao + 1
            continue
        for iblock in range(n_xlayers): # n layers
            x = (iblock - (n_xlayers-1)/2.) * cacao_xy
            PrintBlock(x,y,z,a,b,c)
            n_cacao = n_cacao + 1

################################################
# backward blocks (left, right)
#
def Rings(n_rings = 7):
    global n_cacao
    (a, b, c) = (0, 90, 90)
    slope = -95/160.
    offset = []
    for i in range(2):
        offset.append(0)
    for i in range(1,4):
        offset.append(slope*i*16)
    for i in range(2):
        offset.append(slope*3*16)

    for iring in range(0,n_rings):
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("%% for ring %d" % iring)
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

        xDef = 110
        xOff = offset[iring]
        z = -100.0 - iring * 16

        n_block = 7
        if iring > 1: n_block = 6
        if iring > 3: n_block = 5
        for iblock in range(-n_block+1,n_block):
            y = iblock * 16.
            x = ((xDef+xOff)**2 - y**2)**0.5
            PrintBlock(x,y,z,a,b,-c)
            PrintBlock(-x,y,z,a,b,c)
            n_cacao = n_cacao + 2

def Rings2x2(n_rings = 4):
    global n_cacao
    (a, b, c) = (0, 90, 90)
    slope = -95/160.
    offset = []
    #for i in range(1):   offset.append(0)
    #for i in range(1,2): offset.append(slope*i*32)
    #for i in range(2):   offset.append(slope*3*16)
    offset.append(0)
    offset.append(-20)
    offset.append(-20)
    #offset.append(-12)

    for iring in range(0,len(offset)):
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("%% for ring %d" % iring)
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

        xDef = 110
        xOff = offset[iring]
        z = -100.0 - iring * 32

        n_block = 4
        #if iring > 0: n_block = 3
        for iblock in range(-n_block+1,n_block):
            y = iblock * 32.
            if (xDef+xOff)**2 - 1.05*y**2 > 0: x = ((xDef+xOff)**2 - 1.05*y**2)**0.5
            else: x = 30
            PrintBlock(x,y+8,z,a,b,-c)
            PrintBlock(x,y-8,z,a,b,-c)
            PrintBlock(x,y+8,z-16,a,b,-c)
            PrintBlock(x,y-8,z-16,a,b,-c)
            n_cacao = n_cacao + 4

        for iblock in range(-n_block+1,n_block):
            y = iblock * 32.
            if (xDef+xOff)**2 - 1.05*y**2 > 0: x = ((xDef+xOff)**2 - 1.05*y**2)**0.5
            else: x = 30
            PrintBlock(-x,y+8,z,a,b,c)
            PrintBlock(-x,y-8,z,a,b,c)
            PrintBlock(-x,y+8,z-16,a,b,c)
            PrintBlock(-x,y-8,z-16,a,b,c)
            n_cacao = n_cacao + 4

#######################
# total = 456 blocks
#######################

if __name__=="__main__":

    # dist 변수 설정: 첫 번째 매개변수로부터 받거나, 없으면 105로 기본값 설정
    if len(sys.argv) > 1:
        try:
            dist = float(sys.argv[1])
        except ValueError:
            print("경고: 첫 번째 인자가 숫자가 아닙니다. 기본값 105를 사용합니다.")
            dist = 105
    else:
        dist = 105
    BottomWall(distance=105)
    print("%% Bottom total = %d" % n_cacao)
    SideWall(distance=116)
    print("%% Bottom + Side total = %d" % n_cacao)
    TopWall(distance=105)
    print("%% Bottom + Side + Top total = %d" % n_cacao)
    #Rings()
    #Rings2x2()
    print("%% total = %d" % n_cacao)
