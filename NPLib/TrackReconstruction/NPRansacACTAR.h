#ifndef __RANSACACTAR__
#define __RANSACACTAR__
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Damien THISSE contact address: damien.thisse@cea.fr    *
 *                                                                           *
 * Creation Date   : October 2024                                            *
 * Last update     : October 2024                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class deal with finding all the track event by event                *
 *****************************************************************************/

//C++ Header
#include <stdio.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <limits>
#include <string>
#include <algorithm>

//NPL
#include "NPTrack.h"
#include "NPLinearRansac3D.h"

using namespace NPL;

// ROOT Headers
#include <TH2F.h>
#include <TMath.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TVector.h>

using namespace std;

namespace NPL{
    
    class RansacACTAR
    {
    public:
        RansacACTAR();
        ~RansacACTAR() {};
        
    public:
        void ReadParameterValue(string filename);
        vector<NPL::LinearCluster3D> SimpleRansac(vector<double> & vX, vector<double> & vY, vector<double> & vZ, vector<double> & vQ);
        
    private:
        TRandom3* Rand;
        
    private:
        float fRANSACThreshold;
        float fRANSACPointThreshold;
        float fRANSACChargeThreshold;
        float fRANSACDistance;
        float fRANSACMaxIteration;
        float fMAXBarycenterDistance;
        float fAngleMax;
        int fNumberOfTracksMax;
        int fOriginalCloudSize;
        double fTotalCharge;
        int fNumberOfPadsX;
        int fNumberOfPadsY;
        
    };
}
#endif
