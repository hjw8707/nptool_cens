/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: jungwoo  contact address: phyjics@gmail.com                        *
 *                                                                           *
 * Creation Date  : May 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ATOMX Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TATOMXData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TATOMXData)

//////////////////////////////////////////////////////////////////////
TATOMXData::TATOMXData()
{
}

//////////////////////////////////////////////////////////////////////
TATOMXData::~TATOMXData()
{
}

//////////////////////////////////////////////////////////////////////
void TATOMXData::Clear()
{
    fEnergyLoss.clear();
    fTime.clear();
    fPosition.clear();
}

//////////////////////////////////////////////////////////////////////
void TATOMXData::Dump(int max_points) const
{
    auto n = fEnergyLoss.size();
    cout << "ATOMX event containing " << n << " points [i](e,t|x,y,z)" << endl;
    if (n>max_points)
        n = max_points;
    for (auto i=0; i<n; ++i)
    {
        cout << "  [" << i << "] (" << fEnergyLoss[i] << ", "
            << fTime [i] << " | "
            << fPosition[i].x() << ", "
            << fPosition[i].y() << ", "
            << fPosition[i].z() << ")" << endl;
    }
}
//////////////////////////////////////////////////////////////////////
void TATOMXData::Dump() const
{
    auto n = fEnergyLoss.size();
    cout << "ATOMX event containing " << n << " points [i](e,t|x,y,z)" << endl;
    for (auto i=0; i<n; ++i)
    {
        cout << "  [" << i << "] (" << fEnergyLoss[i] << ", "
            << fTime [i] << " | "
            << fPosition[i].x() << ", "
            << fPosition[i].y() << ", "
            << fPosition[i].z() << ")" << endl;
    }
}
