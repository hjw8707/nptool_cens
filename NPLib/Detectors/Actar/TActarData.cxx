/*****************************************************************************
 * Copyright (C) 2009-2017   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: morfouace@ganil.fr                        *
 *                                                                           *
 * Creation Date  : September 2017                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Actar Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "TActarData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

ClassImp(TActarData)


//////////////////////////////////////////////////////////////////////
TActarData::TActarData() {
}



//////////////////////////////////////////////////////////////////////
TActarData::~TActarData() {
}



//////////////////////////////////////////////////////////////////////
void TActarData::Clear() {
  // Charge
    fActar_PadNumber.clear();
    fActar_PadX.clear();
    fActar_PadY.clear();
    fActar_PadZ.clear();
    fActar_PadCharge.clear();

    fActar_RebinningMap.clear();

    fSilicon_Energy.clear();
    fSilicon_Time.clear();
    fSilicon_DetectorNumber.clear();

    fCsI_Energy.clear();
    fCsI_CrystalNumber.clear();
}



//////////////////////////////////////////////////////////////////////
void TActarData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TActarData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Charge
  size_t mysize = fActar_PadNumber.size();
  cout << "Actar_Mult: " << mysize << endl;

  for (size_t i = 0 ; i < mysize ; i++){
    cout << "Pad Number: " << fActar_PadNumber[i]
         << "Charge: " << fActar_PadCharge[i]
      << " X: " << fActar_PadX[i]
      << " Y: " << fActar_PadY[i]
        << " Z: " << fActar_PadZ[i];
  }

}
//////////////////////////////////////////////////////////////////////
void TActarData::RebinData(){
    fActar_PadNumber.clear();
    fActar_PadX.clear();
    fActar_PadY.clear();
    fActar_PadZ.clear();
    fActar_PadCharge.clear();
    fZ_trigger = 256;
    int min_offset = 1000;
    for(auto it = fActar_RebinningMap.begin(); it != fActar_RebinningMap.end(); ++it)
    {
      int x, y, z;
      int rest_z;
      z = it->first/1000000;
      rest_z = it->first-z*1000000;
      y = rest_z/1000;
      x = rest_z-y*1000;
      if((x > 1 && x < 126 && y == 39) || (x>1 && x < 126 && y == 68)) //58Ni experiment
      {
        int this_offset = fabs(z-256);
        if(this_offset < min_offset){
          fZ_trigger = z;
          min_offset = this_offset;
        }
      }
      //cout << "ID : " << it->first << " x = " << x << " y =  " << y << " z = " << z << " q = " << it->second <<  endl;
      //if(it->second > 100)
      fActar_PadNumber.push_back(x*128+y); fActar_PadX.push_back(x); fActar_PadY.push_back(y); fActar_PadZ.push_back(z); fActar_PadCharge.push_back(it->second);
    }
}
