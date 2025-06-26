/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : September 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold FissionChamber Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TFissionChamberData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TFissionChamberData)


//////////////////////////////////////////////////////////////////////
TFissionChamberData::TFissionChamberData() {
}



//////////////////////////////////////////////////////////////////////
TFissionChamberData::~TFissionChamberData() {
}



//////////////////////////////////////////////////////////////////////
void TFissionChamberData::Clear() {
  fFC_AnodeNbr.clear();
  fFC_Q1.clear();
  fFC_Q2.clear();
  fFC_Qmax.clear();
  fFC_Time.clear();
  fFC_Time_HF.clear();
  fFC_ToF.clear();
  fFC_isFakeFission.clear();
}



//////////////////////////////////////////////////////////////////////
void TFissionChamberData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TFissionChamberData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  size_t mysize = fFC_Q1.size();
  cout << "FissionChamber_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "AnodeNbr: " << fFC_AnodeNbr[i]
         << " Q1: " << fFC_Q1[i]
         << " Q2: " << fFC_Q2[i]
         << " Qmax: " << fFC_Qmax[i]
         << " Time: " << fFC_Time[i];
  }
}
