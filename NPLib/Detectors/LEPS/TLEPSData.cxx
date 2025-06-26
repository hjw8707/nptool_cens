/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : March 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold LEPS Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TLEPSData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TLEPSData)


//////////////////////////////////////////////////////////////////////
TLEPSData::TLEPSData() {
}



//////////////////////////////////////////////////////////////////////
TLEPSData::~TLEPSData() {
}



//////////////////////////////////////////////////////////////////////
void TLEPSData::Clear() {
  // Energy
  fLEPS_E_DetectorNbr.clear();
  fLEPS_Energy.clear();
  // Time
  fLEPS_T_DetectorNbr.clear();
  fLEPS_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TLEPSData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TLEPSData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fLEPS_E_DetectorNbr.size();
  cout << "LEPS_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fLEPS_E_DetectorNbr[i]
         << " Energy: " << fLEPS_Energy[i];
  }
  
  // Time
  mysize = fLEPS_T_DetectorNbr.size();
  cout << "LEPS_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fLEPS_T_DetectorNbr[i]
         << " Time: " << fLEPS_Time[i];
  }
}
