/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: hjw8707@gmail.com                        *
 *                                                                           *
 * Creation Date  : 7월 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Plunger Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TPlungerData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TPlungerData)


//////////////////////////////////////////////////////////////////////
TPlungerData::TPlungerData() {
}



//////////////////////////////////////////////////////////////////////
TPlungerData::~TPlungerData() {
}



//////////////////////////////////////////////////////////////////////
void TPlungerData::Clear() {
  // Energy
  fPlunger_E_DetectorNbr.clear();
  fPlunger_Energy.clear();
  // Time
  fPlunger_T_DetectorNbr.clear();
  fPlunger_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TPlungerData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TPlungerData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fPlunger_E_DetectorNbr.size();
  cout << "Plunger_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fPlunger_E_DetectorNbr[i]
         << " Energy: " << fPlunger_Energy[i];
  }
  
  // Time
  mysize = fPlunger_T_DetectorNbr.size();
  cout << "Plunger_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fPlunger_T_DetectorNbr[i]
         << " Time: " << fPlunger_Time[i];
  }
}
