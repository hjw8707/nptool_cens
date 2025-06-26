/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Léo Plagnol  contact address: leo.plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : January 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Coaxial_Germanium Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TCoaxial_GermaniumData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TCoaxial_GermaniumData)


//////////////////////////////////////////////////////////////////////
TCoaxial_GermaniumData::TCoaxial_GermaniumData() {
}



//////////////////////////////////////////////////////////////////////
TCoaxial_GermaniumData::~TCoaxial_GermaniumData() {
}



//////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumData::Clear() {
  // Energy
  fCoaxial_Germanium_E_DetectorNbr.clear();
  fCoaxial_Germanium_Energy.clear();
  // Time
  fCoaxial_Germanium_T_DetectorNbr.clear();
  fCoaxial_Germanium_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TCoaxial_GermaniumData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fCoaxial_Germanium_E_DetectorNbr.size();
  cout << "Coaxial_Germanium_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fCoaxial_Germanium_E_DetectorNbr[i]
         << " Energy: " << fCoaxial_Germanium_Energy[i];
  }
  
  // Time
  mysize = fCoaxial_Germanium_T_DetectorNbr.size();
  cout << "Coaxial_Germanium_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fCoaxial_Germanium_T_DetectorNbr[i]
         << " Time: " << fCoaxial_Germanium_Time[i];
  }
}
