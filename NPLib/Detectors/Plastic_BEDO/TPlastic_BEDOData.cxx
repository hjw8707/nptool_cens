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
 *  This class hold Plastic_BEDO Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TPlastic_BEDOData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TPlastic_BEDOData)


//////////////////////////////////////////////////////////////////////
TPlastic_BEDOData::TPlastic_BEDOData() {
}



//////////////////////////////////////////////////////////////////////
TPlastic_BEDOData::~TPlastic_BEDOData() {
}



//////////////////////////////////////////////////////////////////////
void TPlastic_BEDOData::Clear() {
  // Energy
  fPlastic_BEDO_E_DetectorNbr.clear();
  fPlastic_BEDO_Energy.clear();
  // Time
  fPlastic_BEDO_T_DetectorNbr.clear();
  fPlastic_BEDO_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TPlastic_BEDOData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TPlastic_BEDOData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fPlastic_BEDO_E_DetectorNbr.size();
  cout << "Plastic_BEDO_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fPlastic_BEDO_E_DetectorNbr[i]
         << " Energy: " << fPlastic_BEDO_Energy[i];
  }
  
  // Time
  mysize = fPlastic_BEDO_T_DetectorNbr.size();
  cout << "Plastic_BEDO_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fPlastic_BEDO_T_DetectorNbr[i]
         << " Time: " << fPlastic_BEDO_Time[i];
  }
}
