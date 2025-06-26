/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : December 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold GeLi Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TGeLiData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TGeLiData)


//////////////////////////////////////////////////////////////////////
TGeLiData::TGeLiData() {
}



//////////////////////////////////////////////////////////////////////
TGeLiData::~TGeLiData() {
}



//////////////////////////////////////////////////////////////////////
void TGeLiData::Clear() {
  // Energy
  fGeLi_E_DetectorNbr.clear();
  fGeLi_Energy.clear();
  // Time
  fGeLi_T_DetectorNbr.clear();
  fGeLi_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TGeLiData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TGeLiData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fGeLi_E_DetectorNbr.size();
  cout << "GeLi_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fGeLi_E_DetectorNbr[i]
         << " Energy: " << fGeLi_Energy[i];
  }
  
  // Time
  mysize = fGeLi_T_DetectorNbr.size();
  cout << "GeLi_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fGeLi_T_DetectorNbr[i]
         << " Time: " << fGeLi_Time[i];
  }
}
