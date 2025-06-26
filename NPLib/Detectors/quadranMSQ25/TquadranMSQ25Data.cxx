/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: leo.plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : January 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold quadranMSQ25 Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TquadranMSQ25Data.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TquadranMSQ25Data)


//////////////////////////////////////////////////////////////////////
TquadranMSQ25Data::TquadranMSQ25Data() {
}



//////////////////////////////////////////////////////////////////////
TquadranMSQ25Data::~TquadranMSQ25Data() {
}



//////////////////////////////////////////////////////////////////////
void TquadranMSQ25Data::Clear() {
  // Energy
  fquadranMSQ25_E_DetectorNbr.clear();
  fquadranMSQ25_Energy.clear();
  // Time
  fquadranMSQ25_T_DetectorNbr.clear();
  fquadranMSQ25_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TquadranMSQ25Data::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TquadranMSQ25Data::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fquadranMSQ25_E_DetectorNbr.size();
  cout << "quadranMSQ25_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fquadranMSQ25_E_DetectorNbr[i]
         << " Energy: " << fquadranMSQ25_Energy[i];
  }
  
  // Time
  mysize = fquadranMSQ25_T_DetectorNbr.size();
  cout << "quadranMSQ25_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fquadranMSQ25_T_DetectorNbr[i]
         << " Time: " << fquadranMSQ25_Time[i];
  }
}
