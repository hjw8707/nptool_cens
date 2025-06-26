/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: P. Morfouace contact address: pierre.morfouace@cea.fr    *
 *                                                                           *
 * Creation Date  : May 2020                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold PISTA Raw data                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TPISTAData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TPISTAData)


//////////////////////////////////////////////////////////////////////
TPISTAData::TPISTAData() {
}



//////////////////////////////////////////////////////////////////////
TPISTAData::~TPISTAData() {
}



//////////////////////////////////////////////////////////////////////
void TPISTAData::Clear() {
  // DE //
  fPISTA_DE_DetectorNbr.clear();
  fPISTA_DE_StripNbr.clear();
  fPISTA_DE_StripEnergy.clear();
  fPISTA_DE_BackEnergy.clear();
  fPISTA_DE_StripTime.clear();
  fPISTA_DE_BackTime.clear();
  fPISTA_DE_BackDetector.clear();
 
  // E //
  fPISTA_E_DetectorNbr.clear();
  fPISTA_E_StripNbr.clear();
  fPISTA_E_StripEnergy.clear();
  fPISTA_E_BackEnergy.clear();
  fPISTA_E_StripTime.clear();
  fPISTA_E_BackTime.clear();
  fPISTA_E_BackDetector.clear();
  
}



//////////////////////////////////////////////////////////////////////
void TPISTAData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TPISTAData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // DE
  size_t mysize = fPISTA_DE_DetectorNbr.size();
  cout << "PISTA_DE_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DE-DetNbr: " << fPISTA_DE_DetectorNbr[i] << endl;
    cout << " DE-Strip: " << fPISTA_DE_StripNbr[i] << endl;
    cout << " DE-Energy: " << fPISTA_DE_StripEnergy[i] << endl;
  }
 
  // E
  mysize = fPISTA_E_DetectorNbr.size();
  cout << "PISTA_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "E-DetNbr: " << fPISTA_E_DetectorNbr[i] << endl;
    cout << " E-Strip: " << fPISTA_E_StripNbr[i] << endl;
    cout << " E-Energy: " << fPISTA_E_BackEnergy[i] << endl;
  }
}
