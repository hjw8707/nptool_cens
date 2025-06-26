/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Pohl  contact address: thomas.pohl@riken.jp                        *
 *                                                                           *
 * Creation Date  : February 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TOGAXSI_SI Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TTOGAXSI_SIData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TTOGAXSI_SIData)


//////////////////////////////////////////////////////////////////////
TTOGAXSI_SIData::TTOGAXSI_SIData() {
}



//////////////////////////////////////////////////////////////////////
TTOGAXSI_SIData::~TTOGAXSI_SIData() {
}



//////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIData::Clear() {
  // Energy X
  fInnerX_E_DetectorNbr.clear();
  fInnerX_E_StripNbr.clear();
  fInnerX_E_Energy.clear();
  // Energy Z
  fInnerZ_E_DetectorNbr.clear();
  fInnerZ_E_StripNbr.clear();
  fInnerZ_E_Energy.clear();

  // Energy X
  fOuterX_E_DetectorNbr.clear();
  fOuterX_E_StripNbr.clear();
  fOuterX_E_Energy.clear();
  // Energy Z
  fOuterZ_E_DetectorNbr.clear();
  fOuterZ_E_StripNbr.clear();
  fOuterZ_E_Energy.clear();

  // Energy X
  fClusterInner_E_DetectorNbr.clear();
  fClusterInner_E_StripNbr.clear();
  fClusterInner_E_Energy.clear();

  // Energy X1
  fClusterX1_E_DetectorNbr.clear();
  fClusterX1_E_StripNbr.clear();
  fClusterX1_E_Energy.clear();

  // Energy Y1
  fClusterY1_E_DetectorNbr.clear();
  fClusterY1_E_StripNbr.clear();
  fClusterY1_E_Energy.clear();

  // Energy X2
  fClusterX2_E_DetectorNbr.clear();
  fClusterX2_E_StripNbr.clear();
  fClusterX2_E_Energy.clear();

  // Energy Y2
  fClusterY2_E_DetectorNbr.clear();
  fClusterY2_E_StripNbr.clear();
  fClusterY2_E_Energy.clear();

}



//////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TTOGAXSI_SIData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fInnerX_E_DetectorNbr.size();
  cout << "TOGAXSI_SI_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fInnerX_E_DetectorNbr[i]
         << " Energy: " << fInnerX_E_Energy[i];
  }
  
}
