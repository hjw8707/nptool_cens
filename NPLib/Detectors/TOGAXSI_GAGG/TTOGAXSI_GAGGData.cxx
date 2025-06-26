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
 *  This class hold TOGAXSI_GAGG Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TTOGAXSI_GAGGData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TTOGAXSI_GAGGData)


//////////////////////////////////////////////////////////////////////
TTOGAXSI_GAGGData::TTOGAXSI_GAGGData() {
}



//////////////////////////////////////////////////////////////////////
TTOGAXSI_GAGGData::~TTOGAXSI_GAGGData() {
}



//////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGData::Clear() {
  // Energy
  fRecoil_E_DetectorNbr.clear();
  fRecoil_Energy.clear();
  // Time
  fRecoil_T_DetectorNbr.clear();
  fRecoil_Time.clear();

  // Energy
  fCluster_E_DetectorNbr.clear();
  fCluster_Energy.clear();
  // Time
  fCluster_T_DetectorNbr.clear();
  fCluster_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGData::Dump() const {
  // This method is very useful for debuging and worth the dev.

  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TTOGAXSI_GAGGData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fRecoil_E_DetectorNbr.size();
  cout << "Recoil GAGG_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fRecoil_E_DetectorNbr[i]
         << " Energy: " << fRecoil_Energy[i];
  }
  
  // Time
  mysize = fRecoil_T_DetectorNbr.size();
  cout << "Recoil GAGG_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fRecoil_T_DetectorNbr[i]
         << " Time: " << fRecoil_Time[i];
  }

}
