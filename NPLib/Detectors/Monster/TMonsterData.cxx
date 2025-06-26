/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Emile Cantacuzène  contact address: emile.cantacuzene@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : August 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Monster Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TMonsterData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TMonsterData)


//////////////////////////////////////////////////////////////////////
TMonsterData::TMonsterData() {
}



//////////////////////////////////////////////////////////////////////
TMonsterData::~TMonsterData() {
}



//////////////////////////////////////////////////////////////////////
void TMonsterData::Clear() {
  // Energy
  fMonster_E_DetectorNbr.clear();
  fMonster_Energy.clear();
  // Time
  fMonster_T_DetectorNbr.clear();
  fMonster_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TMonsterData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TMonsterData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fMonster_E_DetectorNbr.size();
  cout << "Monster_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fMonster_E_DetectorNbr[i]
         << " Energy: " << fMonster_Energy[i];
  }
  
  // Time
  mysize = fMonster_T_DetectorNbr.size();
  cout << "Monster_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fMonster_T_DetectorNbr[i]
         << " Time: " << fMonster_Time[i];
  }
}
