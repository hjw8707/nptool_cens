/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: freddy flavigny  contact: flavigny@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  : OCtober 2022                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold NebulaPlus Raw data                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TNebulaPlusData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TNebulaPlusData)


//////////////////////////////////////////////////////////////////////
TNebulaPlusData::TNebulaPlusData() {
}



//////////////////////////////////////////////////////////////////////
TNebulaPlusData::~TNebulaPlusData() {
}



//////////////////////////////////////////////////////////////////////
void TNebulaPlusData::Clear() {
    // UP // 
    fNebulaPlus_u_ID.clear();
    fNebulaPlus_u_Q.clear();
    fNebulaPlus_u_Q4.clear();
    fNebulaPlus_u_T.clear();
    
    // DOWN // 
    fNebulaPlus_d_ID.clear();
    fNebulaPlus_d_Q.clear();
    fNebulaPlus_d_Q4.clear();
    fNebulaPlus_d_T.clear();
}



//////////////////////////////////////////////////////////////////////
void TNebulaPlusData::Dump() const {
/*  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TNebulaPlusData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fNebulaPlus_E_DetectorNbr.size();
  cout << "NebulaPlus_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fNebulaPlus_E_DetectorNbr[i]
         << " Energy: " << fNebulaPlus_Energy[i];
  }
  
  // Time
  mysize = fNebulaPlus_T_DetectorNbr.size();
  cout << "NebulaPlus_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fNebulaPlus_T_DetectorNbr[i]
         << " Time: " << fNebulaPlus_Time[i];
  }
  */
}
