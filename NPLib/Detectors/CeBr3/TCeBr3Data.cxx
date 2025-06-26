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
 *  This class hold CeBr3 Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TCeBr3Data.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TCeBr3Data)


//////////////////////////////////////////////////////////////////////
TCeBr3Data::TCeBr3Data() {
}



//////////////////////////////////////////////////////////////////////
TCeBr3Data::~TCeBr3Data() {
}



//////////////////////////////////////////////////////////////////////
void TCeBr3Data::Clear() {
  // Energy
  fCeBr3_E_DetectorNbr.clear();
  fCeBr3_Energy.clear();
  // Time
  fCeBr3_T_DetectorNbr.clear();
  fCeBr3_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TCeBr3Data::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TCeBr3Data::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fCeBr3_E_DetectorNbr.size();
  cout << "CeBr3_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fCeBr3_E_DetectorNbr[i]
         << " Energy: " << fCeBr3_Energy[i];
  }
  
  // Time
  mysize = fCeBr3_T_DetectorNbr.size();
  cout << "CeBr3_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fCeBr3_T_DetectorNbr[i]
         << " Time: " << fCeBr3_Time[i];
  }
}
