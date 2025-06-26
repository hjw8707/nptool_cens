/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : February 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Vendeta Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TVendetaData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TVendetaData)


//////////////////////////////////////////////////////////////////////
TVendetaData::TVendetaData() {
}



//////////////////////////////////////////////////////////////////////
TVendetaData::~TVendetaData() {
}



//////////////////////////////////////////////////////////////////////
void TVendetaData::Clear() {
  
  fVendeta_LG_DetectorNbr.clear();
  fVendeta_LG_Q1.clear();
  fVendeta_LG_Q2.clear();
  fVendeta_LG_Time.clear();
  fVendeta_LG_Qmax.clear();
  fVendeta_LG_IsSat.clear();
  fVendeta_LG_IsGamma.clear();
 
  fVendeta_HG_DetectorNbr.clear();
  fVendeta_HG_Q1.clear();
  fVendeta_HG_Q2.clear();
  fVendeta_HG_Time.clear();
  fVendeta_HG_Qmax.clear();
  fVendeta_HG_IsSat.clear();
  fVendeta_HG_IsGamma.clear();
  
}



//////////////////////////////////////////////////////////////////////
void TVendetaData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TVendetaData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  size_t mysize = fVendeta_HG_DetectorNbr.size();
  cout << "Vendeta_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fVendeta_HG_DetectorNbr[i]
         << " Q1: "    << fVendeta_HG_Q1[i]
         << " Q2: "    << fVendeta_HG_Q2[i]
         << " Time: "  << fVendeta_HG_Time[i];
  }
}
