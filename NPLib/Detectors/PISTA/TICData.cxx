/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr
 *                                                                           *
 * Creation Date  : Oct 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold IC Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TICData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TICData)


//////////////////////////////////////////////////////////////////////
TICData::TICData() {
}



//////////////////////////////////////////////////////////////////////
TICData::~TICData() {
}



//////////////////////////////////////////////////////////////////////
void TICData::Clear() {
  fIC_Section.clear();
  fIC_Charge.clear();
  fIC_TS.clear();
}



//////////////////////////////////////////////////////////////////////
void TICData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TICData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  size_t mysize = fIC_Section.size();
  cout << "IC_Mult: " << mysize << endl;
  for(unsigned int i=0; i<mysize; i++){
    cout << "Charge section " << i+1 << " / Q= " << fIC_Charge[i] << endl;
  }
}
