/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold NeuLAND Raw data                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TNeuLANDData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TNeuLANDData)


//////////////////////////////////////////////////////////////////////
TNeuLANDData::TNeuLANDData() {
}



//////////////////////////////////////////////////////////////////////
TNeuLANDData::~TNeuLANDData() {
}



//////////////////////////////////////////////////////////////////////
void TNeuLANDData::Clear() {
  // 0 // 
  fNeuLAND_0_ID.clear();
  fNeuLAND_0_Sam.clear();
  fNeuLAND_0_Gtb.clear();
  fNeuLAND_0_Module.clear();
  fNeuLAND_0_Channel.clear();
  fNeuLAND_0_Cycle.clear();
  fNeuLAND_0_ADC.clear();
  fNeuLAND_0_TAC.clear();
  
  // 1 // 
  fNeuLAND_1_ID.clear();
  fNeuLAND_1_Sam.clear();
  fNeuLAND_1_Gtb.clear();
  fNeuLAND_1_Module.clear();
  fNeuLAND_1_Channel.clear();
  fNeuLAND_1_Cycle.clear();
  fNeuLAND_1_ADC.clear();
  fNeuLAND_1_TAC.clear();

  // Veto // 
  fNeuLAND_QdVETO_ID.clear();
  fNeuLAND_QdVETO_Charge.clear();
  fNeuLAND_TdVETO_ID.clear();
  fNeuLAND_TdVETO_Time.clear();
  fNeuLAND_QuVETO_ID.clear();
  fNeuLAND_QuVETO_Charge.clear();
  fNeuLAND_TuVETO_ID.clear();
  fNeuLAND_TuVETO_Time.clear();


}



//////////////////////////////////////////////////////////////////////
void TNeuLANDData::Dump() const {
   // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TNeuLANDData::Dump()] XXXXXXXXXXXXXXXXX" << endl;
  /*
  // Energy
  size_t mysize = fNeuLAND_E_DetectorNbr.size();
  cout << "NeuLAND_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fNeuLAND_E_DetectorNbr[i]
         << " Energy: " << fNeuLAND_Energy[i];
  }
  
  // Time
  mysize = fNeuLAND_T_DetectorNbr.size();
  cout << "NeuLAND_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fNeuLAND_T_DetectorNbr[i]
         << " Time: " << fNeuLAND_Time[i];
  }
  */
}
