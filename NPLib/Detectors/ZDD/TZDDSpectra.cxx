/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ZDD Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TZDDSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TZDDSpectra::TZDDSpectra() 
   : fNumberOfDetectors(0) {
  SetName("ZDD");
}



////////////////////////////////////////////////////////////////////////////////
TZDDSpectra::TZDDSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TZDDSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("ZDD");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TZDDSpectra::~TZDDSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TZDDSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "ZDD"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "ZDD/RAW");
    // Time 
    name = "ZDD"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "ZDD/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TZDDSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "ZDD"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "ZDD/CAL");
    // Time
    name = "ZDD"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "ZDD/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TZDDSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "ZDD_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "ZDD/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TZDDSpectra::FillRawSpectra(TZDDData* RawData) {
  static string name;
  static string family;

/*
  // Energy 
  unsigned int sizeE = RawData->GetMultEnergy("Plastic");
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "ZDD"+NPL::itoa(RawData->GetE_DetectorNbr("Plastic",i))+"_ENERGY_RAW";
    family = "ZDD/RAW";

    FillSpectra(family,name,RawData->Get_Energy("Plastic", i));
  }

  // Time
  unsigned int sizeT = RawData->GetMultTime("Plastic");
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "ZDD"+NPL::itoa(RawData->GetT_DetectorNbr("Plastic", i))+"_TIME_RAW";
    family = "ZDD/RAW";

    FillSpectra(family,name,RawData->Get_Time("Plastic", i));
  }
*/
}



////////////////////////////////////////////////////////////////////////////////
void TZDDSpectra::FillPreTreatedSpectra(TZDDData* PreTreatedData) {
  static string name;
  static string family;
  
  /*// Energy 
  unsigned int sizeE = PreTreatedData->GetMultEnergy("Plastic");
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "ZDD"+NPL::itoa(PreTreatedData->GetE_DetectorNbr("Plastic", i))+"_ENERGY_CAL";
    family = "ZDD/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Energy("Plastic", i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetMultTime("Plastic");
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "ZDD"+NPL::itoa(PreTreatedData->GetT_DetectorNbr("Plastic", i))+"_TIME_CAL";
    family = "ZDD/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Time("Plastic", i));
  }
*/
}



////////////////////////////////////////////////////////////////////////////////
void TZDDSpectra::FillPhysicsSpectra(TZDDPhysics* Physics) {
  static string name;
  static string family;
  family= "ZDD/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->PL_E.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "ZDD_ENERGY_TIME";
    FillSpectra(family,name,Physics->PL_E[i],Physics->PL_TS[i]);
  }
}

