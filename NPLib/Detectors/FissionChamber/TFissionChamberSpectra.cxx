/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : September 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold FissionChamber Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TFissionChamberSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TFissionChamberSpectra::TFissionChamberSpectra() 
   : fNumberOfDetectors(0) {
  SetName("FissionChamber");
}



////////////////////////////////////////////////////////////////////////////////
TFissionChamberSpectra::TFissionChamberSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TFissionChamberSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("FissionChamber");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TFissionChamberSpectra::~TFissionChamberSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TFissionChamberSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "FissionChamber"+NPL::itoa(i+1)+"_Q1_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "FissionChamber/RAW");
    // Time 
    name = "FissionChamber"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "FissionChamber/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TFissionChamberSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "FissionChamber"+NPL::itoa(i+1)+"_Q1_CAL";
    AddHisto1D(name, name, 500, 0, 25, "FissionChamber/CAL");
    // Time
    name = "FissionChamber"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "FissionChamber/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TFissionChamberSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "FissionChamber_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "FissionChamber/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TFissionChamberSpectra::FillRawSpectra(TFissionChamberData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMultiplicity();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "FissionChamber"+NPL::itoa(RawData->GetAnodeNbr(i))+"_Q1_RAW";
    family = "FissionChamber/RAW";
    FillSpectra(family,name,RawData->GetQ1(i));

    name = "FissionChamber"+NPL::itoa(RawData->GetAnodeNbr(i))+"_TIME_RAW";
    family = "FissionChamber/RAW";
    FillSpectra(family,name,RawData->GetTime(i));

  }
}



////////////////////////////////////////////////////////////////////////////////
void TFissionChamberSpectra::FillPreTreatedSpectra(TFissionChamberData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMultiplicity();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "FissionChamber"+NPL::itoa(PreTreatedData->GetAnodeNbr(i))+"_Q1_CAL";
    family = "FissionChamber/CAL";
    FillSpectra(family,name,PreTreatedData->GetQ1(i));
  
    name = "FissionChamber"+NPL::itoa(PreTreatedData->GetAnodeNbr(i))+"_TIME_CAL";
    family = "FissionChamber/CAL";
    FillSpectra(family,name,PreTreatedData->GetTime(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TFissionChamberSpectra::FillPhysicsSpectra(TFissionChamberPhysics* Physics) {
  static string name;
  static string family;
  family= "FissionChamber/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->Q1.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "FissionChamber_Q1_TIME";
    FillSpectra(family,name,Physics->Q1[i],Physics->Time[i]);
  }
}

