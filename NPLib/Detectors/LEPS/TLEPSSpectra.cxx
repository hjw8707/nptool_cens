/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : March 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold LEPS Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TLEPSSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TLEPSSpectra::TLEPSSpectra() 
   : fNumberOfDetectors(0) {
  SetName("LEPS");
}



////////////////////////////////////////////////////////////////////////////////
TLEPSSpectra::TLEPSSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TLEPSSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("LEPS");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TLEPSSpectra::~TLEPSSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TLEPSSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "LEPS"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "LEPS/RAW");
    // Time 
    name = "LEPS"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "LEPS/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TLEPSSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "LEPS"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "LEPS/CAL");
    // Time
    name = "LEPS"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "LEPS/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TLEPSSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "LEPS_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "LEPS/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TLEPSSpectra::FillRawSpectra(TLEPSData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "LEPS"+NPL::itoa(RawData->GetE_DetectorNbr(i))+"_ENERGY_RAW";
    family = "LEPS/RAW";

    FillSpectra(family,name,RawData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = RawData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "LEPS"+NPL::itoa(RawData->GetT_DetectorNbr(i))+"_TIME_RAW";
    family = "LEPS/RAW";

    FillSpectra(family,name,RawData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TLEPSSpectra::FillPreTreatedSpectra(TLEPSData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "LEPS"+NPL::itoa(PreTreatedData->GetE_DetectorNbr(i))+"_ENERGY_CAL";
    family = "LEPS/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "LEPS"+NPL::itoa(PreTreatedData->GetT_DetectorNbr(i))+"_TIME_CAL";
    family = "LEPS/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TLEPSSpectra::FillPhysicsSpectra(TLEPSPhysics* Physics) {
  static string name;
  static string family;
  family= "LEPS/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->Energy.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "LEPS_ENERGY_TIME";
    FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  }
}

