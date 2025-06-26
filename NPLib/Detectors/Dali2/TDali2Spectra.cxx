/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elidiano Tronchin  contact address: elidiano.tronchin@studenti.unipd.it                        *
 *                                                                           *
 * Creation Date  : septembre 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Dali Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TDali2Spectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TDali2Spectra::TDali2Spectra() 
   : fNumberOfDetectors(0) {
  SetName("Dali2");
}



////////////////////////////////////////////////////////////////////////////////
TDali2Spectra::TDali2Spectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TDali2Spectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("Dali2");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TDali2Spectra::~TDali2Spectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TDali2Spectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Dali2"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Dali2/RAW");
    // Time 
    name = "Dali2"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Dali2/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TDali2Spectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Dali2"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Dali2/CAL");
    // Time
    name = "Dali2"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Dali2/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TDali2Spectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "Dali2_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "Dali2/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TDali2Spectra::FillRawSpectra(TDali2Data* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Dali2"+NPL::itoa(RawData->GetE_DetectorNbr(i))+"_ENERGY_RAW";
    family = "Dali2/RAW";

    FillSpectra(family,name,RawData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = RawData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "Dali2"+NPL::itoa(RawData->GetT_DetectorNbr(i))+"_TIME_RAW";
    family = "Dali2/RAW";

    FillSpectra(family,name,RawData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TDali2Spectra::FillPreTreatedSpectra(TDali2Data* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Dali2"+NPL::itoa(PreTreatedData->GetE_DetectorNbr(i))+"_ENERGY_CAL";
    family = "Dali2/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "Dali2"+NPL::itoa(PreTreatedData->GetT_DetectorNbr(i))+"_TIME_CAL";
    family = "Dali2/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TDali2Spectra::FillPhysicsSpectra(TDali2Physics* Physics) {
  static string name;
  static string family;
  family= "Dali2/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->Energy.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "Dali2_ENERGY_TIME";
    FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  }
}

