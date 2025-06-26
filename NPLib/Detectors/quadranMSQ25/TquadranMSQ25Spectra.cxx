/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: leo.plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : January 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold quadranMSQ25 Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TquadranMSQ25Spectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TquadranMSQ25Spectra::TquadranMSQ25Spectra() 
   : fNumberOfDetectors(0) {
  SetName("quadranMSQ25");
}



////////////////////////////////////////////////////////////////////////////////
TquadranMSQ25Spectra::TquadranMSQ25Spectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TquadranMSQ25Spectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("quadranMSQ25");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TquadranMSQ25Spectra::~TquadranMSQ25Spectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Spectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "quadranMSQ25"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "quadranMSQ25/RAW");
    // Time 
    name = "quadranMSQ25"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "quadranMSQ25/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Spectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "quadranMSQ25"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "quadranMSQ25/CAL");
    // Time
    name = "quadranMSQ25"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "quadranMSQ25/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Spectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "quadranMSQ25_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "quadranMSQ25/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Spectra::FillRawSpectra(TquadranMSQ25Data* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "quadranMSQ25"+NPL::itoa(RawData->GetE_DetectorNbr(i))+"_ENERGY_RAW";
    family = "quadranMSQ25/RAW";

    FillSpectra(family,name,RawData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = RawData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "quadranMSQ25"+NPL::itoa(RawData->GetT_DetectorNbr(i))+"_TIME_RAW";
    family = "quadranMSQ25/RAW";

    FillSpectra(family,name,RawData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Spectra::FillPreTreatedSpectra(TquadranMSQ25Data* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "quadranMSQ25"+NPL::itoa(PreTreatedData->GetE_DetectorNbr(i))+"_ENERGY_CAL";
    family = "quadranMSQ25/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "quadranMSQ25"+NPL::itoa(PreTreatedData->GetT_DetectorNbr(i))+"_TIME_CAL";
    family = "quadranMSQ25/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Spectra::FillPhysicsSpectra(TquadranMSQ25Physics* Physics) {
  static string name;
  static string family;
  family= "quadranMSQ25/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->Energy.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "quadranMSQ25_ENERGY_TIME";
    FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  }
}

