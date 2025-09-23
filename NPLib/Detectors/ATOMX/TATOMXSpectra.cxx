/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: jungwoo  contact address: phyjics@gmail.com                        *
 *                                                                           *
 * Creation Date  : May 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ATOMX Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TATOMXSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TATOMXSpectra::TATOMXSpectra() 
   : fNumberOfDetectors(0) {
  SetName("ATOMX");
}



////////////////////////////////////////////////////////////////////////////////
TATOMXSpectra::TATOMXSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TATOMXSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("ATOMX");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TATOMXSpectra::~TATOMXSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TATOMXSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "ATOMX"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "ATOMX/RAW");
    // Time 
    name = "ATOMX"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "ATOMX/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TATOMXSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "ATOMX"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "ATOMX/CAL");
    // Time
    name = "ATOMX"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "ATOMX/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TATOMXSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "ATOMX_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "ATOMX/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TATOMXSpectra::FillRawSpectra(TATOMXData* RawData) {
  //static string name;
  //static string family;

  //// Energy 
  //unsigned int sizeE = RawData->GetMultEnergy();
  //for (unsigned int i = 0; i < sizeE; i++) {
  //  name = "ATOMX"+NPL::itoa(RawData->GetE_DetectorNbr(i))+"_ENERGY_RAW";
  //  family = "ATOMX/RAW";

  //  FillSpectra(family,name,RawData->Get_Energy(i));
  //}

  //// Time
  //unsigned int sizeT = RawData->GetMultTime();
  //for (unsigned int i = 0; i < sizeT; i++) {
  //  name = "ATOMX"+NPL::itoa(RawData->GetT_DetectorNbr(i))+"_TIME_RAW";
  //  family = "ATOMX/RAW";

  //  FillSpectra(family,name,RawData->Get_Time(i));
  //}
}



////////////////////////////////////////////////////////////////////////////////
void TATOMXSpectra::FillPreTreatedSpectra(TATOMXData* PreTreatedData) {
  //static string name;
  //static string family;
  //
  //// Energy 
  //unsigned int sizeE = PreTreatedData->GetMultEnergy();
  //for (unsigned int i = 0; i < sizeE; i++) {
  //  name = "ATOMX"+NPL::itoa(PreTreatedData->GetE_DetectorNbr(i))+"_ENERGY_CAL";
  //  family = "ATOMX/CAL";

  //  FillSpectra(family,name,PreTreatedData->Get_Energy(i));
  //}

  //// Time
  //unsigned int sizeT = PreTreatedData->GetMultTime();
  //for (unsigned int i = 0; i < sizeT; i++) {
  //  name = "ATOMX"+NPL::itoa(PreTreatedData->GetT_DetectorNbr(i))+"_TIME_CAL";
  //  family = "ATOMX/CAL";

  //  FillSpectra(family,name,PreTreatedData->Get_Time(i));
  //}
}



////////////////////////////////////////////////////////////////////////////////
void TATOMXSpectra::FillPhysicsSpectra(TATOMXPhysics* Physics) {
  //static string name;
  //static string family;
  //family= "ATOMX/PHY";

  //// Energy vs time
  //unsigned int sizeE = Physics->Energy.size();
  //for(unsigned int i = 0 ; i < sizeE ; i++){
  //  name = "ATOMX_ENERGY_TIME";
  //  FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  //}
}

