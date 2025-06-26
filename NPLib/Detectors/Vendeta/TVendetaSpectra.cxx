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
 *  This class hold Vendeta Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TVendetaSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TVendetaSpectra::TVendetaSpectra() 
   : fNumberOfDetectors(0) {
  SetName("Vendeta");
}



////////////////////////////////////////////////////////////////////////////////
TVendetaSpectra::TVendetaSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TVendetaSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("Vendeta");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TVendetaSpectra::~TVendetaSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TVendetaSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Vendeta"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Vendeta/RAW");
    // Time 
    name = "Vendeta"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Vendeta/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TVendetaSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Vendeta"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Vendeta/CAL");
    // Time
    name = "Vendeta"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Vendeta/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TVendetaSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "Vendeta_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "Vendeta/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TVendetaSpectra::FillRawSpectra(TVendetaData* RawData) {
  static string name;
  static string family;

  unsigned int sizeE = RawData->GetLGMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Vendeta"+NPL::itoa(RawData->GetLGDetectorNbr(i))+"_ENERGY_RAW";
    family = "Vendeta/RAW";
    FillSpectra(family,name,RawData->GetLGQ1(i));

    name = "Vendeta"+NPL::itoa(RawData->GetLGDetectorNbr(i))+"_TIME_RAW";
    family = "Vendeta/RAW";
    FillSpectra(family,name,RawData->GetLGTime(i));
 
  }
}



////////////////////////////////////////////////////////////////////////////////
void TVendetaSpectra::FillPreTreatedSpectra(TVendetaData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetLGMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Vendeta"+NPL::itoa(PreTreatedData->GetLGDetectorNbr(i))+"_ENERGY_CAL";
    family = "Vendeta/CAL";
    FillSpectra(family,name,PreTreatedData->GetLGQ1(i));

    name = "Vendeta"+NPL::itoa(PreTreatedData->GetLGDetectorNbr(i))+"_TIME_CAL";
    family = "Vendeta/CAL";
    FillSpectra(family,name,PreTreatedData->GetLGTime(i));
 
  }
}



////////////////////////////////////////////////////////////////////////////////
void TVendetaSpectra::FillPhysicsSpectra(TVendetaPhysics* Physics) {
  static string name;
  static string family;
  family= "Vendeta/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->LG_DetectorNumber.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "Vendeta_Q1_Q2";
    FillSpectra(family,name,Physics->LG_Q1[i],Physics->LG_Q2[i]);

    name = "Vendeta_Q2_Time";
    FillSpectra(family,name,Physics->LG_Time[i],Physics->LG_Q2[i]);
  }
}

