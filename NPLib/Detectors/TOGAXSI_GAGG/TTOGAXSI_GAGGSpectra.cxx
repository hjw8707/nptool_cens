/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Pohl  contact address: thomas.pohl@riken.jp                        *
 *                                                                           *
 * Creation Date  : February 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TOGAXSI_GAGG Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TTOGAXSI_GAGGSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TTOGAXSI_GAGGSpectra::TTOGAXSI_GAGGSpectra() 
//   : fNumberOfRecoilDetectors(0) {
   : fNumberOfDetectors(0) {
  SetName("TOGAXSI_GAGG");
}



////////////////////////////////////////////////////////////////////////////////
//TTOGAXSI_GAGGSpectra::TTOGAXSI_GAGGSpectra(unsigned int NumberOfRecoilDetectors, unsigned int NumberOfClusterDetectors) {
TTOGAXSI_GAGGSpectra::TTOGAXSI_GAGGSpectra(unsigned int NumberOfDetectors) {
    if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TTOGAXSI_GAGGSpectra : Initalizing control spectra for " 
//      << NumberOfRecoilDetectors << "Recoil Detectors" << endl
//      << NumberOfClusterDetectors << "Cluster Detectors" << endl
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("TOGAXSI_GAGG");
//  fNumberOfRecoilDetectors = NumberOfRecoilDetectors;
//  fNumberOfClusterDetectors = NumberOfClusterDetectors;
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
 
}



////////////////////////////////////////////////////////////////////////////////
TTOGAXSI_GAGGSpectra::~TTOGAXSI_GAGGSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGSpectra::InitRawSpectra() {
  static string name;
//  for (unsigned int i = 0; i < fNumberOfRecoilDetectors; i++) { // loop on number of detectors
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "TOGAXSI_GAGG"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "TOGAXSI_GAGG/RAW");
    // Time 
    name = "TOGAXSI_GAGG"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "TOGAXSI_GAGG/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGSpectra::InitPreTreatedSpectra() {
  static string name;
//  for (unsigned int i = 0; i < fNumberOfRecoilDetectors; i++) { // loop on number of detectors
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "TOGAXSI_GAGG"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "TOGAXSI_GAGG/CAL");
    // Time
    name = "TOGAXSI_GAGG"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "TOGAXSI_GAGG/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "TOGAXSI_GAGG_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "TOGAXSI_GAGG/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGSpectra::FillRawSpectra(TTOGAXSI_GAGGData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetRecoilMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "TOGAXSI_GAGG"+NPL::itoa(RawData->GetRecoil_E_DetectorNbr(i))+"_ENERGY_RAW";
    family = "TOGAXSI_GAGG/RAW";

    FillSpectra(family,name,RawData->GetRecoil_Energy(i));
  }

  // Time
  unsigned int sizeT = RawData->GetRecoilMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "TOGAXSI_GAGG"+NPL::itoa(RawData->GetRecoil_T_DetectorNbr(i))+"_TIME_RAW";
    family = "TOGAXSI_GAGG/RAW";

    FillSpectra(family,name,RawData->GetRecoil_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGSpectra::FillPreTreatedSpectra(TTOGAXSI_GAGGData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetRecoilMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "TOGAXSI_GAGG"+NPL::itoa(PreTreatedData->GetRecoil_E_DetectorNbr(i))+"_ENERGY_CAL";
    family = "TOGAXSI_GAGG/CAL";

    FillSpectra(family,name,PreTreatedData->GetRecoil_Energy(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetRecoilMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "TOGAXSI_GAGG"+NPL::itoa(PreTreatedData->GetRecoil_T_DetectorNbr(i))+"_TIME_CAL";
    family = "TOGAXSI_GAGG/CAL";

    FillSpectra(family,name,PreTreatedData->GetRecoil_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGSpectra::FillPhysicsSpectra(TTOGAXSI_GAGGPhysics* Physics) {
  static string name;
  static string family;
  family= "TOGAXSI_GAGG/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->RecoilE.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "TOGAXSI_GAGG_ENERGY_TIME";
    FillSpectra(family,name,Physics->RecoilE[i],Physics->RecoilT[i]);
  }
}

