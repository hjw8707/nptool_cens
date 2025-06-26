/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Léo Plagnol  contact address: leo.plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : January 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Coaxial_Germanium Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TCoaxial_GermaniumSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TCoaxial_GermaniumSpectra::TCoaxial_GermaniumSpectra() 
   : fNumberOfDetectors(0) {
  SetName("Coaxial_Germanium");
}



////////////////////////////////////////////////////////////////////////////////
TCoaxial_GermaniumSpectra::TCoaxial_GermaniumSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TCoaxial_GermaniumSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("Coaxial_Germanium");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TCoaxial_GermaniumSpectra::~TCoaxial_GermaniumSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Coaxial_Germanium"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Coaxial_Germanium/RAW");
    // Time 
    name = "Coaxial_Germanium"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Coaxial_Germanium/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Coaxial_Germanium"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Coaxial_Germanium/CAL");
    // Time
    name = "Coaxial_Germanium"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Coaxial_Germanium/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "Coaxial_Germanium_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "Coaxial_Germanium/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumSpectra::FillRawSpectra(TCoaxial_GermaniumData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Coaxial_Germanium"+NPL::itoa(RawData->GetE_DetectorNbr(i))+"_ENERGY_RAW";
    family = "Coaxial_Germanium/RAW";

    FillSpectra(family,name,RawData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = RawData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "Coaxial_Germanium"+NPL::itoa(RawData->GetT_DetectorNbr(i))+"_TIME_RAW";
    family = "Coaxial_Germanium/RAW";

    FillSpectra(family,name,RawData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumSpectra::FillPreTreatedSpectra(TCoaxial_GermaniumData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Coaxial_Germanium"+NPL::itoa(PreTreatedData->GetE_DetectorNbr(i))+"_ENERGY_CAL";
    family = "Coaxial_Germanium/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "Coaxial_Germanium"+NPL::itoa(PreTreatedData->GetT_DetectorNbr(i))+"_TIME_CAL";
    family = "Coaxial_Germanium/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumSpectra::FillPhysicsSpectra(TCoaxial_GermaniumPhysics* Physics) {
  static string name;
  static string family;
  family= "Coaxial_Germanium/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->Energy.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "Coaxial_Germanium_ENERGY_TIME";
    FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  }
}

