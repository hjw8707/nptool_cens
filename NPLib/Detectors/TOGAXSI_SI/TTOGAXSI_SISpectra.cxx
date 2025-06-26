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
 *  This class hold TOGAXSI_SI Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TTOGAXSI_SISpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TTOGAXSI_SISpectra::TTOGAXSI_SISpectra() 
   : fNumberOfDetectors(0) {
  SetName("TOGAXSI_SI");
}



////////////////////////////////////////////////////////////////////////////////
TTOGAXSI_SISpectra::TTOGAXSI_SISpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TTOGAXSI_SISpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("TOGAXSI_SI");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TTOGAXSI_SISpectra::~TTOGAXSI_SISpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SISpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "TOGAXSI_SI"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "TOGAXSI_SI/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SISpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "TOGAXSI_SI"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "TOGAXSI_SI/CAL"); 
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SISpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "TOGAXSI_SI_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "TOGAXSI_SI/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SISpectra::FillRawSpectra(TTOGAXSI_SIData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetInnerXMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "TOGAXSI_SI"+NPL::itoa(RawData->GetInnerX_E_DetectorNbr(i))+"_ENERGY_RAW";
    family = "TOGAXSI_SI/RAW";

    FillSpectra(family,name,RawData->GetInnerX_E_Energy(i));
  }

}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SISpectra::FillPreTreatedSpectra(TTOGAXSI_SIData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetInnerXMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "TOGAXSI_SI"+NPL::itoa(PreTreatedData->GetInnerX_E_DetectorNbr(i))+"_ENERGY_CAL";
    family = "TOGAXSI_SI/CAL";

    FillSpectra(family,name,PreTreatedData->GetInnerX_E_Energy(i));
  }

}



////////////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SISpectra::FillPhysicsSpectra(TTOGAXSI_SIPhysics* Physics) {
  static string name;
  static string family;
  family= "TOGAXSI_SI/PHY";

  // Energy vs time
  unsigned int sizeXE = Physics->InnerXE.size();
  for(unsigned int i = 0 ; i < sizeXE ; i++){
    //name = "TOGAXSI_SI_ENERGY_TIME";
    //FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  }
}

