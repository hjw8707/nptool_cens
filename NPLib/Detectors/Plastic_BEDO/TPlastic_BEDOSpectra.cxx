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
 *  This class hold Plastic_BEDO Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TPlastic_BEDOSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TPlastic_BEDOSpectra::TPlastic_BEDOSpectra() 
   : fNumberOfDetectors(0) {
  SetName("Plastic_BEDO");
}



////////////////////////////////////////////////////////////////////////////////
TPlastic_BEDOSpectra::TPlastic_BEDOSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TPlastic_BEDOSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("Plastic_BEDO");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TPlastic_BEDOSpectra::~TPlastic_BEDOSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Plastic_BEDO"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Plastic_BEDO/RAW");
    // Time 
    name = "Plastic_BEDO"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Plastic_BEDO/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Plastic_BEDO"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Plastic_BEDO/CAL");
    // Time
    name = "Plastic_BEDO"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Plastic_BEDO/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "Plastic_BEDO_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "Plastic_BEDO/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOSpectra::FillRawSpectra(TPlastic_BEDOData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Plastic_BEDO"+NPL::itoa(RawData->GetE_DetectorNbr(i))+"_ENERGY_RAW";
    family = "Plastic_BEDO/RAW";

    FillSpectra(family,name,RawData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = RawData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "Plastic_BEDO"+NPL::itoa(RawData->GetT_DetectorNbr(i))+"_TIME_RAW";
    family = "Plastic_BEDO/RAW";

    FillSpectra(family,name,RawData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOSpectra::FillPreTreatedSpectra(TPlastic_BEDOData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Plastic_BEDO"+NPL::itoa(PreTreatedData->GetE_DetectorNbr(i))+"_ENERGY_CAL";
    family = "Plastic_BEDO/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "Plastic_BEDO"+NPL::itoa(PreTreatedData->GetT_DetectorNbr(i))+"_TIME_CAL";
    family = "Plastic_BEDO/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOSpectra::FillPhysicsSpectra(TPlastic_BEDOPhysics* Physics) {
  static string name;
  static string family;
  family= "Plastic_BEDO/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->Energy.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "Plastic_BEDO_ENERGY_TIME";
    FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  }
}

