/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold NeuLAND Spectra                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TNeuLANDSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TNeuLANDSpectra::TNeuLANDSpectra() 
   : fNumberOfDetectors(0) {
  SetName("NeuLAND");
}



////////////////////////////////////////////////////////////////////////////////
TNeuLANDSpectra::TNeuLANDSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TNeuLANDSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("NeuLAND");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TNeuLANDSpectra::~TNeuLANDSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TNeuLANDSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "NeuLAND"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "NeuLAND/RAW");
    // Time 
    name = "NeuLAND"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "NeuLAND/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TNeuLANDSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "NeuLAND"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "NeuLAND/CAL");
    // Time
    name = "NeuLAND"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "NeuLAND/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TNeuLANDSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "NeuLAND_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "NeuLAND/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TNeuLANDSpectra::FillRawSpectra(TNeuLANDData* RawData) {
/*  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "NeuLAND"+NPL::itoa(RawData->GetE_DetectorNbr(i))+"_ENERGY_RAW";
    family = "NeuLAND/RAW";

    FillSpectra(family,name,RawData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = RawData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "NeuLAND"+NPL::itoa(RawData->GetT_DetectorNbr(i))+"_TIME_RAW";
    family = "NeuLAND/RAW";

    FillSpectra(family,name,RawData->Get_Time(i));
  }*/
}



////////////////////////////////////////////////////////////////////////////////
void TNeuLANDSpectra::FillPreTreatedSpectra(TNeuLANDData* PreTreatedData) {
/*  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "NeuLAND"+NPL::itoa(PreTreatedData->GetE_DetectorNbr(i))+"_ENERGY_CAL";
    family = "NeuLAND/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "NeuLAND"+NPL::itoa(PreTreatedData->GetT_DetectorNbr(i))+"_TIME_CAL";
    family = "NeuLAND/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Time(i));
  }*/
}



////////////////////////////////////////////////////////////////////////////////
void TNeuLANDSpectra::FillPhysicsSpectra(TNeuLANDPhysics* Physics) {
/*  static string name;
  static string family;
  family= "NeuLAND/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->Energy.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "NeuLAND_ENERGY_TIME";
    FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  }*/
}

