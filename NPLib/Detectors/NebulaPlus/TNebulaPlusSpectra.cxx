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
 *  This class hold NebulaPlus Spectra                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TNebulaPlusSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TNebulaPlusSpectra::TNebulaPlusSpectra() 
   : fNumberOfDetectors(0) {
  SetName("NebulaPlus");
}



////////////////////////////////////////////////////////////////////////////////
TNebulaPlusSpectra::TNebulaPlusSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TNebulaPlusSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("NebulaPlus");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TNebulaPlusSpectra::~TNebulaPlusSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TNebulaPlusSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "NebulaPlus"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "NebulaPlus/RAW");
    // Time 
    name = "NebulaPlus"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "NebulaPlus/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TNebulaPlusSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "NebulaPlus"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "NebulaPlus/CAL");
    // Time
    name = "NebulaPlus"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "NebulaPlus/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TNebulaPlusSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "NebulaPlus_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "NebulaPlus/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TNebulaPlusSpectra::FillRawSpectra(TNebulaPlusData* RawData) {
/*  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "NebulaPlus"+NPL::itoa(RawData->GetE_DetectorNbr(i))+"_ENERGY_RAW";
    family = "NebulaPlus/RAW";

    FillSpectra(family,name,RawData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = RawData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "NebulaPlus"+NPL::itoa(RawData->GetT_DetectorNbr(i))+"_TIME_RAW";
    family = "NebulaPlus/RAW";

    FillSpectra(family,name,RawData->Get_Time(i));
  }*/
}



////////////////////////////////////////////////////////////////////////////////
void TNebulaPlusSpectra::FillPreTreatedSpectra(TNebulaPlusData* PreTreatedData) {
/*  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "NebulaPlus"+NPL::itoa(PreTreatedData->GetE_DetectorNbr(i))+"_ENERGY_CAL";
    family = "NebulaPlus/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "NebulaPlus"+NPL::itoa(PreTreatedData->GetT_DetectorNbr(i))+"_TIME_CAL";
    family = "NebulaPlus/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Time(i));
  }*/
}



////////////////////////////////////////////////////////////////////////////////
void TNebulaPlusSpectra::FillPhysicsSpectra(TNebulaPlusPhysics* Physics) {
/*  static string name;
  static string family;
  family= "NebulaPlus/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->Energy.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "NebulaPlus_ENERGY_TIME";
    FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  }*/
}

