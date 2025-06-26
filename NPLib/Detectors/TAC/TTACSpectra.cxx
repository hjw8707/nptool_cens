/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : July 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TAC Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TTACSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TTACSpectra::TTACSpectra() 
   : fNumberOfDetectors(0) {
  SetName("TAC");
}



////////////////////////////////////////////////////////////////////////////////
TTACSpectra::TTACSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TTACSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("TAC");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TTACSpectra::~TTACSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TTACSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "TAC"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "TAC/RAW");
    // Time 
    name = "TAC"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "TAC/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TTACSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "TAC"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "TAC/CAL");
    // Time
    name = "TAC"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "TAC/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TTACSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "TAC_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "TAC/PHY");
}



////////////////////////////////////////////////////////////////////////////////
/*void TTACSpectra::FillRawSpectra(TTACData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "TAC"+NPL::itoa(RawData->GetE_DetectorNbr(i))+"_ENERGY_RAW";
    family = "TAC/RAW";

    FillSpectra(family,name,RawData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = RawData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "TAC"+NPL::itoa(RawData->GetT_DetectorNbr(i))+"_TIME_RAW";
    family = "TAC/RAW";

    FillSpectra(family,name,RawData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TTACSpectra::FillPreTreatedSpectra(TTACData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "TAC"+NPL::itoa(PreTreatedData->GetE_DetectorNbr(i))+"_ENERGY_CAL";
    family = "TAC/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "TAC"+NPL::itoa(PreTreatedData->GetT_DetectorNbr(i))+"_TIME_CAL";
    family = "TAC/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Time(i));
  }
}
*/


////////////////////////////////////////////////////////////////////////////////
void TTACSpectra::FillPhysicsSpectra(TTACPhysics* Physics) {
  static string name;
  static string family;
  family= "TAC/PHY";

  // Energy vs time
  //unsigned int sizeE = Physics->.size();
  //for(unsigned int i = 0 ; i < sizeE ; i++){
  //  name = "TAC_ENERGY_TIME";
  //  FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  //}
}

