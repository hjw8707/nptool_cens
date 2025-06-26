/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : December 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold GeLi Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TGeLiSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TGeLiSpectra::TGeLiSpectra() 
   : fNumberOfDetectors(0) {
  SetName("GeLi");
}



////////////////////////////////////////////////////////////////////////////////
TGeLiSpectra::TGeLiSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TGeLiSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("GeLi");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TGeLiSpectra::~TGeLiSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TGeLiSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "GeLi"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "GeLi/RAW");
    // Time 
    name = "GeLi"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "GeLi/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TGeLiSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "GeLi"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "GeLi/CAL");
    // Time
    name = "GeLi"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "GeLi/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TGeLiSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "GeLi_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "GeLi/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TGeLiSpectra::FillRawSpectra(TGeLiData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "GeLi"+NPL::itoa(RawData->GetE_DetectorNbr(i))+"_ENERGY_RAW";
    family = "GeLi/RAW";

    FillSpectra(family,name,RawData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = RawData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "GeLi"+NPL::itoa(RawData->GetT_DetectorNbr(i))+"_TIME_RAW";
    family = "GeLi/RAW";

    FillSpectra(family,name,RawData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TGeLiSpectra::FillPreTreatedSpectra(TGeLiData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMultEnergy();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "GeLi"+NPL::itoa(PreTreatedData->GetE_DetectorNbr(i))+"_ENERGY_CAL";
    family = "GeLi/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Energy(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetMultTime();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "GeLi"+NPL::itoa(PreTreatedData->GetT_DetectorNbr(i))+"_TIME_CAL";
    family = "GeLi/CAL";

    FillSpectra(family,name,PreTreatedData->Get_Time(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TGeLiSpectra::FillPhysicsSpectra(TGeLiPhysics* Physics) {
  static string name;
  static string family;
  family= "GeLi/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->Energy.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "GeLi_ENERGY_TIME";
    FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  }
}

