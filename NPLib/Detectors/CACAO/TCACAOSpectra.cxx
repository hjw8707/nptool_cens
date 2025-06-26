/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: jwhwang@ibs.re.kr                        *
 *                                                                           *
 * Creation Date  : 4월 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold CACAO Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TCACAOSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TCACAOSpectra::TCACAOSpectra() 
   : fNumberOfDetectors(0) {
  SetName("CACAO");
}



////////////////////////////////////////////////////////////////////////////////
TCACAOSpectra::TCACAOSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TCACAOSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("CACAO");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TCACAOSpectra::~TCACAOSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TCACAOSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "CACAO"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "CACAO/RAW");
    // Time 
    name = "CACAO"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "CACAO/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TCACAOSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "CACAO"+NPL::itoa(i+1)+"_ENERGY_CAL";
    AddHisto1D(name, name, 500, 0, 25, "CACAO/CAL");
    // Time
    name = "CACAO"+NPL::itoa(i+1)+"_TIME_CAL";
    AddHisto1D(name, name, 500, 0, 25, "CACAO/CAL");

  
  }  // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TCACAOSpectra::InitPhysicsSpectra() {
  static string name;
  // Kinematic Plot 
  name = "CACAO_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "CACAO/PHY");
}



////////////////////////////////////////////////////////////////////////////////
void TCACAOSpectra::FillRawSpectra(TCACAOData* RawData) {
  static string name;
  static string family;

  // Energy 
  unsigned int sizeE = RawData->GetMult();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "CACAO"+NPL::itoa(RawData->GetDetN(i))+"_ENERGY_RAW";
    family = "CACAO/RAW";

    FillSpectra(family,name,RawData->GetE(i));
  }

  // Time
  unsigned int sizeT = RawData->GetMult();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "CACAO"+NPL::itoa(RawData->GetDetN(i))+"_TIME_RAW";
    family = "CACAO/RAW";

    FillSpectra(family,name,RawData->GetT(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TCACAOSpectra::FillPreTreatedSpectra(TCACAOData* PreTreatedData) {
  static string name;
  static string family;
  
  // Energy 
  unsigned int sizeE = PreTreatedData->GetMult();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "CACAO"+NPL::itoa(PreTreatedData->GetDetN(i))+"_ENERGY_CAL";
    family = "CACAO/CAL";

    FillSpectra(family,name,PreTreatedData->GetE(i));
  }

  // Time
  unsigned int sizeT = PreTreatedData->GetMult();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "CACAO"+NPL::itoa(PreTreatedData->GetDetN(i))+"_TIME_CAL";
    family = "CACAO/CAL";

    FillSpectra(family,name,PreTreatedData->GetT(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TCACAOSpectra::FillPhysicsSpectra(TCACAOPhysics* Physics) {
  static string name;
  static string family;
  family= "CACAO/PHY";

  // Energy vs time
  unsigned int sizeE = Physics->nhit;
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "CACAO_ENERGY_TIME";
    FillSpectra(family,name,Physics->E[i],Physics->T[i]);
  }
}

