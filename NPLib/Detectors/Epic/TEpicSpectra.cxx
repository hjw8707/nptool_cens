/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey Chatillon  contact address: audrey.chatillon@cea.fr                        *
 *                                                                           *
 * Creation Date  : décembre 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Epic Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TEpicSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TEpicSpectra::TEpicSpectra() 
   : fNumberOfDetectors(0) {
  SetName("Epic");
}



////////////////////////////////////////////////////////////////////////////////
TEpicSpectra::TEpicSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TEpicSpectra : Initalizing control spectra for " 
      << NumberOfDetectors << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("Epic");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TEpicSpectra::~TEpicSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TEpicSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
    // Energy 
    name = "Epic"+NPL::itoa(i+1)+"_Q1_RAW";
    AddHisto1D(name, name, 5000, 0, 35000, "Epic/RAW");
    name = "Epic"+NPL::itoa(i+1)+"_Qmax_RAW";
    AddHisto1D(name, name, 5000, 0, 35000, "Epic/RAW");
  } // end loop on number of detectors
}



////////////////////////////////////////////////////////////////////////////////
void TEpicSpectra::InitPreTreatedSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TEpicSpectra::InitPhysicsSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TEpicSpectra::FillRawSpectra(TEpicData* RawData) {
  static string name;
  static string family;

  unsigned int sizeE = RawData->GetMultiplicity();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Epic"+NPL::itoa(RawData->GetAnodeNbr(i))+"_Q1_RAW";
    family = "Epic/RAW";
    FillSpectra(family,name,RawData->GetQ1(i));

    name = "Epic"+NPL::itoa(RawData->GetAnodeNbr(i))+"_Qmax_RAW";
    family = "Epic/RAW";
    FillSpectra(family,name,RawData->GetQmax(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TEpicSpectra::FillPreTreatedSpectra(TEpicData* PreTreatedData) {
}



////////////////////////////////////////////////////////////////////////////////
void TEpicSpectra::FillPhysicsSpectra(TEpicPhysics* Physics) {
}

