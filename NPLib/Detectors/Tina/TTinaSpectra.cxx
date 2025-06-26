/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Serge Franchoo  contact address: franchoo@ipno.in2p3.fr  *
 *                                                                           *
 * Creation Date  : February 2018                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Tina Spectra                                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// class header 
#include "TTinaSpectra.h"

// STL
#include <iostream>  
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"

////////////////////////////////////////////////////////////////////////////////
TTinaSpectra::TTinaSpectra() 
   : fNumberOfDetectors(0) {
  SetName("Tina");
}

////////////////////////////////////////////////////////////////////////////////
TTinaSpectra::TTinaSpectra(unsigned int NumberOfDetectors) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TTinaSpectra : Initialising control spectra for " 
      << NumberOfDetectors << " detectors" << endl
      << "************************************************" << endl ;
  SetName("Tina");
  fNumberOfDetectors = NumberOfDetectors;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}

////////////////////////////////////////////////////////////////////////////////
TTinaSpectra::~TTinaSpectra() {
}

////////////////////////////////////////////////////////////////////////////////
void TTinaSpectra::InitRawSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { 
    // energy 
    name = "Tina"+NPL::itoa(i+1)+"_ENERGY_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Tina/RAW");
    // time 
    name = "Tina"+NPL::itoa(i+1)+"_TIME_RAW";
    AddHisto1D(name, name, 4096, 0, 16384, "Tina/RAW");
  }
}

////////////////////////////////////////////////////////////////////////////////
void TTinaSpectra::InitPreTreatedSpectra() {
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetectors; i++) { 
    // energy 
    name = "Tina"+NPL::itoa(i+1)+"_Energy_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Tina/CAL");
    // time
    name = "Tina"+NPL::itoa(i+1)+"_Time_CAL";
    AddHisto1D(name, name, 500, 0, 25, "Tina/CAL");  
  }  
}

////////////////////////////////////////////////////////////////////////////////
void TTinaSpectra::InitPhysicsSpectra() {
  static string name;
  // kinematic plot 
  name = "Tina_ENERGY_TIME";
  AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "Tina/PHY");
}

////////////////////////////////////////////////////////////////////////////////
void TTinaSpectra::FillRawSpectra(TTinaData* RawData) {
  static string name;
  static string family;

  // YY1
  /* unsigned int sizeE = RawData->GetYY1_E_Mult();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Tina"+NPL::itoa(RawData->GetYY1_E_DetNbr(i))+"_Energy_RAW";
    family = "Tina/RAW";
    // this routine originates from .../NPLib/include/NPVSpectra.h    
    FillSpectra(family,name,RawData->GetYY1_Energy(i));
  }
  unsigned int sizeT = RawData->GetYY1_T_Mult();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "Tina"+NPL::itoa(RawData->GetYY1_T_DetNbr(i))+"_Time_RAW";
    family = "Tina/RAW";
    FillSpectra(family,name,RawData->GetYY1_Time(i));
  } */
}

////////////////////////////////////////////////////////////////////////////////
void TTinaSpectra::FillPreTreatedSpectra(TTinaData* PreTreatedData) {
  static string name;
  static string family;
  
  // YY1
  /* unsigned int sizeE = PreTreatedData->GetYY1_E_Mult();
  for (unsigned int i = 0; i < sizeE; i++) {
    name = "Tina"+NPL::itoa(PreTreatedData->GetYY1_E_DetNbr(i))+"_Energy_CAL";
    family = "Tina/CAL";
    FillSpectra(family,name,PreTreatedData->GetYY1_Energy(i));
  }
  unsigned int sizeT = PreTreatedData->GetYY1_T_Mult();
  for (unsigned int i = 0; i < sizeT; i++) {
    name = "Tina"+NPL::itoa(PreTreatedData->GetYY1_T_DetNbr(i))+"_Time_CAL";
    family = "Tina/CAL";
    FillSpectra(family,name,PreTreatedData->GetYY1_Time(i));
  } */
}

////////////////////////////////////////////////////////////////////////////////
void TTinaSpectra::FillPhysicsSpectra(TTinaPhysics* Physics) {
  static string name;
  static string family;
  family= "Tina/PHY";

  // energy vs time
  unsigned int sizeE = Physics->Energy.size();
  for(unsigned int i = 0 ; i < sizeE ; i++){
    name = "Tina_Energy_Time";
    FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
  }
}
