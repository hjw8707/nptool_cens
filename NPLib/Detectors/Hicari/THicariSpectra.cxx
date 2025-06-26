/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : march 2011                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for Hicari                 *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    + first version (not complete yet)                                     *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <string>
using namespace std;

// NPL
#include "THicariSpectra.h"
#include "NPOptionManager.h"
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif


////////////////////////////////////////////////////////////////////////////////
THicariSpectra::THicariSpectra(){
  SetName("Hicari");
  fnbinsRaw=16384;
  fbinRawMin=0;
  fbinRawMax=16384;
  fnbinsCal=5000;
  fbinCalMin=0;
  fbinCalMax=5000;

}

////////////////////////////////////////////////////////////////////////////////
THicariSpectra::THicariSpectra(unsigned int NumberOfClover){
/*
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "THicariSpectra : Initalising control spectra for " 
      << NumberOfClover << " Clover" << endl
      << "************************************************" << endl ;
  SetName("Hicari");
  fNumberOfClover = NumberOfClover;
  fNumberOfSegments=4;
  fNumberOfCores=4;
  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
  fnbinsRaw=4096;
  fbinRawMin=0;
  fbinRawMax=16384;
  fnbinsCal=5000;
  fbinCalMin=0;
  fbinCalMax=5000;
*/

}

////////////////////////////////////////////////////////////////////////////////
THicariSpectra::~THicariSpectra(){
}

////////////////////////////////////////////////////////////////////////////////
void THicariSpectra::InitRawSpectra(){
/*
  string name;
  for (unsigned int i = 0; i < fNumberOfClover; i++) { // loop on number of detectors
   for (unsigned int j = 0; j < fNumberOfCores; j++) { // loop on number of cores

    	name = Form("HicariEnergyRaw_Clover%d_ECC%d", i+1, j+1);
    	AddHisto1D(name, name, 16384, 0, 16384, "Hicari/ERAW/ECC");
    	name = Form("HicariTimeRaw_Clover%d_ECC%d", i+1, j+1);
    	AddHisto1D(name, name, 16384, 0, 16384, "Hicari/TRAW/ECC");

  	for (unsigned int k = 0; k < fNumberOfSegments; k++) { // loop on number of segments
    	name = Form("HicariEnergyRaw_Clover%d_ECC%d_GOCCE%d", i+1, j+1, k+1);
    	AddHisto1D(name, name, 16384, 0, 16384, "Hicari/ERAW/GOCCE");
    }
   }
  } // end loop on number of detectors
  */
}

////////////////////////////////////////////////////////////////////////////////
void THicariSpectra::InitPreTreatedSpectra()
{
/*
  string name;
  for (unsigned int i = 0; i < fNumberOfClover; i++) { // loop on number of detectors
   for (unsigned int j = 0; j < fNumberOfCores; j++) { // loop on number of cores
    	name = Form("HicariEnergyCal_Clover%d_ECC%d", i+1, j+1);
    	AddHisto1D(name, name, 5000, 0, 5000, "Hicari/ECal/ECC");
    	name = Form("HicariTimeCal_Clover%d_ECC%d", i+1, j+1);
    	AddHisto1D(name, name, 5000, 0, 5000, "Hicari/TCal/ECC");

  	for (unsigned int k = 0; k < fNumberOfSegments; k++) { // loop on number of segments
    	name = Form("HicariEnergyCal_Clover%d_ECC%d_GOCCE%d", i+1, j+1, k+1);
    	AddHisto1D(name, name, 5000, 0, 5000, "Hicari/ECal/GOCCE");

    }
   }
  } // end loop on number of detectors
*/
}

////////////////////////////////////////////////////////////////////////////////
void THicariSpectra::InitPhysicsSpectra(){
/*
  string name;
  name = "HicariEnergyAddBack";
  AddHisto1D(name, name, 5000, 0, 5000, "Hicari/DC");
*/
}



////////////////////////////////////////////////////////////////////////////////
void THicariSpectra::FillRawSpectra(THicariData* RawData){
/*
  string name;
  string family;

  // Energy 
  for (unsigned int i = 0; i < RawData->GetECCEMult(); i++){
    name = Form("HicariEnergyRaw_Clover%d_ECC%d", RawData->GetECCEClover(i)+1,RawData->GetECCECristal(i)+1);
    family = "Hicari/ERAW/ECC";
    FillSpectra(family,name
      ,RawData->GetECCEEnergy(i));
   }

  for (unsigned int i = 0; i < RawData->GetGOCCEEMult(); i++){
    name = Form("HicariEnergyRaw_Clover%d_ECC%d_GOCCE%d", RawData->GetGOCCEEClover(i)+1,RawData->GetGOCCEECristal(i)+1,RawData->GetGOCCEESegment(i)+1);
    family = "Hicari/ERAW/GOCCE";
    
  FillSpectra(family,name
      ,RawData->GetGOCCEEEnergy(i));
    }

  // Time
  for (unsigned int i = 0; i < RawData->GetECCTMult(); i++){
   name = Form("HicariTimeRaw_Clover%d_ECC%d", RawData->GetECCTClover(i)+1,RawData->GetECCTCristal(i)+1);
    family = "Hicari/RAW";

    FillSpectra(family,name
      ,RawData->GetECCTTime(i));
  }
*/
}

////////////////////////////////////////////////////////////////////////////////
void THicariSpectra::FillPreTreatedSpectra(THicariData* PreTreatedData){
/*
  string name ;
  string family;
  // Energy 
  for (unsigned int i = 0; i < PreTreatedData->GetECCEMult(); i++) {
    name = Form("HicariEnergyCal_Clover%d_ECC%d", PreTreatedData->GetECCEClover(i)+1,PreTreatedData->GetECCECristal(i)+1);
    family = "Hicari/ECal/ECC";

    FillSpectra(family,name
      ,PreTreatedData->GetECCEEnergy(i));
  }

  for (unsigned int i = 0; i < PreTreatedData->GetGOCCEEMult(); i++) {   
    name = Form("HicariEnergyCal_Clover%d_ECC%d_GOCCE%d", PreTreatedData->GetGOCCEEClover(i)+1,PreTreatedData->GetGOCCEECristal(i)+1,PreTreatedData->GetGOCCEESegment(i)+1);
    family = "Hicari/ECal/GOCCE";

    FillSpectra(family,name
      ,PreTreatedData->GetGOCCEEEnergy(i));
    }

  for (unsigned int i = 0; i < PreTreatedData->GetECCTMult(); i++) {
    name = Form("HicariTimeCal_Clover%d_ECC%d", PreTreatedData->GetECCTClover(i)+1,PreTreatedData->GetECCTCristal(i)+1);
    family = "Hicari/TCal/ECC";

    FillSpectra(family,name
      ,PreTreatedData->GetECCTTime(i));

  }
  */ 
}

////////////////////////////////////////////////////////////////////////////////
void THicariSpectra::FillPhysicsSpectra(THicariPhysics* Physics){
/*
  string name;
  string family= "Hicari/PHY";
  // Doppler Correct and Add Back
  name = "HicariEnergyAddBack";
  family = "Hicari/DC";
  for (unsigned int i = 0; i < Physics->DopplerCorrectedEnergy.size(); i++) {
    FillSpectra(family,name
      ,Physics->DopplerCorrectedEnergy[i]);
  }
*/
}
