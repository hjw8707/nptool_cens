/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Louis Heitz  contact address: louis.heitz@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : mars 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold EdinburghDSSD Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// class header
#include "TEdinburghDSSDSpectra.h"
#include "NPOptionManager.h"
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

// STL
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <string>
#include <map>

using namespace std;

// NPTool header
#include "NPOptionManager.h"



////////////////////////////////////////////////////////////////////////////////
TEdinburghDSSDSpectra::TEdinburghDSSDSpectra(){
  SetName("EdinburghDSSD");
  fStripX=16;
  fStripY=16;
  fStripSecondLayer=16;
}



////////////////////////////////////////////////////////////////////////////////
TEdinburghDSSDSpectra::TEdinburghDSSDSpectra(std::map<int,int> DetectorIndex) {
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
  std::cout << "************************************************" << std::endl
    << "TEdinburghDSSDSpectra : Initalising control spectra for "
    << DetectorIndex.size() << " Detectors" << std::endl
    << "************************************************" << std::endl ;

  SetName("EdinburghDSSD");
  for(std::map<int,int>::iterator it=DetectorIndex.begin() ; it!=DetectorIndex.end() ; it++){
    fDetectorToIndex[it->second-1]=it->first;
    }

    fIndexToDetector=DetectorIndex;
    fNumberOfDetector=fDetectorToIndex.size();
    fStripX=16;
    fStripY=16;
    fStripSecondLayer=16;


  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}



////////////////////////////////////////////////////////////////////////////////
TEdinburghDSSDSpectra::~TEdinburghDSSDSpectra() {
}



////////////////////////////////////////////////////////////////////////////////
void TEdinburghDSSDSpectra::InitRawSpectra() {

  static string name;
  for (unsigned int i = 0; i < fDetectorToIndex.size(); i++) { // loop on number of detectors
    TString nbr = NPL::itoa(fDetectorToIndex[i]);

    // STRX_E_RAW
    name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_STRX_E_RAW";
    // AddHisto2D(name, name, fStripX, 1, fStripX+1, 512, 8192, 16384,  "EDIN/RAW/STRXE");
    AddHisto2D(name, name, fStripX, 1, fStripX+1, 10000, 7000, 17000,  "EDIN/RAW/STRXE");

    // STRY_E_RAW
    name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_STRY_E_RAW";
    // AddHisto2D(name, name, fStripY, 1, fStripY+1, 512, 0, 8192, "EDIN/RAW/STRYE");
    AddHisto2D(name, name, fStripY, 1, fStripY+1,10000, 0, 10000, "EDIN/RAW/STRYE");

    // STRX_T_RAW
    name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_STRX_T_RAW";
    AddHisto2D(name, name, fStripX, 1, fStripX+1, 512, 0, 8192, "EDIN/RAW/STRXT");

    // STRY_T_RAW
    name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_STRY_T_RAW";
    AddHisto2D(name, name, fStripY, 1, fStripY+1, 512, 0, 8192, "EDIN/RAW/STRYT");

    // SDLR_E_RAW
    name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_SDLR_E_RAW";
    AddHisto2D(name, name, fStripSecondLayer, 1, fStripSecondLayer+1, 512, 0, 8192, "EDIN/RAW/SDLRE");

    // SDLR_T_RAW
    name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_SDLR_T_RAW";
    AddHisto2D(name, name, fStripSecondLayer, 1, fStripSecondLayer+1, 512, 0, 8192, "EDIN/RAW/SDLRT");
  }
}



////////////////////////////////////////////////////////////////////////////////
void TEdinburghDSSDSpectra::InitPreTreatedSpectra() {

    string name;
    for (unsigned int i = 0; i < fNumberOfDetector; i++) { // loop on number of detectors
      // STRX_E_CAL
      name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_STRX_E_CAL";
      AddHisto2D(name, name, fStripX, 1, fStripX+1, 10000, 0, 50, "EDIN/CAL/STRXE");

      // STRY_E_CAL
      name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_STRY_E_CAL";
      AddHisto2D(name, name, fStripY, 1, fStripY+1, 10000, 0, 50, "EDIN/CAL/STRYE");

      // STR X-Y Correlation
      name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_STRXY_CORR_CAL";
      AddHisto2D(name, name, 500, 0, 50, 500, 0, 50, "EDIN/CAL/STRXY");

      // STRX_T_CAL
      name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_STRX_T_CAL";
      AddHisto2D(name, name, fStripX, 1, fStripX+1, 1000, 0, 1000, "EDIN/CAL/STRXT");

      // STRY_T_CAL
      name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_STRY_T_CAL";
      AddHisto2D(name, name, fStripY, 1, fStripY+1, 1000, 0, 1000, "EDIN/CAL/STRYT");

      // SDLR_E_CAL
      name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_SDLR_E_CAL";
      AddHisto2D(name, name, fStripSecondLayer, 1, fStripSecondLayer+1, 500, 0, 50, "EDIN/CAL/SDLRE");

      // SDLR_T_CAL
      name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_SDLR_T_CAL";
      AddHisto2D(name, name, fStripSecondLayer, 1, fStripSecondLayer+1, 500, 0, 50, "EDIN/CAL/SDLRT");

    }  // end loop on number of detectors

}



////////////////////////////////////////////////////////////////////////////////
void TEdinburghDSSDSpectra::InitPhysicsSpectra() {

  static string name;
  // X-Y Impact Matrix
  name = "MG_IMPACT_MATRIX";
  AddHisto2D(name, name,500,-150,150,500,-150,150, "EDIN/PHY");

  // X-Y Impact Matrix
  name = "MG_THETA_E";
  AddHisto2D(name, name,360,0,180,500,0,50,"EDIN/PHY");

  // ID Plot
  // E-TOF:
  name = "MG_E_TOF";
  AddHisto2D(name, name,500,0,50,500,200,1200,"EDIN/PHY");

  // SDLRE-DE:
  name = "MG_SDLRE_E";
  AddHisto2D(name, name,500,0,200,500,0,50,"EDIN/PHY");

  // Etot-DE:
  name = "MG_Etot_E";
  AddHisto2D(name, name,500,0,500,500,0,50,"EDIN/PHY");


  // ID plot detector by detector
  for (unsigned int i = 0; i < fNumberOfDetector; i++) { // loop on number of detectors
    // E-TOF:
    name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_E_TOF";
    AddHisto2D(name, name,500,0,50,500,200,1200,"EDIN/PHY");

    // SDLRE-DE:
    name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_SDLRE_E";
    AddHisto2D(name, name,500,0,200,500,0,50,"EDIN/PHY");

    // Etot-DE:
    name = "MG"+NPL::itoa(fDetectorToIndex[i])+"_Etot_E";
    AddHisto2D(name, name,500,0,500,500,0,50,"EDIN/PHY");
  }

}



////////////////////////////////////////////////////////////////////////////////
void TEdinburghDSSDSpectra::FillRawSpectra(TEdinburghDSSDData* RawData) {

    static string name;
    static string family;

    // STRX_E
    unsigned int size = RawData->GetDSSDXEMult();
    for (unsigned int i = 0; i < size ; i++) {
      name = "MG"+NPL::itoa( RawData->GetDSSDXEDetectorNbr(i))+"_STRX_E_RAW";
      family = "EDIN/RAW/STRXE";

      FillSpectra(family,name
          ,RawData->GetDSSDXEStripNbr(i),
          RawData->GetDSSDXEEnergy(i));
    }

    // STRY_E
    size = RawData->GetDSSDYEMult();
    for (unsigned int i = 0; i < size ; i++) {
      name = "MG"+NPL::itoa( RawData->GetDSSDYEDetectorNbr(i) )+"_STRY_E_RAW";
      family = "EDIN/RAW/STRYE";

      FillSpectra(family,name
          ,RawData->GetDSSDYEStripNbr(i),
          RawData->GetDSSDYEEnergy(i));
    }

    // STRX_T
    size = RawData->GetDSSDXTMult();
    for (unsigned int i = 0; i < size; i++) {
      name = "MG"+NPL::itoa(RawData->GetDSSDXTDetectorNbr(i))+"_STRX_T_RAW";
      family = "EDIN/RAW/STRXT";

      FillSpectra(family,name
          ,RawData->GetDSSDXTStripNbr(i),
          RawData->GetDSSDXTTime(i));
    }

    // STRY_T
    size = RawData->GetDSSDYTMult();
    for (unsigned int i = 0; i < size ; i++) {
      name = "MG"+NPL::itoa( RawData->GetDSSDYTDetectorNbr(i))+"_STRY_T_RAW";
      family = "EDIN/RAW/STRYT";

      FillSpectra(family,name
          ,RawData->GetDSSDYTStripNbr(i),
          RawData->GetDSSDYTTime(i));
    }

}



////////////////////////////////////////////////////////////////////////////////
void TEdinburghDSSDSpectra::FillPreTreatedSpectra(TEdinburghDSSDData* PreTreatedData) {

  static string name ;
  static string family;
  // STRX_E
  unsigned int sizeX = PreTreatedData->GetDSSDXEMult();
  for (unsigned int i = 0; i < sizeX ; i++) {
    name = "MG"+NPL::itoa( PreTreatedData->GetDSSDXEDetectorNbr(i))+"_STRX_E_CAL";
    family = "EDIN/CAL/STRXE";

    FillSpectra(family,name
        ,PreTreatedData->GetDSSDXEStripNbr(i),
        PreTreatedData->GetDSSDXEEnergy(i));
  }
  // STRY_E
  unsigned int sizeY = PreTreatedData->GetDSSDYEMult();
  for (unsigned int i = 0; i < sizeY; i++) {
    name = "MG"+NPL::itoa( PreTreatedData->GetDSSDYEDetectorNbr(i))+"_STRY_E_CAL";
    family = "EDIN/CAL/STRYE";

    FillSpectra(family,name
        ,PreTreatedData->GetDSSDYEStripNbr(i),
        PreTreatedData->GetDSSDYEEnergy(i));
  }

  // STR XY Correlation
  for (unsigned int i = 0; i < sizeX; i++) {
    for (unsigned int j = 0; j < sizeY; j++) {
      if(PreTreatedData->GetDSSDXEDetectorNbr(i)==PreTreatedData->GetDSSDYEDetectorNbr(j))
    name = "MG"+NPL::itoa( PreTreatedData->GetDSSDXEDetectorNbr(i) )+"_STRXY_CORR_CAL";
    family = "EDIN/CAL/STRXY";
    FillSpectra(family,name
        ,PreTreatedData->GetDSSDXEEnergy(i),
        PreTreatedData->GetDSSDYEEnergy(j));
    }
  }


  // STRX_T
  unsigned int size = PreTreatedData->GetDSSDXTMult();
  for (unsigned int i = 0; i < size; i++) {
    name = "MG"+NPL::itoa( PreTreatedData->GetDSSDXTDetectorNbr(i))+"_STRX_T_CAL";
    family = "EDIN/CAL/STRXT";

    FillSpectra(family,name
        ,PreTreatedData->GetDSSDXTStripNbr(i),
        PreTreatedData->GetDSSDXTTime(i));
  }
  // STRY_T
  size = PreTreatedData->GetDSSDYTMult();
  for (unsigned int i = 0; i < size; i++) {
    name = "MG"+NPL::itoa( PreTreatedData->GetDSSDYTDetectorNbr(i))+"_STRY_T_CAL";
    family = "EDIN/CAL/STRYT";

    FillSpectra(family,name
        ,PreTreatedData->GetDSSDYTStripNbr(i),
        PreTreatedData->GetDSSDYTTime(i));
  }
}



////////////////////////////////////////////////////////////////////////////////
void TEdinburghDSSDSpectra::FillPhysicsSpectra(TEdinburghDSSDPhysics* Physics) {

  static string name;
  static string family= "EDIN/PHY";
  // X-Y Impact Matrix

  for(unsigned int i = 0 ; i < Physics->DSSD_E.size(); i++){
    name = "MG_IMPACT_MATRIX";
    double x = Physics->GetPositionOfInteraction(i).x();
    double y = Physics->GetPositionOfInteraction(i).y();
    FillSpectra(family,name,x,y);

    name = "MG_THETA_E";
    double Theta = Physics->GetPositionOfInteraction(i).Angle(TVector3(0,0,1));
    Theta = Theta/deg;
    FillSpectra(family,name,Theta,Physics->DSSD_E[i]);


  }

}
