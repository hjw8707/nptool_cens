/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 ****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 * Author: M. Labiche                     address: marc.labiche@stfc.ac.uk   *
 *                                                                           *
 * Creation Date  : July 2019                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for Iss                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// NPL
#include "TIssSpectra.h"
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
using namespace std;


////////////////////////////////////////////////////////////////////////////////
TIssSpectra::TIssSpectra(){
  SetName("ISS");
  fNumberOfDetector = 0;
  fStripFront=128;
  fStripBack=22;
}

////////////////////////////////////////////////////////////////////////////////
TIssSpectra::TIssSpectra(unsigned int NumberOfDetector){
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TIssSpectra : Initalising control spectra for " 
      << NumberOfDetector << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("ISS");
  fNumberOfDetector = NumberOfDetector;
  fStripFront=128;
  fStripBack=22;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}

////////////////////////////////////////////////////////////////////////////////
TIssSpectra::~TIssSpectra(){
}

////////////////////////////////////////////////////////////////////////////////
void TIssSpectra::InitRawSpectra(){

  static string name;
  for (unsigned int i = 0; i < fNumberOfDetector; i++) { // loop on number of detectors
    name = "ISSRaw"+NPL::itoa(i+1);
    // STR_FRONT_E_RAW
    name = "ISS"+NPL::itoa(i+1)+"_STR_FRONT_E_RAW";
    AddHisto2D(name, name, fStripFront, 1, fStripFront+1, 5000, 0, 1.5e6, "ISS/RAW/STR_FRONT_E")->Draw("colz");

    // STR_BACK_E_RAW
    name = "ISS"+NPL::itoa(i+1)+"_STR_BACK_E_RAW";
    AddHisto2D(name, name, fStripBack, 1, fStripBack+1, 5000, 0, 1.5e6, "ISS/RAW/STR_BACK_E")->Draw("colz");

    // STR_FRONT_EMAX_RAW
    name = "ISS"+NPL::itoa(i+1)+"_STR_FRONT_EMAX_RAW";
    AddHisto2D(name, name, fStripFront, 1, fStripFront+1, 5000, 0, 1.5e6, "ISS/RAW/STR_FRONT_EMAX");

    // STR_BACK_EMAX_Raw
    name = "ISS"+NPL::itoa(i+1)+"_STR_BACK_EMAX_RAW";
    AddHisto2D(name, name, fStripBack, 1, fStripBack+1, 5000, 0, 1.5e6, "ISS/RAW/STR_BACK_EMAX");

    // PAD_E_RAW
    name = "ISS"+NPL::itoa(i+1)+"_PAD_E_RAW";
    AddHisto1D(name, name, 500, 0, 2500, "ISS/RAW/PAD_E")->Draw("");

    // STR_FRONT_RAW_MULT
    name = "ISS"+NPL::itoa(i+1)+"_STR_FRONT_RAW_MULT";
    AddHisto1D(name, name, fStripFront, 1, fStripFront+1, "ISS/RAW/MULT")->Draw("");
    gPad->SetLogy();

    // STR_BACK_RAW_MULT
    name = "ISS"+NPL::itoa(i+1)+"_STR_BACK_RAW_MULT";
    AddHisto1D(name, name, fStripFront, 1, fStripFront+1, "ISS/RAW/MULT")->Draw("");
    gPad->SetLogy();

    // PAD_RAW_MULT
    name = "ISS"+NPL::itoa(i+1)+"_PAD_RAW_MULT";
    AddHisto1D(name, name, fNumberOfDetector, 1, fNumberOfDetector+1, "ISS/RAW/MULT")->Draw("");
    gPad->SetLogy();

  } // end loop on number of detectors

  // STR_PAD_DetN_MAP : useful for mapping issue
  name = "ISS_STR_PAD_DetN_RAW";
  AddHisto2D(name, name, fNumberOfDetector, 1, fNumberOfDetector+1, fNumberOfDetector, 1, fNumberOfDetector+1, "ISS/RAW/MAP")->Draw("colz");

}

////////////////////////////////////////////////////////////////////////////////
void TIssSpectra::InitPreTreatedSpectra(){
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetector; i++) { // loop on number of detectors
    // STR_FRONT_E_CAL
    name = "ISS"+NPL::itoa(i+1)+"_STR_FRONT_E_CAL";
    AddHisto2D(name, name, fStripFront, 1, fStripFront+1, 500, 0, 25, "ISS/CAL/STR_FRONT_E");

    // STR_BACK_E_CAL
    name = "ISS"+NPL::itoa(i+1)+"_STR_BACK_E_CAL";
    AddHisto2D(name, name, fStripBack, 1, fStripBack+1, 500, 0, 25, "ISS/CAL/STR_BACK_E");

    // PAD_E_CAL
    name = "ISS"+NPL::itoa(i+1)+"_PAD_E_CAL";
    AddHisto1D(name, name, 100, 0, 50, "ISS/CAL/PAD_E");

    // STR_FRONT_CAL_MULT
    name = "ISS"+NPL::itoa(i+1)+"_STR_FRONT_CAL_MULT";
    AddHisto1D(name, name, fStripFront, 1, fStripFront+1, "ISS/CAL/MULT");

    // STR_BACK_CAL_MULT
    name = "ISS"+NPL::itoa(i+1)+"_STR_BACK_CAL_MULT";
    AddHisto1D(name, name, fStripFront, 1, fStripFront+1, "ISS/CAL/MULT");

    // PAD_CAL_MULT
    name = "ISS"+NPL::itoa(i+1)+"_PAD_CAL_MULT";
    AddHisto1D(name, name, fNumberOfDetector, 1, fNumberOfDetector+1, "ISS/CAL/MULT");

    // PAD_CAL_ID 
    name = "ISS"+NPL::itoa(i+1)+"_PAD_CAL_ID";
    AddHisto2D(name, name,100,0,50,500,0,50, "ISS/CAL/ID");

    // Front-Back Energy Correlation
      name = "ISS"+NPL::itoa(i+1)+"_FB_COR";
      AddHisto2D(name, name,500,0,25,500,0,25, "ISS/CAL/FB"); 

  }  // end loop on number of detectors

  // STR_PAD_DetN_MAP : useful for mapping issue
  name = "ISS_STR_PAD_DetN_CAL";
  AddHisto2D(name, name, fNumberOfDetector, 1, fNumberOfDetector+1, fNumberOfDetector, 1, fNumberOfDetector+1, "ISS/CAL/MAP")->Draw("colz");

}

////////////////////////////////////////////////////////////////////////////////
void TIssSpectra::InitPhysicsSpectra(){
  static string name;
  // Kinematic Plot 
  name = "ISS_THETA_E";
  AddHisto2D(name, name,360,0,180,500,0,50,"ISS/PHY");

  // ID Plot
  // PAD-DE:
  name = "ISS_PAD_E_E";
  AddHisto1D(name, name,500,0,25,"ISS/PHY");

  for (unsigned int i = 0; i < fNumberOfDetector; i++) { // loop on number of detectors
    // PAD-DE:
    name = "ISS"+NPL::itoa(i+1)+"_PAD_E_E";
    AddHisto2D(name, name,100,0,100,500,0,25,"ISS/PHY");
  }
}



////////////////////////////////////////////////////////////////////////////////
void TIssSpectra::FillRawSpectra(TIssData* RawData){
  static string index;
  static string name;

  // STR_FRONT_E 
  unsigned int mysize = RawData->GetMultiplicityFront();
  double EFMAX = 0 ;
  int SFMAX = 0;
  int DFMAX = 0 ;
  for (unsigned int i = 0; i < mysize; i++) {
    index = "ISS/RAW/STR_FRONT_E/ISS"+NPL::itoa(RawData->GetFront_DetectorNbr(i))+"_STR_FRONT_E_RAW";
    if(RawData->GetFront_Energy(i) > EFMAX){
      EFMAX = RawData->GetFront_Energy(i);
      SFMAX = RawData->GetFront_StripNbr(i);
      DFMAX = RawData->GetFront_DetectorNbr(i);
    }
    
    FillSpectra(index
      ,RawData->GetFront_StripNbr(i), 
          RawData->GetFront_Energy(i));
  }
 
  if(DFMAX!=0){
    index = "ISS/RAW/STR_FRONT_EMAX/ISS"+NPL::itoa(DFMAX)+"_STR_FRONT_EMAX_RAW";
    FillSpectra(index,SFMAX, EFMAX);
  }
 
  // STR_BACK_E
  mysize = RawData->GetMultiplicityBack();
  double EBMAX = 0 ;
  int SBMAX = 0;
  int DBMAX = 0 ;
 
  for (unsigned int i = 0; i < mysize; i++) {
     index = "ISS/RAW/STR_BACK_E/ISS"+NPL::itoa( RawData->GetBack_DetectorNbr(i) )+"_STR_BACK_E_RAW";
    if(RawData->GetBack_Energy(i) > EBMAX){
      EBMAX = RawData->GetBack_Energy(i);
      SBMAX = RawData->GetBack_StripNbr(i);
      DBMAX = RawData->GetBack_DetectorNbr(i);
    }
   
    FillSpectra(index
      ,RawData->GetBack_StripNbr(i),
          RawData->GetBack_Energy(i));
  }
 
  if(DBMAX!=0){
    index = "ISS/RAW/STR_BACK_EMAX/ISS"+NPL::itoa(DBMAX)+"_STR_BACK_EMAX_RAW";
    FillSpectra(index,SBMAX, EBMAX);
  }


  // STR_FRONT MULT
  int myMULT[fNumberOfDetector];
  for( unsigned int i = 0; i < fNumberOfDetector; i++)
    myMULT[i] = 0 ; 

  for(unsigned int i = 0 ; i < RawData->GetMultiplicityFront();i++){
    myMULT[RawData->GetFront_DetectorNbr(i)-1] += 1;  
  }

  for( unsigned int i = 0; i < fNumberOfDetector; i++){
    index = "ISS/RAW/MULT/ISS"+NPL::itoa(i+1)+"_STR_FRONT_RAW_MULT";
    FillSpectra(index
      ,myMULT[i]);
  }

  // STR_BACK MULT
  for( unsigned int i = 0; i < fNumberOfDetector; i++)
    myMULT[i] = 0 ; 

  mysize = RawData->GetMultiplicityBack();
  for(unsigned int i = 0 ; i < mysize;i++){
    myMULT[RawData->GetBack_DetectorNbr(i)-1] += 1;  
  }

  for( unsigned int i = 0; i < fNumberOfDetector; i++){
    index= "ISS/RAW/MULT/ISS"+NPL::itoa(i+1)+"_STR_BACK_RAW_MULT";
    
    FillSpectra(index
      ,myMULT[i]);
  }



}

////////////////////////////////////////////////////////////////////////////////
void TIssSpectra::FillPreTreatedSpectra(TIssData* PreTreatedData){
  static string index;
  static string name;

  // Front-Back
  unsigned int mysizeF = PreTreatedData->GetMultiplicityFront();
  unsigned int mysizeB = PreTreatedData->GetMultiplicityBack();
//  unsigned int mysizePAD = PreTreatedData->GetMultiplicityPAD(); 

  for (unsigned int i = 0; i < mysizeF; i++) {
    for (unsigned int j = 0; j < mysizeB; j++) {
      if(PreTreatedData->GetFront_DetectorNbr(i)==PreTreatedData->GetBack_DetectorNbr(j)){
        index="ISS/CAL/FB/";
        name="ISS"+NPL::itoa(PreTreatedData->GetFront_DetectorNbr(i))+"_FB_COR";

      FillSpectra(index, name
        ,PreTreatedData->GetFront_Energy(i),
                PreTreatedData->GetBack_Energy(j) );
      }
    }
  } 

  // STR_FRONT_E
  unsigned int mysize = PreTreatedData->GetMultiplicityFront();
  for (unsigned int i = 0; i < mysize; i++) {
    index = "ISS/CAL/STR_FRONT_E";
    name="ISS"+NPL::itoa(PreTreatedData->GetFront_DetectorNbr(i))+"_STR_FRONT_E_CAL";

    FillSpectra(index,name
      ,PreTreatedData->GetFront_StripNbr(i), 
          PreTreatedData->GetFront_Energy(i));
  }
  // STR_BACK_E
  mysize = PreTreatedData->GetMultiplicityBack();
  for (unsigned int i = 0; i < mysize; i++) {
   index = "ISS/CAL/STR_BACK_E";
   string name = "ISS"+NPL::itoa( PreTreatedData->GetBack_DetectorNbr(i))+"_STR_BACK_E_CAL";

    FillSpectra(index,name
      ,PreTreatedData->GetBack_StripNbr(i), 
          PreTreatedData->GetBack_Energy(i));
  }


  // STR_FRONT MULT
  int myMULT[fNumberOfDetector];
  for( unsigned int i = 0; i < fNumberOfDetector; i++)
    myMULT[i] = 0 ; 

  mysize = PreTreatedData->GetMultiplicityFront(); 
  for(unsigned int i = 0 ; i < mysize ;i++){
    myMULT[PreTreatedData->GetFront_DetectorNbr(i)-1] += 1;  
  }

  for( unsigned int i = 0; i < fNumberOfDetector; i++){
    index= "ISS/CAL/MULT";
    name= "ISS"+NPL::itoa(i+1)+"_STR_FRONT_CAL_MULT";
    FillSpectra(index,name,myMULT[i]);
  }

  // STR_BACK MULT
  for( unsigned int i = 0; i < fNumberOfDetector; i++)
    myMULT[i] = 0 ; 

  mysize = PreTreatedData->GetMultiplicityBack();
  for(unsigned int i = 0 ; i < mysize ;i++){
    myMULT[PreTreatedData->GetBack_DetectorNbr(i)-1] += 1;  
  }

  for( unsigned int i = 0; i < fNumberOfDetector; i++){
    index= "ISS/CAL/MULT";
    name = "ISS"+NPL::itoa(i+1)+"_STR_BACK_CAL_MULT";
    FillSpectra(index,name
      ,myMULT[i]);
  }




}

////////////////////////////////////////////////////////////////////////////////
void TIssSpectra::FillPhysicsSpectra(TIssPhysics* Physics){
  static string index="ISS/PHY";
  static string name;
  // Kine plot
  unsigned int mysize = Physics->Strip_E.size();
  for(unsigned int i = 0 ; i < mysize ; i++){
    double Theta = Physics->GetPositionOfInteraction(i).Angle(TVector3(0,0,1));
    Theta = Theta/deg;
    double Etot=Physics->Strip_E[i];

    name = "ISS_THETA_E";
    FillSpectra(index,name,Theta,Etot);
  }
}

