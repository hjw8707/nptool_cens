/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 * Author: M. Labiche                     address: marc.labiche@stfc.ac.uk   *
*                                                                            *
 * Creation Date  : July 2019                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Iss Raw data                                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

#include "TIssData.h"

ClassImp(TIssData)

TIssData::TIssData()
{
   // Default constructor

   // DSSD
  fIss_StripFront_DetectorNbr.clear();
  fIss_StripFront_StripNbr.clear();
  fIss_StripFront_Energy.clear();
  fIss_StripFront_TimeCFD.clear();
  fIss_StripFront_TimeLED.clear();
  fIss_StripFront_Time.clear();

  fIss_StripBack_DetectorNbr.clear();
  fIss_StripBack_StripNbr.clear();
  fIss_StripBack_Energy.clear();
  fIss_StripBack_TimeCFD.clear();
  fIss_StripBack_TimeLED.clear();
  fIss_StripBack_Time.clear();

/*


   // (X,E)
   fMM_StripXE_DetectorNbr.clear();
   fMM_StripXE_StripNbr.clear();
   fMM_StripXE_Energy.clear();
   // (X,T)
   fMM_StripXT_DetectorNbr.clear();
   fMM_StripXT_StripNbr.clear();
   fMM_StripXT_Time.clear();
   // (Y,E)
   fMM_StripYE_DetectorNbr.clear();   
   fMM_StripYE_StripNbr.clear();
   fMM_StripYE_Energy.clear();
   // (Y,T)
   fMM_StripYT_DetectorNbr.clear();
   fMM_StripYT_StripNbr.clear();
   fMM_StripYT_Time.clear();
*/
}

TIssData::~TIssData()
{}

void TIssData::Clear()
{
   // DSSD

  fIss_StripFront_DetectorNbr.clear();
  fIss_StripFront_StripNbr.clear();
  fIss_StripFront_Energy.clear();
  fIss_StripFront_TimeCFD.clear();
  fIss_StripFront_TimeLED.clear();
  fIss_StripFront_Time.clear();

  fIss_StripBack_DetectorNbr.clear();
  fIss_StripBack_StripNbr.clear();
  fIss_StripBack_Energy.clear();
  fIss_StripBack_TimeCFD.clear();
  fIss_StripBack_TimeLED.clear();
  fIss_StripBack_Time.clear();

/*
   // (X,E)
   fMM_StripXE_DetectorNbr.clear();
   fMM_StripXE_StripNbr.clear();
   fMM_StripXE_Energy.clear();
   // (X,T)
   fMM_StripXT_DetectorNbr.clear();
   fMM_StripXT_StripNbr.clear();
   fMM_StripXT_Time.clear();
   // (Y,E)
   fMM_StripYE_DetectorNbr.clear();
   fMM_StripYE_StripNbr.clear();
   fMM_StripYE_Energy.clear();
   // (Y,T)
   fMM_StripYT_DetectorNbr.clear();
   fMM_StripYT_StripNbr.clear();
   fMM_StripYT_Time.clear();
*/
}



void TIssData::Dump() const
{
   cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event XXXXXXXXXXXXXXXXX" << endl;

   // DSSD
  // Front
  cout << "Iss Strip Front Mult = " << fIss_StripFront_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fIss_StripFront_DetectorNbr.size(); i++){
    cout << "DetNbr (Front): " << fIss_StripFront_DetectorNbr[i]
         << "   Strip: " << fIss_StripFront_StripNbr[i]
         << "   Energy: " << fIss_StripFront_Energy[i]
         << "   Time CFD: " << fIss_StripFront_TimeCFD[i]
         << "   Time LED: " << fIss_StripFront_TimeLED[i] << endl;
  }

  // Back  
  cout << "Iss Strip Back Mult  = " << fIss_StripBack_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fIss_StripBack_DetectorNbr.size(); i++){
    cout << "DetNbr (Back): " << fIss_StripBack_DetectorNbr[i]
    << "   Strip: " << fIss_StripBack_StripNbr[i]
    << "   Energy: " << fIss_StripBack_Energy[i]
    << "   Time CFD: " << fIss_StripBack_TimeCFD[i]
    << "   Time LED: " << fIss_StripBack_TimeLED[i] << endl;
  }


/*
   // (X,E)
   cout << "MM_StripXE_Mult = " << fMM_StripXE_DetectorNbr.size() << endl;
   for (UShort_t i = 0; i < fMM_StripXE_DetectorNbr.size(); i++)
      cout << "DetNbr: " << fMM_StripXE_DetectorNbr[i] << " Strip: " << fMM_StripXE_StripNbr[i] << " Energy: " << fMM_StripXE_Energy[i] << endl;
   // (X,T)
   cout << "MM_StripXT_Mult = " << fMM_StripXT_DetectorNbr.size() << endl;
   for (UShort_t i = 0; i < fMM_StripXT_DetectorNbr.size(); i++)
      cout << "DetNbr: " << fMM_StripXT_DetectorNbr[i] << " Strip: " << fMM_StripXT_StripNbr[i] << " Time: " << fMM_StripXT_Time[i] << endl;
   // (Y,E)
   cout << "MM_StripYE_Mult = " << fMM_StripYE_DetectorNbr.size() << endl;
   for (UShort_t i = 0; i < fMM_StripYE_DetectorNbr.size(); i++)
      cout << "DetNbr: " << fMM_StripYE_DetectorNbr[i] << " Strip: " << fMM_StripYE_StripNbr[i] << " Energy: " << fMM_StripYE_Energy[i] << endl;
   // (Y,T)
   cout << "MM_StripYT_Mult = " << fMM_StripYT_DetectorNbr.size() << endl;
   for (UShort_t i = 0; i < fMM_StripYT_DetectorNbr.size(); i++)
      cout << "DetNbr: " << fMM_StripYT_DetectorNbr[i] << " Strip: " << fMM_StripYT_StripNbr[i] << " Time: " << fMM_StripYT_Time[i] << endl;
*/

}
