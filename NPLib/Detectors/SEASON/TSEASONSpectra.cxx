/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Emmanuel Rey-herme                                       *
 * contact address: marine.vandebrouck@cea.fr                                *
 *                                                                           *
 * Creation Date  : septembre 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SEASON Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
 
 // class header 
 #include "TSEASONSpectra.h"

 // STL
 #include <iostream>  
 #include <string>
 using namespace std;

 // NPTool header
 #include "NPOptionManager.h"



 ////////////////////////////////////////////////////////////////////////////////
 TSEASONSpectra::TSEASONSpectra() {
   SetName("SEASON");
   fNumberOfDetectors = 0;
   fNumberOfStripsX = 0;
   fNumberOfStripsY = 0;
 }



 ////////////////////////////////////////////////////////////////////////////////
 TSEASONSpectra::TSEASONSpectra(unsigned int NumberOfDetectors, unsigned int NumberOfStripsX, unsigned int NumberOfStripsY) {
   if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
     cout << "************************************************" << endl
       << "TSEASONSpectra : Initalizing control spectra for " 
       << NumberOfDetectors << " Detectors" << endl
       << "************************************************" << endl ;
   SetName("SEASON");
   fNumberOfDetectors = NumberOfDetectors;
   fNumberOfStripsX = NumberOfStripsX;
   fNumberOfStripsY = NumberOfStripsY;

   InitRawSpectra();
   InitPreTreatedSpectra();
   InitPhysicsSpectra();
 }



 ////////////////////////////////////////////////////////////////////////////////
 TSEASONSpectra::~TSEASONSpectra() {
 }



 ////////////////////////////////////////////////////////////////////////////////
 void TSEASONSpectra::InitRawSpectra() {
   static string name;
   for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
     // X Energy 
     name = "SEASON"+NPL::itoa(i+1)+"_XENERGY_RAW";
     AddHisto2D(name, name,fNumberOfStripsX,1,fNumberOfStripsX+1, 4096, 0, 16384, "SEASON/RAWXE");
     // X Time 
     name = "SEASON"+NPL::itoa(i+1)+"_XTIME_RAW";
     AddHisto2D(name, name,fNumberOfStripsX,1,fNumberOfStripsY+1, 4096, 0, 16384, "SEASON/RAWXT");
     // Y Energy 
     name = "SEASON"+NPL::itoa(i+1)+"_YENERGY_RAW";
     AddHisto2D(name, name,fNumberOfStripsY,1,fNumberOfStripsY+1, 4096, 0, 16384, "SEASON/RAWYE");
     // Y Time 
     name = "SEASON"+NPL::itoa(i+1)+"_YTIME_RAW";
     AddHisto2D(name, name,fNumberOfStripsY,1,fNumberOfStripsY+1, 4096, 0, 16384, "SEASON/RAWYT");
   } // end loop on number of detectors
 }



 ////////////////////////////////////////////////////////////////////////////////
 void TSEASONSpectra::InitPreTreatedSpectra() {
   static string name;
   for (unsigned int i = 0; i < fNumberOfDetectors; i++) { // loop on number of detectors
     // X Energy 
     name = "SEASON"+NPL::itoa(i+1)+"_XENERGY_CAL";
     AddHisto2D(name, name,fNumberOfStripsX, 1, fNumberOfStripsX+1, 4096, 0, 25, "SEASON/CALXE");
     // X Time
     name = "SEASON"+NPL::itoa(i+1)+"_XTIME_CAL";
     AddHisto2D(name, name,fNumberOfStripsX,1,fNumberOfStripsX+1, 4096, 0, 50, "SEASON/CALXT");
     // Y Energy 
     name = "SEASON"+NPL::itoa(i+1)+"_YENERGY_CAL";
     AddHisto2D(name, name,fNumberOfStripsY, 1, fNumberOfStripsY+1, 4096, 0, 25, "SEASON/CALYE");
     // Y Time
     name = "SEASON"+NPL::itoa(i+1)+"_YTIME_CAL";
     AddHisto2D(name, name,fNumberOfStripsY,1,fNumberOfStripsY+1, 4096, 0, 50, "SEASON/CALYT");
   }  // end loop on number of detectors
 }



 ////////////////////////////////////////////////////////////////////////////////
 void TSEASONSpectra::InitPhysicsSpectra() {
   static string name;
   // Kinematic Plot 
   name = "SEASON_ENERGY_TIME";
   AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "SEASON/PHY");
 }



 ////////////////////////////////////////////////////////////////////////////////
 void TSEASONSpectra::FillRawSpectra(TSEASONData* RawData) {
   static string name;
   static string family;

   // X Energy 
   unsigned int sizeXE = RawData->GetXMultEnergy();
   for (unsigned int i = 0; i < sizeXE; i++) {
     name = "SEASON"+NPL::itoa(RawData->GetXE_DetectorNbr(i))+"_XENERGY_RAW";
     family = "SEASON/RAW";

     FillSpectra(family,name,RawData->GetXE_StripNbr(i),RawData->GetX_Energy(i));
   }

   // X Time
   unsigned int sizeXT = RawData->GetXMultTime();
   for (unsigned int i = 0; i < sizeXT; i++) {
     name = "SEASON"+NPL::itoa(RawData->GetXT_DetectorNbr(i))+"_XTIME_RAW";
     family = "SEASON/RAW";

     FillSpectra(family,name,RawData->GetXT_StripNbr(i),RawData->GetX_Time(i));
   }

   // Y Energy 
   unsigned int sizeYE = RawData->GetYMultEnergy();
   for (unsigned int i = 0; i < sizeYE; i++) {
     name = "SEASON"+NPL::itoa(RawData->GetYE_DetectorNbr(i))+"_YENERGY_RAW";
     family = "SEASON/RAW";

     FillSpectra(family,name,RawData->GetYE_StripNbr(i),RawData->GetY_Energy(i));
   }

   // Y Time
   unsigned int sizeYT = RawData->GetYMultTime();
   for (unsigned int i = 0; i < sizeYT; i++) {
     name = "SEASON"+NPL::itoa(RawData->GetYT_DetectorNbr(i))+"_YTIME_RAW";
     family = "SEASON/RAW";

     FillSpectra(family,name,RawData->GetYT_StripNbr(i),RawData->GetY_Time(i));
   }
 }



 ////////////////////////////////////////////////////////////////////////////////
 void TSEASONSpectra::FillPreTreatedSpectra(TSEASONData* PreTreatedData) {
   static string name;
   static string family;
   
   // X Energy 
   unsigned int sizeXE = PreTreatedData->GetXMultEnergy();
   for (unsigned int i = 0; i < sizeXE; i++) {
     name = "SEASON"+NPL::itoa(PreTreatedData->GetXE_DetectorNbr(i))+"_XENERGY_CAL";
     family = "SEASON/CAL";

     FillSpectra(family,name,PreTreatedData->GetXE_StripNbr(i),PreTreatedData->GetX_Energy(i));
   }

   // X Time
   unsigned int sizeXT = PreTreatedData->GetXMultTime();
   for (unsigned int i = 0; i < sizeXT; i++) {
     name = "SEASON"+NPL::itoa(PreTreatedData->GetXT_DetectorNbr(i))+"_XTIME_CAL";
     family = "SEASON/CAL";

     FillSpectra(family,name,PreTreatedData->GetXT_StripNbr(i),PreTreatedData->GetX_Time(i));
   }

   // Y Energy 
   unsigned int sizeYE = PreTreatedData->GetYMultEnergy();
   for (unsigned int i = 0; i < sizeYE; i++) {
     name = "SEASON"+NPL::itoa(PreTreatedData->GetYE_DetectorNbr(i))+"_YENERGY_CAL";
     family = "SEASON/CAL";

     FillSpectra(family,name,PreTreatedData->GetYE_StripNbr(i),PreTreatedData->GetY_Energy(i));
   }

   // Y Time
   unsigned int sizeYT = PreTreatedData->GetYMultTime();
   for (unsigned int i = 0; i < sizeYT; i++) {
     name = "SEASON"+NPL::itoa(PreTreatedData->GetYT_DetectorNbr(i))+"_YTIME_CAL";
     family = "SEASON/CAL";

     FillSpectra(family,name,PreTreatedData->GetYT_StripNbr(i),PreTreatedData->GetY_Time(i));
   }
 }



 ////////////////////////////////////////////////////////////////////////////////
 void TSEASONSpectra::FillPhysicsSpectra(TSEASONPhysics* Physics) {
   static string name;
   static string family;
   family= "SEASON/PHY";

   // Energy vs time
   unsigned int sizeE = Physics->Energy.size();
   for(unsigned int i = 0 ; i < sizeE ; i++){
     name = "SEASON_ENERGY_TIME";
     FillSpectra(family,name,Physics->Energy[i],Physics->Time[i]);
   }
 }

