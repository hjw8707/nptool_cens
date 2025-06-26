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
 *  This class hold SEASON Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
 #include "TSEASONData.h"

 #include <iostream>
 #include <fstream>
 #include <sstream>
 #include <string>
 using namespace std; 

 ClassImp(TSEASONData)


 //////////////////////////////////////////////////////////////////////
 TSEASONData::TSEASONData() {
 }



 //////////////////////////////////////////////////////////////////////
 TSEASONData::~TSEASONData() {
 }



 //////////////////////////////////////////////////////////////////////
 void TSEASONData::Clear() {
   // X Energy
   fSEASONX_E_DetectorNbr.clear();
   fSEASONX_E_StripNbr.clear();
   fSEASONX_Energy.clear();
   
   // X Time
   fSEASONX_T_DetectorNbr.clear();
   fSEASONX_T_StripNbr.clear();
   fSEASONX_Time.clear();
   
   // X Detected Particle ID (ID of primary vertex particle)
   fSEASONX_ParticleID.clear();

   // Y Energy
   fSEASONY_E_DetectorNbr.clear();
   fSEASONY_E_StripNbr.clear();
   fSEASONY_Energy.clear();
   
   // Y Time
   fSEASONY_T_DetectorNbr.clear();
   fSEASONY_T_StripNbr.clear();
   fSEASONY_Time.clear();
   
   // Y Detected Particle ID (ID of primary vertex particle)
   fSEASONY_ParticleID.clear();
   
   // Energy
   fSEASON_E_DetectorNbr.clear();
   fSEASON_E_StripNbrX.clear();
   fSEASON_E_StripNbrY.clear();
   fSEASON_Energy.clear();
   
   // Time
   fSEASON_T_DetectorNbr.clear();
   fSEASON_T_StripNbrX.clear();
   fSEASON_T_StripNbrY.clear();
   fSEASON_Time.clear();
 }



 //////////////////////////////////////////////////////////////////////
 void TSEASONData::Dump() const {
   // This method is very useful for debuging and worth the dev.
   cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TSEASONData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

   // X Energy
   size_t mysizeX = fSEASONX_E_DetectorNbr.size();
   cout << "SEASONX_E_Mult: " << mysizeX << endl;
  
   for (size_t i = 0 ; i < mysizeX ; i++){
     cout << "XDetNbr: " << fSEASONX_E_DetectorNbr[i]
 	 << " XStrip: " << fSEASONX_E_StripNbr[i]
          << " XEnergy: " << fSEASONX_Energy[i]
 	 << " XPartID: " << fSEASONX_ParticleID[i]
 	 << endl;
   }
   
   // X Time
   mysizeX = fSEASONX_T_DetectorNbr.size();
   cout << "SEASONX_T_Mult: " << mysizeX << endl;
  
   for (size_t i = 0 ; i < mysizeX ; i++){
     cout << "XDetNbr: " << fSEASONX_T_DetectorNbr[i]
 	 << " XStrip: " << fSEASONX_T_StripNbr[i]
          << " Time: " << fSEASONX_Time[i]
 	 << " XPartID: " << fSEASONX_ParticleID[i]
 	 << endl;
   }

   // Y Energy
   size_t mysizeY = fSEASONY_E_DetectorNbr.size();
   cout << "SEASONY_E_Mult: " << mysizeY << endl;
  
   for (size_t i = 0 ; i < mysizeY ; i++){
     cout << "YDetNbr: " << fSEASONY_E_DetectorNbr[i]
 	 << " YStrip: " << fSEASONY_E_StripNbr[i]
          << " YEnergy: " << fSEASONY_Energy[i]
 	 << " YPartID: " << fSEASONY_ParticleID[i]
 	 << endl;
   }
   
   // Y Time
   mysizeY = fSEASONY_T_DetectorNbr.size();
   cout << "SEASONY_T_Mult: " << mysizeY << endl;
  
   for (size_t i = 0 ; i < mysizeY ; i++){
     cout << "YDetNbr: " << fSEASONY_T_DetectorNbr[i]
 	 << " YStrip: " << fSEASONY_T_StripNbr[i]
          << " Time: " << fSEASONY_Time[i]
 	 << " YPartID: " << fSEASONY_ParticleID[i]
 	 << endl;
   }
 }
