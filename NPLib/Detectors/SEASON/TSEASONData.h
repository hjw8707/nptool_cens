#ifndef __SEASONDATA__
#define __SEASONDATA__
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
 
 // STL
 #include <vector>
 using namespace std;

 // ROOT
 #include "TObject.h"

 class TSEASONData : public TObject {
   //////////////////////////////////////////////////////////////
   // data members are hold into vectors in order 
   // to allow multiplicity treatment
 private:
   // DSSD X
   
   // Energy
   vector<UShort_t>   fSEASONX_E_DetectorNbr;
   vector<UShort_t>   fSEASONX_E_StripNbr;
   vector<Double_t>   fSEASONX_Energy;
   
   // Time
   vector<UShort_t>   fSEASONX_T_DetectorNbr;
   vector<UShort_t>   fSEASONX_T_StripNbr;
   vector<Double_t>   fSEASONX_Time;
   
   // Detected Particle ID (ID of primary vertex particle)
   vector<UInt_t>     fSEASONX_ParticleID;
   
   // DSSD Y
   
   // Energy
   vector<UShort_t>   fSEASONY_E_DetectorNbr;
   vector<UShort_t>   fSEASONY_E_StripNbr;
   vector<Double_t>   fSEASONY_Energy;
   
   // Time
   vector<UShort_t>   fSEASONY_T_DetectorNbr;
   vector<UShort_t>   fSEASONY_T_StripNbr;
   vector<Double_t>   fSEASONY_Time;
   
   // Detected Particle ID (ID of primary vertex particle)
   vector<UInt_t>     fSEASONY_ParticleID;
   
   // DSSD Final
   
   // Energy
   vector<UShort_t>   fSEASON_E_DetectorNbr;
   vector<UShort_t>   fSEASON_E_StripNbrX;
   vector<UShort_t>   fSEASON_E_StripNbrY;
   vector<Double_t>   fSEASON_Energy;
   
   // Time
   vector<UShort_t>   fSEASON_T_DetectorNbr;
   vector<UShort_t>   fSEASON_T_StripNbrX;
   vector<UShort_t>   fSEASON_T_StripNbrY;
   vector<Double_t>   fSEASON_Time;
   
   
   
   //////////////////////////////////////////////////////////////
   // Constructor and destructor
 public: 
   TSEASONData();
   ~TSEASONData();
   
   
   //////////////////////////////////////////////////////////////
   // Inherited from TObject and overriden to avoid warnings
 public:
   void Clear();
   void Clear(const Option_t*) {};
   void Dump() const;
   
   
   //////////////////////////////////////////////////////////////
   // Getters and Setters
   // Prefer inline declaration to avoid unnecessary called of 
   // frequently used methods
   // add //! to avoid ROOT creating dictionnary for the methods
 public:
   //////////////////////    SETTERS    ////////////////////////
   // X Energy
   inline void SetXEnergy(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Energy){
     fSEASONX_E_DetectorNbr.push_back(DetNbr);
     fSEASONX_E_StripNbr.push_back(StripNbr);
     fSEASONX_Energy.push_back(Energy);
   };//!
   
   // X Time
   inline void SetXTime(const UShort_t& DetNbr,const UShort_t& StripNbr, const Double_t& Time)	{
     fSEASONX_T_DetectorNbr.push_back(DetNbr);
     fSEASONX_T_StripNbr.push_back(StripNbr);
     fSEASONX_Time.push_back(Time);
   };//!
   
   // X Detected Particle ID
   inline void SetXParticleID(const UInt_t& PartID){
     fSEASONX_ParticleID.push_back(PartID);
   };//!
   
   // Add Energy and time to the last detected event
   inline void AddToLastX(const Double_t& Energy,const Double_t& Time){
     fSEASONX_Energy[fSEASONX_Energy.size()-1] += Energy;
     fSEASONX_Time[fSEASONX_Time.size()-1] += Time;
   };//!
   
   
   // Y Energy
   inline void SetYEnergy(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Energy){
     fSEASONY_E_DetectorNbr.push_back(DetNbr);
     fSEASONY_E_StripNbr.push_back(StripNbr);
     fSEASONY_Energy.push_back(Energy);
   };//!
   
   // Y Time
   inline void SetYTime(const UShort_t& DetNbr,const UShort_t& StripNbr, const Double_t& Time)	{
     fSEASONY_T_DetectorNbr.push_back(DetNbr);
     fSEASONY_T_StripNbr.push_back(StripNbr);
     fSEASONY_Time.push_back(Time);
   };//!
   
   // Y Detected Particle ID
   inline void SetYParticleID(const UInt_t& PartID){
     fSEASONY_ParticleID.push_back(PartID);
   };//!
   
   // Add Energy and time to the last detected event
   inline void AddToLastY(const Double_t& Energy,const Double_t& Time){
     fSEASONY_Energy[fSEASONY_Energy.size()-1] += Energy;
     fSEASONY_Time[fSEASONY_Time.size()-1] += Time;
   };//!
   
   // Energy
   inline void SetEnergy(const UShort_t& DetNbr,const UShort_t& StripNbrX,const UShort_t& StripNbrY,const Double_t& Energy){
     fSEASON_E_DetectorNbr.push_back(DetNbr);
     fSEASON_E_StripNbrX.push_back(StripNbrX);
     fSEASON_E_StripNbrY.push_back(StripNbrY);
     fSEASON_Energy.push_back(Energy);
   };//!
   
   // Time
   inline void SetTime(const UShort_t& DetNbr,const UShort_t& StripNbrX,const UShort_t& StripNbrY, const Double_t& Time)	{
     fSEASON_T_DetectorNbr.push_back(DetNbr);
     fSEASON_T_StripNbrX.push_back(StripNbrX);
     fSEASON_T_StripNbrY.push_back(StripNbrY);
     fSEASON_Time.push_back(Time);
   };//!
   
   
   //////////////////////    GETTERS    ////////////////////////
   // X Energy
   inline UShort_t GetXMultEnergy() const
   {return fSEASONX_E_DetectorNbr.size();}
   
   inline UShort_t GetXE_DetectorNbr(const unsigned int &i) const 
   {return fSEASONX_E_DetectorNbr[i];}//!
   
   inline vector<UShort_t> VGetXE_DetectorNbr() const 
   {return fSEASONX_E_DetectorNbr;}//!
   
   inline UShort_t GetXE_StripNbr(const unsigned int &i) const
   {return fSEASONX_E_StripNbr[i];}//!
   
   inline vector<UShort_t> VGetXE_StripNbr() const
   {return fSEASONX_E_StripNbr;}//!
   
   inline Double_t GetX_Energy(const unsigned int &i) const 
   {return fSEASONX_Energy[i];}//!
   
   inline vector<Double_t> VGetX_Energy() const 
   {return fSEASONX_Energy;}//!
   
   // X Time
   
   inline UShort_t GetXMultTime() const
   {return fSEASONX_T_DetectorNbr.size();}
   
   inline UShort_t GetXT_DetectorNbr(const unsigned int &i) const 
   {return fSEASONX_T_DetectorNbr[i];}//!
   
   inline vector<UShort_t> VGetXT_DetectorNbr() const 
   {return fSEASONX_T_DetectorNbr;}//!
   
   inline UShort_t GetXT_StripNbr(const unsigned int &i) const
   {return fSEASONX_T_StripNbr[i];}//!
   
   inline vector<UShort_t> VGetXT_StripNbr() const
   {return fSEASONX_T_StripNbr;}//!
   
   inline Double_t GetX_Time(const unsigned int &i) const 
   {return fSEASONX_Time[i];}//!
   
   inline vector<Double_t> VGetX_Time() const 
   {return fSEASONX_Time;}//!
   
   // X Particle ID
   inline UInt_t GetX_ParticleID(const unsigned int &i) const
   {return fSEASONX_ParticleID[i];}//!
   
   // Y Energy
   inline UShort_t GetYMultEnergy() const
   {return fSEASONY_E_DetectorNbr.size();}
   
   inline UShort_t GetYE_DetectorNbr(const unsigned int &i) const 
   {return fSEASONY_E_DetectorNbr[i];}//!
   
   inline vector<UShort_t> VGetYE_DetectorNbr() const 
   {return fSEASONY_E_DetectorNbr;}//!
   
   inline UShort_t GetYE_StripNbr(const unsigned int &i) const
   {return fSEASONY_E_StripNbr[i];}//!
   
   inline vector<UShort_t> VGetYE_StripNbr() const
   {return fSEASONY_E_StripNbr;}//!
   
   inline Double_t GetY_Energy(const unsigned int &i) const 
   {return fSEASONY_Energy[i];}//!
   
   inline vector<Double_t> VGetY_Energy() const 
   {return fSEASONY_Energy;}//!
   
   // Y Time
   
   inline UShort_t GetYMultTime() const
   {return fSEASONY_T_DetectorNbr.size();}
   
   inline UShort_t GetYT_DetectorNbr(const unsigned int &i) const 
   {return fSEASONY_T_DetectorNbr[i];}//!
   
   inline vector<UShort_t> VGetYT_DetectorNbr() const 
   {return fSEASONY_T_DetectorNbr;}//!
   
   inline UShort_t GetYT_StripNbr(const unsigned int &i) const
   {return fSEASONY_T_StripNbr[i];}//!
   
   inline vector<UShort_t> VGetYT_StripNbr() const
   {return fSEASONY_T_StripNbr;}//!
   
   inline Double_t GetY_Time(const unsigned int &i) const 
   {return fSEASONY_Time[i];}//!
   
   inline vector<Double_t> VGetY_Time() const 
   {return fSEASONY_Time;}//!
   
   // Y Particle ID
   inline UInt_t GetY_ParticleID(const unsigned int &i) const
   {return fSEASONY_ParticleID[i];}//!
   
   // Energy
   inline UShort_t GetMultEnergy() const
   {return fSEASON_E_DetectorNbr.size();}
   
   inline UShort_t GetE_DetectorNbr(const unsigned int &i) const 
   {return fSEASON_E_DetectorNbr[i];}//!
   
   inline vector<UShort_t> VGetE_DetectorNbr() const 
   {return fSEASON_E_DetectorNbr;}//!
   
   inline UShort_t GetE_StripNbrX(const unsigned int &i) const
   {return fSEASON_E_StripNbrX[i];}//!
   
   inline vector<UShort_t> VGetE_StripNbrX() const
   {return fSEASON_E_StripNbrX;}//!
   
   inline UShort_t GetE_StripNbrY(const unsigned int &i) const
   {return fSEASON_E_StripNbrY[i];}//
   
   inline vector<UShort_t> VGetE_StripNbrY() const
   {return fSEASON_E_StripNbrY;}//!
   
   inline Double_t Get_Energy(const unsigned int &i) const 
   {return fSEASON_Energy[i];}//!
   
   inline vector<Double_t> VGet_Energy() const 
   {return fSEASON_Energy;}//!
   
   // Time
   
   inline UShort_t GetMultTime() const
   {return fSEASON_T_DetectorNbr.size();}
   
   inline UShort_t GetT_DetectorNbr(const unsigned int &i) const 
   {return fSEASON_T_DetectorNbr[i];}//!
   
   inline vector<UShort_t> VGetT_DetectorNbr() const 
   {return fSEASON_T_DetectorNbr;}//!
   
   inline UShort_t GetT_StripNbrX(const unsigned int &i) const
   {return fSEASON_T_StripNbrX[i];}//!
   
   inline vector<UShort_t> VGetT_StripNbrX() const
   {return fSEASON_T_StripNbrX;}//!
   
   inline UShort_t GetT_StripNbrY(const unsigned int &i) const
   {return fSEASON_T_StripNbrY[i];}//
   
   inline vector<UShort_t> VGetT_StripNbrY() const
   {return fSEASON_T_StripNbrY;}//!
   
   inline Double_t Get_Time(const unsigned int &i) const 
   {return fSEASON_Time[i];}//!
   
   inline vector<Double_t> VGet_Time() const 
   {return fSEASON_Time;}//!

   
   //////////////////////////////////////////////////////////////
   // Required for ROOT dictionnary
   ClassDef(TSEASONData,1)  // SEASONData structure
 };

 #endif
