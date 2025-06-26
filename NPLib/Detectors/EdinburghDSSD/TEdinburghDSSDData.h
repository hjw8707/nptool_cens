#ifndef __EdinburghDSSDDATA__
#define __EdinburghDSSDDATA__
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
 *  This class hold EdinburghDSSD Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>
#include <map>
using namespace std;

// ROOT
#include "TObject.h"

class TEdinburghDSSDData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order
  // to allow multiplicity treatment
private:
   // First Layer
   // X strips
   // Energy
   std::vector<unsigned short>   fEdin_DSSDXE_DetectorNbr; // Ce qui apparait
   std::vector<unsigned short>   fEdin_DSSDXE_StripNbr; // Dans l'arbre root
   std::vector<double>           fEdin_DSSDXE_Energy;
   // Time
   std::vector<unsigned short>   fEdin_DSSDXT_DetectorNbr;
   std::vector<unsigned short>   fEdin_DSSDXT_StripNbr;
   std::vector<double>           fEdin_DSSDXT_Time;

   std::vector<double>           fEdin_DSSDX_TimeStamp;
   std::vector<bool>             fEdin_DSSDX_IsInterstrip;
   // Y strips
   // Energy
   std::vector<unsigned short>   fEdin_DSSDYE_DetectorNbr;
   std::vector<unsigned short>   fEdin_DSSDYE_StripNbr;
   std::vector<double>           fEdin_DSSDYE_Energy;

   // Time
   std::vector<unsigned short>   fEdin_DSSDYT_DetectorNbr;
   std::vector<unsigned short>   fEdin_DSSDYT_StripNbr;
   std::vector<double>           fEdin_DSSDYT_Time;

   std::vector<double>           fEdin_DSSDY_TimeStamp;
   std::vector<bool>             fEdin_DSSDY_IsInterstrip;


 private:
    std::map<unsigned int, unsigned int> fEDIN_MapX;//!   // Pour éviter d'écrire dans l'abre ROOT
    std::map<unsigned int, unsigned int> fEDIN_MapY;//!

/////////////////////////////////
  // Constructor and destructor
  public:
    TEdinburghDSSDData();
    ~TEdinburghDSSDData();


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
  // (X,E)
      public:
      inline void   SetDSSDXE(const bool Map, const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Energy){
      if(Map)
        SetDSSDXE(DetNbr,fEDIN_MapX[StripNbr],Energy);
      else
        SetDSSDXE(DetNbr,StripNbr,Energy);
      }
      private:
      inline void   SetDSSDXE(const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Energy){
        fEdin_DSSDXE_DetectorNbr.push_back(DetNbr);
        fEdin_DSSDXE_StripNbr.push_back(StripNbr);
        fEdin_DSSDXE_Energy.push_back(Energy);
      }

      // (X,T)
      public:
      inline void   SetDSSDXT(const bool Map, const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Time){
      if(Map)
        SetDSSDXT(DetNbr,fEDIN_MapX[StripNbr],Time);
      else
        SetDSSDXT(DetNbr,StripNbr,Time);

      }
      private:
      inline void   SetDSSDXT(const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Time){
        fEdin_DSSDXT_DetectorNbr.push_back(DetNbr);
        fEdin_DSSDXT_StripNbr.push_back(StripNbr);
        fEdin_DSSDXT_Time.push_back(Time);
      }


      //(X timestamp )
      public :
      inline void   SetDSSDX_Tstamp(const double &Tstamp){
          fEdin_DSSDX_TimeStamp.push_back(Tstamp);
        }
      inline void   SetDSSDX_Interstrip(const bool &inter){
              fEdin_DSSDX_IsInterstrip.push_back(inter);
            }

      // (Y,E)
      public:
      inline void   SetDSSDYE(const bool Map, const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Energy){
      if(Map)
        SetDSSDYE(DetNbr,fEDIN_MapY[StripNbr],Energy);
      else
        SetDSSDYE(DetNbr,StripNbr,Energy);

      }
      private:
      inline void   SetDSSDYE(const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Energy){
        fEdin_DSSDYE_DetectorNbr.push_back(DetNbr);
        fEdin_DSSDYE_StripNbr.push_back(StripNbr);
        fEdin_DSSDYE_Energy.push_back(Energy);
      }

      // (Y,T)
      public:
      inline void   SetDSSDYT(const bool Map, const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Time){
      if(Map)
        SetDSSDYT(DetNbr,fEDIN_MapY[StripNbr],Time);
      else
        SetDSSDYT(DetNbr,StripNbr,Time);

      }
      private:
      inline void   SetDSSDYT(const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Time){
        fEdin_DSSDYT_DetectorNbr.push_back(DetNbr);
        fEdin_DSSDYT_StripNbr.push_back(StripNbr);
        fEdin_DSSDYT_Time.push_back(Time);
      }


      //(Y timestamp )
      public :
      inline void   SetDSSDY_Tstamp(const double &Tstamp){
          fEdin_DSSDY_TimeStamp.push_back(Tstamp);
        }
      inline void   SetDSSDY_Interstrip(const bool &inter){
            fEdin_DSSDY_IsInterstrip.push_back(inter);
          }


      public:
            /////////////////////           GETTERS           ////////////////////////
            // DSSD
            // (X,E)
            inline unsigned short   GetDSSDXEMult()                      const {return fEdin_DSSDXE_DetectorNbr.size();}
            inline unsigned short   GetDSSDXEDetectorNbr(const int& i)   const {return fEdin_DSSDXE_DetectorNbr[i];}
            inline unsigned short   GetDSSDXEStripNbr(const int& i)      const {return fEdin_DSSDXE_StripNbr[i];}
            inline double           GetDSSDXEEnergy(const int& i)        const {return fEdin_DSSDXE_Energy[i];}
            // (X,T)
            inline unsigned short   GetDSSDXTMult()                      const {return fEdin_DSSDXT_DetectorNbr.size();}
            inline unsigned short   GetDSSDXTDetectorNbr(const int& i)   const {return fEdin_DSSDXT_DetectorNbr[i];}
            inline unsigned short   GetDSSDXTStripNbr(const int& i)      const {return fEdin_DSSDXT_StripNbr[i];}
            inline double           GetDSSDXTTime(const int& i)          const {return fEdin_DSSDXT_Time[i];}
            // (X, Timestamp)
            inline double           GetDSSDX_TStamp(const int& i)          const {return fEdin_DSSDX_TimeStamp[i];}
            inline double           GetDSSDX_Interstrip(const int& i)      const {return fEdin_DSSDX_IsInterstrip[i];}
            // (Y,E)
            inline unsigned short   GetDSSDYEMult()                      const {return fEdin_DSSDYE_DetectorNbr.size();}
            inline unsigned short   GetDSSDYEDetectorNbr(const int& i)   const {return fEdin_DSSDYE_DetectorNbr[i];}
            inline unsigned short   GetDSSDYEStripNbr(const int& i)       const {return fEdin_DSSDYE_StripNbr[i];}
            inline double           GetDSSDYEEnergy(const int& i)        const {return fEdin_DSSDYE_Energy[i];}
            // (Y,T)
            inline unsigned short   GetDSSDYTMult()                      const {return fEdin_DSSDYT_DetectorNbr.size();}
            inline unsigned short   GetDSSDYTDetectorNbr(const int& i)   const {return fEdin_DSSDYT_DetectorNbr[i];}
            inline unsigned short   GetDSSDYTStripNbr(const int& i)       const {return fEdin_DSSDYT_StripNbr[i];}
            inline double           GetDSSDYTTime(const int& i)          const {return fEdin_DSSDYT_Time[i];}
            // (Y, Timestamp)
            inline double           GetDSSDY_TStamp(const int& i)          const {return fEdin_DSSDY_TimeStamp[i];}
            inline double           GetDSSDY_Interstrip(const int& i)     const {return fEdin_DSSDY_IsInterstrip[i];}

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TEdinburghDSSDData,1)  // EdinburghDSSDData structure
};

#endif
