#ifndef __MUSETTDATA__
#define __MUSETTDATA__
/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : June 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold MUSETT Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>
#include <map>
#include <iostream>

//ROOT
#include "TObject.h"

class TMUSETTData : public TObject {
   private:
      // First Layer
      // X strips
      // Energy
      std::vector<unsigned short>   fMUMU_DSSDXE_DetectorNbr; // Ce qui apparait
      std::vector<unsigned short>   fMUMU_DSSDXE_StripNbr; // Dans l'arbre root
      std::vector<double>           fMUMU_DSSDXE_Energy;
      // Time
      std::vector<unsigned short>   fMUMU_DSSDXT_DetectorNbr;
      std::vector<unsigned short>   fMUMU_DSSDXT_StripNbr;
      std::vector<double>           fMUMU_DSSDXT_Time;

      std::vector<double>           fMUMU_DSSDX_TimeStamp;
      std::vector<bool>             fMUMU_DSSDX_IsInterstrip;
      // Y strips
      // Energy
      std::vector<unsigned short>   fMUMU_DSSDYE_DetectorNbr;
      std::vector<unsigned short>   fMUMU_DSSDYE_StripNbr;
      std::vector<double>           fMUMU_DSSDYE_Energy;

      // Time
      std::vector<unsigned short>   fMUMU_DSSDYT_DetectorNbr;
      std::vector<unsigned short>   fMUMU_DSSDYT_StripNbr;
      std::vector<double>           fMUMU_DSSDYT_Time;

      unsigned short                fMUMU_TRIGGER;

      std::vector<double>           fMUMU_DSSDY_TimeStamp;
      std::vector<bool>             fMUMU_DSSDY_IsInterstrip;


    private:
     std::map<unsigned int, unsigned int> fMUMU_MapX;//!
     std::map<unsigned int, unsigned int> fMUMU_MapY;//!


   public:
      TMUSETTData();
      virtual ~TMUSETTData();

      void   Clear();
      void   Clear(const Option_t*) {};
      void   Dump() const;

      /////////////////////           SETTERS           ////////////////////////
      // FirstLayer

      // (X,E)
      public:
      inline void SetTRIGGER(const unsigned int Trig)
      {
        fMUMU_TRIGGER = Trig;
      }
      inline void   SetDSSDXE(const bool Map, const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Energy){
      if(Map)
        SetDSSDXE(DetNbr,fMUMU_MapX[StripNbr],Energy);
      else
        SetDSSDXE(DetNbr,StripNbr,Energy);
      }
      private:
      inline void   SetDSSDXE(const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Energy){
        fMUMU_DSSDXE_DetectorNbr.push_back(DetNbr);
        fMUMU_DSSDXE_StripNbr.push_back(StripNbr);
        fMUMU_DSSDXE_Energy.push_back(Energy);
      }

      // (X,T)
      public:
      inline void   SetDSSDXT(const bool Map, const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Time){
      if(Map)
        SetDSSDXT(DetNbr,fMUMU_MapX[StripNbr],Time);
      else
        SetDSSDXT(DetNbr,StripNbr,Time);

      }
      private:
      inline void   SetDSSDXT(const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Time){
        fMUMU_DSSDXT_DetectorNbr.push_back(DetNbr);
        fMUMU_DSSDXT_StripNbr.push_back(StripNbr);
        fMUMU_DSSDXT_Time.push_back(Time);
      }


      //(X timestamp )
      public :
      inline void   SetDSSDX_Tstamp(const double &Tstamp){
          fMUMU_DSSDX_TimeStamp.push_back(Tstamp);
        }
      inline void   SetDSSDX_Interstrip(const bool &inter){
              fMUMU_DSSDX_IsInterstrip.push_back(inter);
            }

      // (Y,E)
      public:
      inline void   SetDSSDYE(const bool Map, const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Energy){
      if(Map)
        SetDSSDYE(DetNbr,fMUMU_MapY[StripNbr],Energy);
      else
        SetDSSDYE(DetNbr,StripNbr,Energy);

      }
      private:
      inline void   SetDSSDYE(const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Energy){
        fMUMU_DSSDYE_DetectorNbr.push_back(DetNbr);
        fMUMU_DSSDYE_StripNbr.push_back(StripNbr);
        fMUMU_DSSDYE_Energy.push_back(Energy);
      }

      // (Y,T)
      public:
      inline void   SetDSSDYT(const bool Map, const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Time){
      if(Map)
        SetDSSDYT(DetNbr,fMUMU_MapY[StripNbr],Time);
      else
        SetDSSDYT(DetNbr,StripNbr,Time);

      }
      private:
      inline void   SetDSSDYT(const unsigned short& DetNbr, const unsigned short& StripNbr, const double& Time){
        fMUMU_DSSDYT_DetectorNbr.push_back(DetNbr);
        fMUMU_DSSDYT_StripNbr.push_back(StripNbr);
        fMUMU_DSSDYT_Time.push_back(Time);
      }


      //(Y timestamp, interstrip [for simulation] )
      public :
      inline void   SetDSSDY_Tstamp(const double &Tstamp){
          fMUMU_DSSDY_TimeStamp.push_back(Tstamp);
        }
      inline void   SetDSSDY_Interstrip(const bool &inter){
            fMUMU_DSSDY_IsInterstrip.push_back(inter);
          }

public:
      /////////////////////           GETTERS           ////////////////////////
      inline unsigned short   GetDSSDTRIG()                      const {return fMUMU_TRIGGER;}
      // DSSD
      // (X,E)
      inline unsigned short   GetDSSDXEMult()                      const {return fMUMU_DSSDXE_DetectorNbr.size();}
      inline unsigned short   GetDSSDXEDetectorNbr(const int& i)   const {return fMUMU_DSSDXE_DetectorNbr[i];}
      inline unsigned short   GetDSSDXEStripNbr(const int& i)      const {return fMUMU_DSSDXE_StripNbr[i];}
      inline double           GetDSSDXEEnergy(const int& i)        const {return fMUMU_DSSDXE_Energy[i];}
      // (X,T)
      inline unsigned short   GetDSSDXTMult()                      const {return fMUMU_DSSDXT_DetectorNbr.size();}
      inline unsigned short   GetDSSDXTDetectorNbr(const int& i)   const {return fMUMU_DSSDXT_DetectorNbr[i];}
      inline unsigned short   GetDSSDXTStripNbr(const int& i)      const {return fMUMU_DSSDXT_StripNbr[i];}
      inline double           GetDSSDXTTime(const int& i)          const {return fMUMU_DSSDXT_Time[i];}
      // (X, Timestamp/Interstrip)
      inline double           GetDSSDX_TStamp(const int& i)          const {return fMUMU_DSSDX_TimeStamp[i];}
      inline double           GetDSSDX_Interstrip(const int& i)      const {return fMUMU_DSSDX_IsInterstrip[i];}
      // (Y,E)
      inline unsigned short   GetDSSDYEMult()                      const {return fMUMU_DSSDYE_DetectorNbr.size();}
      inline unsigned short   GetDSSDYEDetectorNbr(const int& i)   const {return fMUMU_DSSDYE_DetectorNbr[i];}
      inline unsigned short   GetDSSDYEStripNbr(const int& i)       const {return fMUMU_DSSDYE_StripNbr[i];}
      inline double           GetDSSDYEEnergy(const int& i)        const {return fMUMU_DSSDYE_Energy[i];}
      // (Y,T)
      inline unsigned short   GetDSSDYTMult()                      const {return fMUMU_DSSDYT_DetectorNbr.size();}
      inline unsigned short   GetDSSDYTDetectorNbr(const int& i)   const {return fMUMU_DSSDYT_DetectorNbr[i];}
      inline unsigned short   GetDSSDYTStripNbr(const int& i)       const {return fMUMU_DSSDYT_StripNbr[i];}
      inline double           GetDSSDYTTime(const int& i)          const {return fMUMU_DSSDYT_Time[i];}
      // (Y, Timestamp/Interstrip)
      inline double           GetDSSDY_TStamp(const int& i)          const {return fMUMU_DSSDY_TimeStamp[i];}
      inline double           GetDSSDY_Interstrip(const int& i)     const {return fMUMU_DSSDY_IsInterstrip[i];}



      ClassDef(TMUSETTData,1)  // MUSETTData structure
};

#endif
