#ifndef __ISSDATA__
#define __ISSDATA__
/*****************************************************************************
 * Copyright (C) 2009-2019    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
* Author: M. Labiche                     address: marc.labiche@stfc.ac.uk    *
 *                                                                           *
 * Creation Date  : Jul 2019                                                 *
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


#include <stdlib.h>
#include <vector>
#include <map>
using namespace std ;

// ROOT
//#include "TNamed.h"
#include "TObject.h"

class TIssData : public TObject {
//class TIssData : public TNamed {
   private:
      // DSSD

     vector<UShort_t>   fIss_StripFront_DetectorNbr;
     vector<UShort_t>   fIss_StripFront_StripNbr;
     vector<Double_t>   fIss_StripFront_Energy;
     vector<Double_t>   fIss_StripFront_TimeCFD;
     vector<Double_t>   fIss_StripFront_TimeLED;
     vector<Double_t>   fIss_StripFront_Time;
  
     vector<UShort_t>   fIss_StripBack_DetectorNbr;
     vector<UShort_t>   fIss_StripBack_StripNbr;
     vector<Double_t>   fIss_StripBack_Energy;
     vector<Double_t>   fIss_StripBack_TimeCFD;
     vector<Double_t>   fIss_StripBack_TimeLED;
     vector<Double_t>   fIss_StripBack_Time;



/*

      // X strips
      // Energy
      vector<UShort_t>   fMM_StripXE_DetectorNbr;
      vector<UShort_t>   fMM_StripXE_StripNbr;
      vector<Double_t>   fMM_StripXE_Energy;
      // Time
      vector<UShort_t>   fMM_StripXT_DetectorNbr;
      vector<UShort_t>   fMM_StripXT_StripNbr;
      vector<Double_t>   fMM_StripXT_Time;
      // Y strips
      // Energy
      vector<UShort_t>   fMM_StripYE_DetectorNbr;
      vector<UShort_t>   fMM_StripYE_StripNbr;
      vector<Double_t>   fMM_StripYE_Energy;
      // Time
      vector<UShort_t>   fMM_StripYT_DetectorNbr;
      vector<UShort_t>   fMM_StripYT_StripNbr;
      vector<Double_t>   fMM_StripYT_Time;
*/

   public:
      TIssData();
      virtual ~TIssData();

      void   Clear();
      void  Clear(const Option_t*) {};
      void   Dump() const;

      /////////////////////           SETTERS           ////////////////////////
      // DSSD

      inline void SetFront_DetectorNbr(const UShort_t& DetNbr)  {fIss_StripFront_DetectorNbr.push_back(DetNbr);}
      inline void SetFront_StripNbr(const UShort_t& StripNbr)   {fIss_StripFront_StripNbr.push_back(StripNbr);}
      inline void SetFront_Energy(const Double_t& Energy)       {fIss_StripFront_Energy.push_back(Energy);}
      inline void SetFront_TimeCFD(const Double_t& TimeCFD)     {fIss_StripFront_TimeCFD.push_back(TimeCFD);}
      inline void SetFront_TimeLED(const Double_t& TimeLED)     {fIss_StripFront_TimeLED.push_back(TimeLED);}
      inline void SetFront_Time(const Double_t& Time)           {fIss_StripFront_Time.push_back(Time);}

      inline void SetBack_DetectorNbr(const UShort_t& DetNbr)   {fIss_StripBack_DetectorNbr.push_back(DetNbr);}
      inline void SetBack_StripNbr(const UShort_t& StripNbr)    {fIss_StripBack_StripNbr.push_back(StripNbr);}
      inline void SetBack_Energy(const Double_t& Energy)        {fIss_StripBack_Energy.push_back(Energy);}
      inline void SetBack_TimeCFD(const Double_t& TimeCFD)      {fIss_StripBack_TimeCFD.push_back(TimeCFD);}
      inline void SetBack_TimeLED(const Double_t& TimeLED)      {fIss_StripBack_TimeLED.push_back(TimeLED);}
      inline void SetBack_Time(const Double_t& Time)            {fIss_StripBack_Time.push_back(Time);}

      inline void SetFront(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Energy,const Double_t& TimeCFD,const    Double_t& TimeLED,const Double_t& Time = 0)	{
		SetFront_DetectorNbr(DetNbr);
		SetFront_StripNbr(StripNbr);
		SetFront_Energy(Energy);
		SetFront_TimeCFD(TimeCFD);
		SetFront_TimeLED(TimeLED);
		SetFront_Time(Time);
	};
	inline void SetBack(const UShort_t &DetNbr,const UShort_t &StripNbr,const Double_t &Energy,const Double_t &TimeCFD,const Double_t &TimeLED,const Double_t &Time = 0)	{
		SetBack_DetectorNbr(DetNbr);
		SetBack_StripNbr(StripNbr);
		SetBack_Energy(Energy);
		SetBack_TimeCFD(TimeCFD);
		SetBack_TimeLED(TimeLED);
		SetBack_Time(Time);
	};

/*
      // (X,E)
      inline void   SetStripXE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
        fMM_StripXE_DetectorNbr.push_back(DetNbr);
        fMM_StripXE_StripNbr.push_back(StripNbr);
        fMM_StripXE_Energy.push_back(Energy);
      }
      
      // (X,T)
     inline void   SetStripXT(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Time){   
      fMM_StripXT_DetectorNbr.push_back(DetNbr);  
      fMM_StripXT_StripNbr.push_back(StripNbr);       
      fMM_StripXT_Time.push_back(Time);  
     } 
       // (Y,E)
      inline void   SetStripYE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Energy){
        fMM_StripYE_DetectorNbr.push_back(DetNbr);
        fMM_StripYE_StripNbr.push_back(StripNbr);
        fMM_StripYE_Energy.push_back(Energy);
      }
      
      // (Y,T)
     inline void   SetStripYT(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& Time){   
      fMM_StripYT_DetectorNbr.push_back(DetNbr);  
      fMM_StripYT_StripNbr.push_back(StripNbr);       
      fMM_StripYT_Time.push_back(Time);  
     } 
      
 */
 
      /////////////////////           GETTERS           ////////////////////////
      // DSSD

      inline UShort_t GetFront_DetectorNbr(const unsigned int &i) const {return fIss_StripFront_DetectorNbr[i];}//!
      inline UShort_t GetFront_StripNbr(const unsigned int &i)    const {return fIss_StripFront_StripNbr[i];}//!
      inline Double_t GetFront_Energy(const unsigned int &i)      const {return fIss_StripFront_Energy[i];}//!
      inline Double_t GetFront_TimeCFD(const unsigned int &i)     const {return fIss_StripFront_TimeCFD[i];}//!
      inline Double_t GetFront_TimeLED(const unsigned int &i)     const {return fIss_StripFront_TimeLED[i];}//!
      inline Double_t GetFront_Time(const unsigned int &i)        const {return fIss_StripFront_Time[i];}//!


      inline UShort_t GetBack_DetectorNbr(const unsigned int &i) const {return fIss_StripBack_DetectorNbr[i];}//!
      inline UShort_t GetBack_StripNbr(const unsigned int &i)    const {return fIss_StripBack_StripNbr[i];}//!
      inline Double_t GetBack_Energy(const unsigned int &i)      const {return fIss_StripBack_Energy[i];}//!
      inline Double_t GetBack_TimeCFD(const unsigned int &i)     const {return fIss_StripBack_TimeCFD[i];}//!
      inline Double_t GetBack_TimeLED(const unsigned int &i)     const {return fIss_StripBack_TimeLED[i];}//!
      inline Double_t GetBack_Time(const unsigned int &i)        const {return fIss_StripBack_Time[i];}//!



      inline unsigned int GetMultiplicityFront() const {return fIss_StripFront_DetectorNbr.size();}//!
      inline unsigned int GetMultiplicityBack()  const {return fIss_StripBack_DetectorNbr.size();}//!




/*

      // (X,E)
      UShort_t   GetMMStripXEMult()                      const {return fMM_StripXE_DetectorNbr.size();}
      UShort_t   GetMMStripXEDetectorNbr(const Int_t& i) const {return fMM_StripXE_DetectorNbr[i];}
      UShort_t   GetMMStripXEStripNbr(const Int_t& i)    const {return fMM_StripXE_StripNbr[i];}
      Double_t   GetMMStripXEEnergy(const Int_t& i)      const {return fMM_StripXE_Energy[i];}
      // (X,T)
      UShort_t   GetMMStripXTMult()                      const {return fMM_StripXT_DetectorNbr.size();}
      UShort_t   GetMMStripXTDetectorNbr(const Int_t& i) const {return fMM_StripXT_DetectorNbr[i];}
      UShort_t   GetMMStripXTStripNbr(const Int_t& i)    const {return fMM_StripXT_StripNbr[i];}
      Double_t   GetMMStripXTTime(const Int_t& i)        const {return fMM_StripXT_Time[i];}
      // (Y,E)
      UShort_t   GetMMStripYEMult()                      const {return fMM_StripYE_DetectorNbr.size();}
      UShort_t   GetMMStripYEDetectorNbr(const Int_t& i) const {return fMM_StripYE_DetectorNbr[i];}
      UShort_t   GetMMStripYEStripNbr(const Int_t& i)    const {return fMM_StripYE_StripNbr[i];}
      Double_t   GetMMStripYEEnergy(const Int_t& i)      const {return fMM_StripYE_Energy[i];}
      // (Y,T)
      UShort_t   GetMMStripYTMult()                      const {return fMM_StripYT_DetectorNbr.size();}
      UShort_t   GetMMStripYTDetectorNbr(const Int_t& i) const {return fMM_StripYT_DetectorNbr[i];}
      UShort_t   GetMMStripYTStripNbr(const Int_t& i)    const {return fMM_StripYT_StripNbr[i];}
      Double_t   GetMMStripYTTime(const Int_t& i)        const {return fMM_StripYT_Time[i];}

*/

      ClassDef(TIssData,1)  // IssData structure
};

#endif
