#ifndef NPAnalysis_h 
#define NPAnalysis_h
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : march 2015                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Class describing the property of an Analysis object                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include"NPDetectorManager.h"
#include"TTreeReader.h"
namespace NPL{
  class VAnalysis{
    public:
     VAnalysis();
     ~VAnalysis();

    public: 
      virtual void TreatEvent(){};
      virtual bool UnallocateBeforeBuild(){return true;};
      virtual bool UnallocateBeforeTreat(){return true;};
      virtual bool FillOutputCondition(){return true;};
      virtual void Init(){};
      virtual void InitTreeReader(TTreeReader* TreeReader){};
      virtual void End(){};
      void SetDetectorManager(NPL::DetectorManager* det ) {m_DetectorManager=det;};
      //virtual TObject RawDataData(){}; 

    protected:
      NPL::DetectorManager* m_DetectorManager;
  };
}
#endif
