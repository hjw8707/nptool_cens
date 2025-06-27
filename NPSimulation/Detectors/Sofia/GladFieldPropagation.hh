/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace                                         *
 * contact address: pierre.morfouace@cea.fr                                  *
 *                                                                           *
 * Creation Date  : January 2023                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription: Based On SamuraiFieldPropagation                              *
 * Use to kill the beam track and replace it with the reaction product       *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#ifndef GladFieldPropagation_h
#define GladFieldPropagation_h

#include "G4VFastSimulationModel.hh"
#include "G4Abla.hh"
#include "G4AblaInterface.hh"
#include "G4Fragment.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Region.hh"

#include "GladFieldMap.h"

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

namespace NPS{

  class GladFieldPropagation : public G4VFastSimulationModel{
    public:
      GladFieldPropagation (G4String, G4Region*);
      GladFieldPropagation (G4String);
      ~GladFieldPropagation ();

    public:
      G4bool IsApplicable(const G4ParticleDefinition&);
      G4bool ModelTrigger(const G4FastTrack &);
      void DoIt(const G4FastTrack&, G4FastStep&);

    private:
      bool m_Initialized; 
      double m_StepSize;  
      double m_Rmax;    
      double m_Current;
      double m_GladTiltAngle;
      TVector3 m_GladEntrance;
      TVector3 m_MW3_POS;
      double m_MW3_CentralTheta;
      string m_FieldMap;  
      GladFieldMap* m_Map; 

      void RungeKuttaPropagation (const G4FastTrack& fastTrack, G4FastStep& fastStep);

    public:
      void AttachReactionConditions();

      void SetStepSize(double step){m_StepSize=step;};
      void SetFieldMap(string fieldmap){m_FieldMap=fieldmap;};
      void SetRmax(double r_max){m_Rmax=r_max;};
      void SetCurrent(double val){m_Current=val;}
      void SetGladTiltAngle(double val){m_GladTiltAngle=val;}
      void SetGladEntrance(double x, double y, double z){m_GladEntrance=TVector3(x,y,z);}
      void Set_MWPC3_Position(double x, double y, double z){m_MW3_POS=TVector3(x,y,z);}
      void Set_MWPC3_CentralTheta(double val){m_MW3_CentralTheta=val;}
  };
}


#endif 
