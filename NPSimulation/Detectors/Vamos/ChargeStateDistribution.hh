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

#ifndef ChargeStateDistribution_h
#define ChargeStateDistribution_h

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

  class ChargeStateDistribution : public G4VFastSimulationModel{
    public:
      ChargeStateDistribution (G4String, G4Region*);
      ChargeStateDistribution (G4String);
      ~ChargeStateDistribution ();

    public:
      G4bool IsApplicable(const G4ParticleDefinition&);
      G4bool ModelTrigger(const G4FastTrack &);
      void DoIt(const G4FastTrack&, G4FastStep&);

      void RungeKuttaPropagation (const G4FastTrack& fastTrack, G4FastStep& fastStep);

    public:
      void SetStepSize(double step){m_StepSize=step;};

    private:
      double m_StepSize;
      bool m_shoot;
  };
}


#endif 
