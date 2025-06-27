/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : Octobre 2017                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Use to kill the beam track and replace it with the reaction product       *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#ifndef ThreeBody_h
#define ThreeBody_h

#include "G4VFastSimulationModel.hh"
#include "NPReaction.h"
#include "PhysicsList.hh"
#include "TReactionConditions.h"

#include "TF1.h"
#include "TFile.h"

class G4VPhysicalVolume;
namespace NPS {
class ThreeBody : public G4VFastSimulationModel {
public:
  ThreeBody(G4String, G4Region*);
  ThreeBody(G4String);
  ~ThreeBody();

public:
  void   ReadConfiguration();
  G4bool IsApplicable(const G4ParticleDefinition&);
  G4bool ModelTrigger(const G4FastTrack&);
  void   DoIt(const G4FastTrack&, G4FastStep&);

private:
  NPL::Reaction m_Reaction;
  string        m_BeamName;
  double        m_PreviousEnergy;
  double        m_PreviousLength;
  bool          m_active; // is the process active
  double        m_StepSize;

private:
  TReactionConditions* m_ReactionConditions;

public:
  void AttachReactionConditions();
  void SetStepSize(double step) { m_StepSize = step; };

public:
  bool         m_end;
  double       m_previousE;
  std::string  m_type;
  std::string  m_Beam;
  std::string  m_TargetNucleus;
  double       m_CompoundNucleusEnergy;
  double       m_CompoundNucleusStateWidth;
  std::string  m_IntermediaryNucleus;
  double       m_IntermediaryState;
  double       m_IntermediaryStateWidth;
  std::string  m_Light1;
  std::string  m_Heavy;
  std::string  m_Light2;
  double       m_ExcitationEnergyHeavy;
  int          m_ShootLight1;
  int          m_ShootLight2;
  int          m_ShootHeavy;
  NPL::Particle m_N1;
  NPL::Particle m_N2;

  NPL::Particle m_N3;
  NPL::Particle m_N4;

  NPL::Particle m_N5;
  NPL::Particle m_N6;

  bool        m_UserEnergyCS;
  std::string m_EnergyCSPath;
  std::string m_EnergyCSName;
  TFile*      m_FileEnergyCS;
  TF1*        m_EnergyCS;
  double      m_EnergyCSMax;
  TF1*        fIntermediaryStateEnergy;

public: // Managing the decay
        // Set everything for the decay
  void GetDecay(TVector3 p, double Etot, unsigned int nDecay, double* masses,
                double* T, double* theta, double* phi);
  TF1* GetIntermediaryStateEnergy(double Er, double Width);
};
} // namespace NPS

#endif
