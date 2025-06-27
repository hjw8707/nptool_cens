/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Valerian Alcindor                                        * 
 * contact address: valcindor@ikp.tu-darmstadt.de                            *
 *                                                                           *
 * Creation Date  : 2019                                                     * 
 * Last update    : August 2020                                              *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Generation of sequential three body decay from resonant states            *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include "ThreeBody.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "NPFunction.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "Particle.hh"
#include "RootOutput.h"
#include "TGenPhaseSpace.h"
#include <iostream>
#include <string>
////////////////////////////////////////////////////////////////////////////////
NPS::ThreeBody::ThreeBody(G4String modelName, G4Region* envelope)
    : G4VFastSimulationModel(modelName, envelope) {
  ReadConfiguration();
  m_PreviousEnergy = 0;
  m_PreviousLength = 0;
}

////////////////////////////////////////////////////////////////////////////////
NPS::ThreeBody::ThreeBody(G4String modelName)
    : G4VFastSimulationModel(modelName) {}

////////////////////////////////////////////////////////////////////////////////
NPS::ThreeBody::~ThreeBody() {}

////////////////////////////////////////////////////////////////////////////////
void NPS::ThreeBody::AttachReactionConditions() {
  // Reasssigned the branch address
  if (RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
    RootOutput::getInstance()->GetTree()->SetBranchAddress(
        "ReactionConditions", &m_ReactionConditions);
}

////////////////////////////////////////////////////////////////////////////////
void NPS::ThreeBody::ReadConfiguration() {
  NPL::InputParser input(NPOptionManager::getInstance()->GetReactionFile());

  vector<NPL::InputBlock*> blocks = input.GetAllBlocksWithToken("ThreeBody");
  if (blocks.size() > 0) {
    m_active             = true;
    m_ReactionConditions = new TReactionConditions();
    AttachReactionConditions();
    if (!RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
      RootOutput::getInstance()->GetTree()->Branch(
          "ReactionConditions", "TReactionConditions", &m_ReactionConditions);

    for (unsigned int i = 0; i < blocks.size(); i++) {
      if (blocks[i]->HasToken("Type"))
        m_type = blocks[i]->GetString("Type");

      if (blocks[i]->HasToken("Beam"))
        m_Beam = blocks[i]->GetString("Beam");

      if (blocks[i]->HasToken("Target"))
        m_TargetNucleus = blocks[i]->GetString("Target");

      if (blocks[i]->HasToken("CompoundNucleusEnergy"))
        m_CompoundNucleusEnergy
            = blocks[i]->GetDouble("CompoundNucleusEnergy", "MeV");

      if (blocks[i]->HasToken("CompoundNucleusStateWidth"))
        m_CompoundNucleusStateWidth
            = blocks[i]->GetDouble("CompoundNucleusStateWidth", "MeV");

      if (blocks[i]->HasToken("IntermediaryNucleus"))
        m_IntermediaryNucleus = blocks[i]->GetString("IntermediaryNucleus");

      if (blocks[i]->HasToken("IntermediaryState"))
        m_IntermediaryState = blocks[i]->GetDouble("IntermediaryState", "MeV");

      if (blocks[i]->HasToken("IntermediaryStateWidth"))
        m_IntermediaryStateWidth
            = blocks[i]->GetDouble("IntermediaryStateWidth", "MeV");

      if (blocks[i]->HasToken("Light1"))
        m_Light1 = blocks[i]->GetString("Light1");

      if (blocks[i]->HasToken("Heavy"))
        m_Heavy = blocks[i]->GetString("Heavy");

      if (blocks[i]->HasToken("Light2"))
        m_Light2 = blocks[i]->GetString("Light2");

      if (blocks[i]->HasToken("ExcitationEnergyHeavy"))
        m_ExcitationEnergyHeavy
            = blocks[i]->GetDouble("ExcitationEnergyHeavy", "MeV");

      if (blocks[i]->HasToken("ShootLight1"))
        m_ShootLight1 = blocks[i]->GetInt("ShootLight1");

      if (blocks[i]->HasToken("ShootLight2"))
        m_ShootLight2 = blocks[i]->GetInt("ShootLight2");

      if (blocks[i]->HasToken("ShootHeavy"))
        m_ShootHeavy = blocks[i]->GetInt("ShootHeavy");

      m_UserEnergyCS = false;
      if (blocks[i]->HasToken("UserEnergyCS")) {
        m_UserEnergyCS = blocks[i]->GetInt("UserEnergyCS");

        cout << m_UserEnergyCS << endl;

        if (blocks[i]->HasToken("EnergyCSPath"))
          m_EnergyCSPath = blocks[i]->GetString("EnergyCSPath");

        if (blocks[i]->HasToken("EnergyCSName"))
          m_EnergyCSName = blocks[i]->GetString("EnergyCSName");

        if (m_UserEnergyCS) {
          m_FileEnergyCS = new TFile(m_EnergyCSPath.c_str(), "open");
          m_EnergyCS
              = (TF1*)m_FileEnergyCS->FindObjectAny(m_EnergyCSName.c_str());
          m_EnergyCSMax = m_EnergyCS->GetMaximum();
          m_EnergyCS->SetNpx(1e6);
          m_FileEnergyCS->Close();
        }
      } else {
        cout << "ERROR: check your input file formatting \033[0m" << endl;
        exit(1);
      }
    }
    m_BeamName = NPL::ChangeNameToG4Standard(m_Beam);
    m_N1       = NPL::Particle(m_Beam);
    m_N2       = NPL::Particle(m_TargetNucleus);

    m_N3 = NPL::Particle(m_Light1);
    m_N4 = NPL::Particle(m_IntermediaryNucleus);

    m_N5                     = NPL::Particle(m_Light2);
    m_N6                     = NPL::Particle(m_Heavy);
    m_end                    = false;
    m_previousE              = -1000;
    fIntermediaryStateEnergy = GetIntermediaryStateEnergy(
        m_IntermediaryState, m_IntermediaryStateWidth);
  } else
    m_active = false;
}

////////////////////////////////////////////////////////////////////////////////
G4bool NPS::ThreeBody::IsApplicable(const G4ParticleDefinition& particleType) {
  NPL::InputParser input(NPOptionManager::getInstance()->GetReactionFile());

  vector<NPL::InputBlock*> blocks = input.GetAllBlocksWithToken("ThreeBody");
  if (!m_active)
    return false;
  if (!m_active)
    return false;
  else {
    static std::string particleName;
    particleName = particleType.GetParticleName();
    if (particleName.find(m_BeamName) != std::string::npos
        && blocks.size() > 0) {
      return true;
    }
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
G4bool NPS::ThreeBody::ModelTrigger(const G4FastTrack& fastTrack) {
  static bool    shoot        = false;
  static double  rand         = 0;
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();
  G4ThreeVector  V            = PrimaryTrack->GetMomentum().unit();
  G4ThreeVector  P            = fastTrack.GetPrimaryTrackLocalPosition();
  G4VSolid*      solid        = fastTrack.GetPrimaryTrack()
                        ->GetVolume()
                        ->GetLogicalVolume()
                        ->GetSolid();
  double in    = solid->DistanceToOut(P, V);
  double out   = solid->DistanceToOut(P, -V);
  double ratio = in / (out + in);

  if (out == 0) { // first step of current event
    rand             = G4RandFlat::shoot();
    m_PreviousLength = m_StepSize;
    m_PreviousEnergy = PrimaryTrack->GetKineticEnergy();
    // Clear Previous Event
    m_ReactionConditions->Clear();
    shoot = true;
  } else if (in == 0) { // last step
    return true;
  }

  // If the condition is met, the event is generated
  if (ratio < rand) {
    // Reset the static for next event
    if (shoot) {
      // && m_Reaction.IsAllowed(PrimaryTrack->GetKineticEnergy())) {
      shoot = false;
      return true;
    } else {
      return false;
    }
  }

  // Record the situation of the current step
  // so it can be used in the next one
  if (!PrimaryTrack->GetStep()->IsLastStepInVolume()) {
    m_PreviousLength = PrimaryTrack->GetStep()->GetStepLength();
    m_PreviousEnergy = PrimaryTrack->GetKineticEnergy();
  }

  return false;
}

////////////////////////////////////////////////////////////////////////////////
void NPS::ThreeBody::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) {
  // Get the track info
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();
  G4ThreeVector  pdirection   = PrimaryTrack->GetMomentum().unit();
  G4ThreeVector  localdir     = fastTrack.GetPrimaryTrackLocalDirection();

  G4ThreeVector worldPosition = PrimaryTrack->GetPosition();
  G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();

  double energy = PrimaryTrack->GetKineticEnergy();
  double time   = PrimaryTrack->GetGlobalTime();

  /////////////////////////
  /// Beam Parameters /////
  /////////////////////////
  m_ReactionConditions->SetBeamParticleName(
      PrimaryTrack->GetParticleDefinition()->GetParticleName());

  G4ThreeVector U(1, 0, 0);
  G4ThreeVector V(0, 1, 0);
  G4ThreeVector ZZ(0, 0, 1);
  m_ReactionConditions->SetBeamEmittanceTheta(
      PrimaryTrack->GetMomentumDirection().theta() / deg);
  m_ReactionConditions->SetBeamEmittancePhi(
      PrimaryTrack->GetMomentumDirection().phi() / deg);
  m_ReactionConditions->SetBeamEmittanceThetaX(
      PrimaryTrack->GetMomentumDirection().angle(U) / deg);
  m_ReactionConditions->SetBeamEmittancePhiY(
      PrimaryTrack->GetMomentumDirection().angle(V) / deg);

  // Randomize within the step
  // Assume energy loss is linear within the step
  // Assume no scattering
  double rand   = G4RandFlat::shoot();
  double length = rand * (m_PreviousLength);
  energy -= (1 - rand) * (m_PreviousEnergy - energy);
  G4ThreeVector ldir = pdirection;
  ldir *= length;
  localPosition = localPosition - ldir;

  m_ReactionConditions->SetBeamReactionEnergy(energy);
  m_ReactionConditions->SetVertexPositionX(localPosition.x());
  m_ReactionConditions->SetVertexPositionY(localPosition.y());
  m_ReactionConditions->SetVertexPositionZ(localPosition.z());

  // Set the end of the step conditions
  fastStep.SetPrimaryTrackFinalKineticEnergyAndDirection(0, pdirection);
  fastStep.SetPrimaryTrackFinalPosition(worldPosition);
  fastStep.SetTotalEnergyDeposited(0);
  fastStep.SetPrimaryTrackFinalTime(time);
  fastStep.KillPrimaryTrack();
  fastStep.SetPrimaryTrackPathLength(0.0);

  double Ecm = sqrt(m_N1.Mass() * m_N1.Mass() + m_N2.Mass() * m_N2.Mass()
                    + 2 * m_N2.Mass() * (energy + m_N1.Mass()))
               - m_N1.Mass() - m_N2.Mass();
  if (m_UserEnergyCS == 0) {
    Ecm = m_CompoundNucleusEnergy;
  }

  bool Allowed = false;
  if (Ecm < m_IntermediaryState)
    Ecm = 0;
  if (m_UserEnergyCS) {
    double cs_value = m_EnergyCS->Eval(Ecm) / m_EnergyCSMax;
    double rand_cs  = G4RandFlat::shoot();
    if (rand_cs < cs_value) {
      Allowed = true;
    }
  } else if (Ecm > m_IntermediaryState) {
    Allowed = true;
  }

  if (Allowed) {
    /////////////////////////////////////////////////
    //////Define the kind of particle to shoot////////
    //////////////////////////////////////////////////

    // Build the decaying particle four momenta vector:
    double BeamEnergy = energy;

    G4ThreeVector BeamMomentumG4V = pdirection;

    // Before reaction
    double m1 = m_N1.Mass();
    double m2 = m_N2.Mass();

    // After first emission
    double m3 = m_N3.Mass();

    double IntermediaryStateEnergy = 0.;
    double SpSn = m_N4.GetBindingEnergy() - m_N6.GetBindingEnergy();
    if (m_IntermediaryStateWidth == 0) {
      IntermediaryStateEnergy = m_IntermediaryState;
    } else {
      while (true) {
        IntermediaryStateEnergy = fIntermediaryStateEnergy->GetRandom();
        if (IntermediaryStateEnergy > SpSn && IntermediaryStateEnergy < Ecm)
          break;
      }
    }

    m_N4.SetExcitationEnergy(IntermediaryStateEnergy);
    double m4 = m_N4.Mass();

    // After second emission
    double m5 = m_N5.Mass();
    double m6 = m_N6.Mass();

    std::vector<NPL::Particle> DaughterNuclei;
    DaughterNuclei.push_back(m_N3);
    DaughterNuclei.push_back(m_N5);
    DaughterNuclei.push_back(m_N6);

    double T1 = BeamEnergy;
    if (m_UserEnergyCS == 0) {
      T1 = (pow((Ecm + m_N1.Mass() + m_N2.Mass()), 2.)
            - (m_N1.Mass() * m_N1.Mass() + m_N2.Mass() * m_N2.Mass()))
               / (2. * m_N2.Mass())
           - m_N1.Mass();
    }

    TVector3 p1(BeamMomentumG4V.x(), BeamMomentumG4V.y(), BeamMomentumG4V.z());

    double gamma1 = T1 / m1 + 1.;
    double Beta1  = sqrt(1. - 1. / (gamma1 * gamma1));
    p1.SetMag(gamma1 * m1 * Beta1);
    TLorentzVector LV_1(p1.x(), p1.y(), p1.z(), T1 + m1);

    TVector3       p2(0, 0, 0);
    TLorentzVector LV_2(p2.x(), p2.y(), p2.z(), m2);

    TLorentzVector LV_CN = LV_1 + LV_2;

    /////////////////////////
    // Sequential emission //
    /////////////////////////

    if (m_type == "sequential") {
      // First particle emission
      double       T3, T4 = 0;
      unsigned int nDecay   = 2;
      double*      masses34 = new double[nDecay];
      masses34[0]           = m3 / 1000.;
      masses34[1]           = m4 / 1000.;
      double* T34           = new double[nDecay];
      double* theta34       = new double[nDecay];
      double* phi34         = new double[nDecay];

      GetDecay(LV_CN.Vect(), LV_CN.E(), 2, masses34, T34, theta34, phi34);
      TVector3 p3(0, 0, 1);
      p3.SetTheta(theta34[0]);
      p3.SetPhi(phi34[0]);
      T3            = T34[0];
      double gamma3 = T3 / m3 + 1;
      double Beta3  = sqrt(1 - 1 / (gamma3 * gamma3));
      p3.SetMag(gamma3 * m3 * Beta3);
      G4ThreeVector p3G4 = G4ThreeVector(p3.x(), p3.y(), p3.z());

      TVector3 p4(0, 0, 1);
      p4.SetTheta(theta34[1]);
      p4.SetPhi(phi34[1]);
      T4            = T34[1];
      double gamma4 = T4 / m4 + 1;
      double Beta4  = sqrt(1 - 1 / (gamma4 * gamma4));
      p4.SetMag(gamma4 * m4 * Beta4);
      TLorentzVector LV_4(p4.x(), p4.y(), p4.z(), T4 + m4);

      // Second particle emission
      double T5, T6 = 0;
      nDecay           = 2;
      double* masses56 = new double[nDecay];
      masses56[0]      = m5 / 1000.;
      masses56[1]      = m6 / 1000.;
      double* T56      = new double[nDecay];
      double* theta56  = new double[nDecay];
      double* phi56    = new double[nDecay];
      GetDecay(LV_4.Vect(), LV_4.E(), 2, masses56, T56, theta56, phi56);
      TVector3 p5(0, 0, 1);
      p5.SetTheta(theta56[0]);
      p5.SetPhi(phi56[0]);
      T5            = T56[0];
      double gamma5 = T5 / m5 + 1;
      double Beta5  = sqrt(1 - 1 / (gamma5 * gamma5));
      p5.SetMag(gamma5 * m5 * Beta5);
      G4ThreeVector p5G4 = G4ThreeVector(p5.x(), p5.y(), p5.z());

      TVector3 p6(0, 0, 1);
      p6.SetTheta(theta56[1]);
      p6.SetPhi(phi56[1]);
      T6            = T56[1];
      double gamma6 = T6 / m6 + 1;
      double Beta6  = sqrt(1 - 1 / (gamma6 * gamma6));
      p6.SetMag(gamma6 * m6 * Beta6);
      G4ThreeVector p6G4 = G4ThreeVector(p6.x(), p6.y(), p6.z());

      G4ParticleDefinition* Light1
          = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(
              m_N3.GetZ(), m_N3.GetA(), 0);

      G4ParticleDefinition* Light2
          = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(
              m_N5.GetZ(), m_N5.GetA(), 0);

      G4ParticleDefinition* Heavy
          = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(
              m_N6.GetZ(), m_N6.GetA(), m_ExcitationEnergyHeavy);

      TLorentzVector LVproton1(p3.x(), p3.y(), p3.z(), T3 + m3);
      TLorentzVector LVproton2(p5.x(), p5.y(), p5.z(), T5 + m5);
      TLorentzVector LVHeavy(p6.x(), p6.y(), p6.z(), T6 + m6);

      TLorentzVector All = LVproton1 + LVproton2 + LVHeavy;

      // Emitt secondary
      if (m_ShootLight1) {
        G4DynamicParticle particleLight1(Light1, p3G4.unit(), T3);
        fastStep.CreateSecondaryTrack(particleLight1, localPosition, time);
      }

      if (m_ShootLight2) {
        G4DynamicParticle particleLight2(Light2, p5G4.unit(), T5);
        fastStep.CreateSecondaryTrack(particleLight2, localPosition, time);
      }
      if (m_ShootHeavy) {
        G4DynamicParticle particleHeavy(Heavy, p6G4.unit(), T6);
        fastStep.CreateSecondaryTrack(particleHeavy, localPosition, time);
      }
      ///////////////////////////////////////
      ///// Emitted particle Parameters /////
      ///////////////////////////////////////
      m_ReactionConditions->SetParticleName(Light1->GetParticleName());
      m_ReactionConditions->SetKineticEnergy(T3);
      m_ReactionConditions->SetTheta(theta34[0] / deg);
      if ((phi34[0] + pi) / deg > 360)
        m_ReactionConditions->SetPhi((phi34[0] - pi) / deg);
      else
        m_ReactionConditions->SetPhi((phi34[0] + pi) / deg);

      m_ReactionConditions->SetParticleName(Light2->GetParticleName());
      m_ReactionConditions->SetKineticEnergy(T5);
      m_ReactionConditions->SetTheta(theta56[0] / deg);
      if ((phi56[0] + pi) / deg > 360)
        m_ReactionConditions->SetPhi((phi56[0] - pi) / deg);
      else
        m_ReactionConditions->SetPhi((phi56[0] + pi) / deg);

      m_ReactionConditions->SetParticleName(Heavy->GetParticleName());
      m_ReactionConditions->SetKineticEnergy(T6);
      m_ReactionConditions->SetTheta(theta56[1] / deg);
      if ((phi56[1] + pi) / deg > 360)
        m_ReactionConditions->SetPhi((phi56[1] - pi) / deg);
      else
        m_ReactionConditions->SetPhi((phi56[1] + pi) / deg);
    }
  }

  // Reinit for next event
  m_PreviousEnergy = 0;
  m_PreviousLength = 0;
}

void NPS::ThreeBody::GetDecay(TVector3 p1, double Etot1, unsigned int nDecay,
                              double* masses, double* T, double* theta,
                              double* phi) {
  TGenPhaseSpace PhaseSpace;
  TLorentzVector NucleiToDecay(p1.x() / 1000., p1.y() / 1000., p1.z() / 1000.,
                               Etot1 / 1000.);

  if (!PhaseSpace.SetDecay(NucleiToDecay, nDecay, masses)) {
    cout << "FORBIDDEN PHASE SPACE CHECK ENERGIES" << endl;
    for (unsigned int i = 0; i < nDecay; i++) {
      T[i]     = (PhaseSpace.GetDecay(i)->E() - masses[i]) * 1000.;
      theta[i] = PhaseSpace.GetDecay(i)->Theta();
      phi[i]   = PhaseSpace.GetDecay(i)->Phi();
    }
    exit(1);
  } else {
    PhaseSpace.Generate();
    for (unsigned int i = 0; i < nDecay; i++) {
      T[i]     = (PhaseSpace.GetDecay(i)->E() - masses[i]) * 1000.;
      theta[i] = PhaseSpace.GetDecay(i)->Theta();
      phi[i]   = PhaseSpace.GetDecay(i)->Phi();
    }
  }
}

TF1* NPS::ThreeBody::GetIntermediaryStateEnergy(double Er, double Width) {

  TF1* f
      = new TF1("f", "1. / (pow((x - [0]), 2.) + pow(([1] / 2.), 2.))", 0, 100);
  f->SetParameter(0, Er);
  f->SetParameter(1, Width);
  f->SetNpx(1e6);
  return f;
}
