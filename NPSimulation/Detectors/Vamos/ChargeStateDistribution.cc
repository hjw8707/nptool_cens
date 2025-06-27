/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elia Pilotto, Omar Nasr                                  *
 * contact address: pilottoelia@gmail.com, omar.nasr@etu.unicaen.fr          *
 *                                                                           *
 * Creation Date  : September 2021                                           *
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

//C++ libraries
#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>
#include <Randomize.hh>
#include <fstream>
//G4 libraries
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4EmCalculator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4IonTable.hh"
#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
//nptool libraries
#include "NPFunction.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "NPSFunction.hh"
//other
#include "ChargeStateDistribution.hh"
#include "RootOutput.h"
#include "TLorentzVector.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
NPS::ChargeStateDistribution::ChargeStateDistribution(G4String modelName, G4Region* envelope)
  : G4VFastSimulationModel(modelName, envelope) {

    m_StepSize = 5*mm;
}

////////////////////////////////////////////////////////////////////////////////
NPS::ChargeStateDistribution::ChargeStateDistribution(G4String modelName)
  : G4VFastSimulationModel(modelName) {}

////////////////////////////////////////////////////////////////////////////////
NPS::ChargeStateDistribution::~ChargeStateDistribution() {}

////////////////////////////////////////////////////////////////////////////////
G4bool NPS::ChargeStateDistribution::IsApplicable(const G4ParticleDefinition& particleType) {
  if (particleType.GetPDGCharge() < 10) return false;
  else return true;
}

////////////////////////////////////////////////////////////////////////////////
G4bool NPS::ChargeStateDistribution::ModelTrigger(const G4FastTrack& fastTrack) {

  cout << "////////////////// MODEL TRIGGER ///////////////////" << endl;
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();

  G4ThreeVector V = PrimaryTrack->GetMomentum().unit();
  G4ThreeVector P = fastTrack.GetPrimaryTrackLocalPosition();
  G4VSolid* solid = fastTrack.GetPrimaryTrack()->GetVolume()->GetLogicalVolume()->GetSolid();
  double to_exit = solid->DistanceToOut(P, V);
  double to_entrance = solid->DistanceToOut(P, -V);
  bool is_first = (to_entrance == 0);
  bool is_end = (to_exit == 0);

  if (is_first && m_shoot) {
    m_shoot = false;
  }

  if (is_first) {
    m_shoot = true;
  }

  cout << "m_StepSize= " << m_StepSize << endl;
  cout << "to_entrance= " << to_entrance << endl;
  cout << "to_exit= " << to_exit << endl;
  cout << "is_first= " << is_first << endl;
  cout << "m_shoot= " << m_shoot << endl;
   
  cout << "Model trigger Q = " << PrimaryTrack->GetDynamicParticle()->GetCharge() << endl;
  return m_shoot;
}

////////////////////////////////////////////////////////////////////////////////
void NPS::ChargeStateDistribution::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) {

  m_shoot = false;
  cout << "////////////////// DO IT ///////////////////" << endl;
  G4ThreeVector localPosition  = fastTrack.GetPrimaryTrackLocalPosition();
  G4ThreeVector localDirection = fastTrack.GetPrimaryTrackLocalDirection();
  const G4Track* track = fastTrack.GetPrimaryTrack();

  G4double time                  = track->GetGlobalTime();
  G4ThreeVector pdirection       = track->GetMomentum().unit();
  G4ThreeVector Momentum         = track->GetMomentum();
  G4ThreeVector worldPosition    = track->GetPosition();
  G4double originalKineticEnergy = track->GetKineticEnergy();
  
  G4int Z = track->GetParticleDefinition()->GetAtomicNumber();
  G4int A = track->GetParticleDefinition()->GetAtomicMass();
  G4double Q = track->GetDynamicParticle()->GetCharge();

  //G4DynamicParticle* dynamicParticle = const_cast<G4DynamicParticle*>(fastTrack.GetPrimaryTrack()->GetDynamicParticle());
  //dynamicParticle->SetCharge(Q + 5);  // Modify the particle's charge

  //G4double newQ = track->GetDynamicParticle()->GetCharge();
  
  //const G4DynamicParticle* dynamicParticle = fastTrack.GetPrimaryTrack()->GetDynamicParticle();
  //G4ParticleDefinition* particleDef = const_cast<G4ParticleDefinition*>(dynamicParticle->GetDefinition());

  /*
  fastStep.SetPrimaryTrackFinalKineticEnergy(originalKineticEnergy);
  fastStep.ProposePrimaryTrackFinalMomentumDirection(localDirection,true);
  fastStep.SetPrimaryTrackFinalPosition(localPosition+G4ThreeVector(0,0,1),true);
  fastStep.ProposePrimaryTrackFinalTime(time);
  fastStep.SetTotalEnergyDeposited(0);
  fastStep.ProposePrimaryTrackPathLength(m_StepSize);
  */

  // Set the end of the step conditions
  fastStep.SetPrimaryTrackFinalKineticEnergyAndDirection(0, localDirection);
  fastStep.SetPrimaryTrackFinalPosition(localPosition+G4ThreeVector(0,0,1),true);
  fastStep.SetTotalEnergyDeposited(0);
  fastStep.SetPrimaryTrackFinalTime(time); // FIXME
  fastStep.KillPrimaryTrack();
  fastStep.SetPrimaryTrackPathLength(0.0);

  G4int newCharge = Q+5;
  static G4IonTable* IonTable = G4ParticleTable::GetParticleTable()->GetIonTable();
  G4ParticleDefinition* particleDef;
  particleDef = IonTable->GetIon(Z,A);

  G4DynamicParticle dynamicParticle(particleDef,Momentum.unit(),originalKineticEnergy);
  fastStep.CreateSecondaryTrack(dynamicParticle,localPosition+G4ThreeVector(0,0,2),time);

  cout << Q << " " << dynamicParticle.GetCharge() << endl;

    return;
}


