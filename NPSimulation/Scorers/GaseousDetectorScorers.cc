/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : February 2013                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  File old the scorer specific to the Sharc Detector                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This new type of scorer is aim to become the standard for DSSD,SSSD and   *
 * PAD detector (any Silicon Detector)                                       *
 *****************************************************************************/
#include "GaseousDetectorScorers.hh"
#include "G4UnitsTable.hh"
using namespace GaseousDetectorScorers;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
unsigned int
GaseousDetectorData::CalculateIndex(const vector<unsigned int>& level) {

  unsigned int size       = level.size();
  unsigned int result     = 0;
  unsigned int multiplier = 1;
  for (unsigned int i = 0; i < size; i++) {
    result += level[i] * multiplier;
    multiplier *= 1000;
  }
  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
vector<GaseousDetectorData>::iterator
GaseousDetectorDataVector::find(const unsigned int& index) {
  for (vector<GaseousDetectorData>::iterator it = m_Data.begin();
       it != m_Data.end(); it++) {
    if ((*it).GetIndex() == index) // comparison of the current index with the argument index
      return it;                   // return immediatly the iterator
  }
  return m_Data.end();             // return m_Data.end() if index is not found
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_GaseousDetector::PS_GaseousDetector(G4String name, vector<G4int> NestingLevel,
                               G4int depth)
    : G4VPrimitiveScorer(name, depth) {
  m_NestingLevel = NestingLevel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_GaseousDetector::~PS_GaseousDetector() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_GaseousDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*) {

  // Contain Particle Id, Energy, Time and as many copy number as nested volume
  unsigned int mysize = m_NestingLevel.size();
  string particlename = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
  G4ThreeVector posVertex = aStep->GetTrack()->GetVertexPosition();
  double x = posVertex.x();
  double y = posVertex.y();
  double z = posVertex.z();
  int particleA = aStep->GetTrack()->GetParticleDefinition()->GetAtomicMass();
  int particleZ = aStep->GetTrack()->GetParticleDefinition()->GetAtomicNumber();

  int trackID  = aStep->GetTrack()->GetTrackID();
  double step_posZ = aStep->GetPostStepPoint()->GetPosition().z(); // z coordinate at the end of the step.
  t_Energy = aStep->GetTotalEnergyDeposit();  
  t_Time   = 0.5*(aStep->GetPreStepPoint()->GetGlobalTime()+aStep->GetPostStepPoint()->GetGlobalTime());  // DeltaT [ns] since the begining of the simulated event up to the average time of the step
  t_ParticleName.push_back(particlename);
  t_ParticleA.push_back(particleA);
  t_ParticleZ.push_back(particleZ);
  t_TrackID.push_back(trackID);
  t_StepPosZ.push_back(step_posZ);
  t_EnergyLossPerStep.push_back(t_Energy);
  t_StepTime.push_back(t_Time);
  
  t_Level.clear();
  for (unsigned int i = 0; i < mysize; i++) {
    t_Level.push_back(
        aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(
            m_NestingLevel[i]));
  }
  

  // Check if the particle has interact before, if yes, add up the energies and update the track
  vector<GaseousDetectorData>::iterator it;
  it = m_Data.find(GaseousDetectorData::CalculateIndex(t_Level));
  if (it != m_Data.end()) {
    // update the m_Data vector at each step corresponding at the same sensitive volume
    it->Add(t_Energy);
    it->SetParticleName(particlename);
    it->SetParticleA(particleA);
    it->SetParticleZ(particleZ);
    it->SetTrackID(trackID);
    it->SetStepPosZ(step_posZ);
    it->SetEnergyLossPerStep(t_Energy);
    it->SetStepTime(t_Time);
  } 
  else {
    // first step: initialize a new m_Data vector<GaseousDetectorData>
    // t_Time [ns] is the time at first step
    m_Data.Set(t_Energy, t_Time, x, y, z, t_Level, t_ParticleName, t_ParticleA, t_ParticleZ, t_TrackID, t_StepPosZ, t_EnergyLossPerStep, t_StepTime);
  }


  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_GaseousDetector::Initialize(G4HCofThisEvent*) { clear(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_GaseousDetector::EndOfEvent(G4HCofThisEvent*) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_GaseousDetector::clear() {
  m_Data.clear();
  t_Level.clear();
  t_ParticleName.clear();
  t_ParticleA.clear();
  t_ParticleZ.clear();
  t_TrackID.clear();
  t_StepPosZ.clear();
  t_EnergyLossPerStep.clear();
  t_StepTime.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_GaseousDetector::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_GaseousDetector::PrintAll() {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
