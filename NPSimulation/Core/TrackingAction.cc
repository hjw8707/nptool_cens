#include "TrackingAction.hh"

#include "NPOptionManager.h"
#include "TrackRecording.hh"

TrackingAction::TrackingAction() {
  m_record_track = NPOptionManager::getInstance()->GetRecordTrack();
  m_cut_parent_id = NPOptionManager::getInstance()->GetCutParentID();
}

TrackingAction::~TrackingAction() {}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack) {
  //////////////////////////////////////////////////////////////////////////////////////////////
  // recording track if --record-track option is on and aTrack->GetParentID() < m_cut_parent_id
  if (m_record_track && aTrack->GetParentID() < m_cut_parent_id)
    TrackRecording::GetInstance()->RecordTrack(aTrack);
  //////////////////////////////////////////////////////////////////////////////////////////////
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack) {
  AddSimulatedParticle(aTrack->GetParticleDefinition()->GetParticleName());
}

void TrackingAction::AddSimulatedParticle(const std::string& particle_name) {
  m_simulated_particles.insert(particle_name);
}

bool TrackingAction::IsSimulatedParticle(const std::string& particle_name) const {
  return m_simulated_particles.find(particle_name) != m_simulated_particles.end();
}

std::set<std::string> TrackingAction::GetSimulatedParticles() const { return m_simulated_particles; }