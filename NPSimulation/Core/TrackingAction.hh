#ifndef TrackingAction_h
#define TrackingAction_h 1

#include <set>
#include <string>

// G4 header
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4UserTrackingAction.hh"

class TrackingAction : public G4UserTrackingAction {
 public:
  TrackingAction();
  ~TrackingAction();

 public:
  void PreUserTrackingAction(const G4Track* aTrack);
  void PostUserTrackingAction(const G4Track* aTrack);

  void AddSimulatedParticle(const std::string& particle_name);
  bool IsSimulatedParticle(const std::string& particle_name) const;
  std::set<std::string> GetSimulatedParticles() const;

 private:
  bool m_record_track;
  int m_cut_parent_id;

  // to write DEDX table, save all simulated particles
  std::set<std::string> m_simulated_particles;
};
#endif