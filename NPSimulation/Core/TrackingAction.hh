#ifndef TrackingAction_h
#define TrackingAction_h 1

// G4 header defining G4 types
#include "globals.hh"
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

   private:
    bool m_record_track;
    int m_cut_parent_id;
};
#endif