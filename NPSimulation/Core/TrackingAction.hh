#ifndef TrackingAction_h
#define TrackingAction_h 1

// G4 header defining G4 types
#include "globals.hh"

// NPL
#include "TTrackInfo.h"

// G4 header
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4UserTrackingAction.hh"

// Root
#include "TTree.h"

class TrackingAction : public G4UserTrackingAction {
   public:
    TrackingAction();
    ~TrackingAction();

   public:
    void PreUserTrackingAction(const G4Track* aTrack);
    void PostUserTrackingAction(const G4Track* aTrack);

   public:
    void TrackRecording(const G4Track* aTrack);

    std::string ParticleName;
    double KineticEnergy;
    double Theta;
    double Phi;
    double Mass;
    double Charge;
    double A;
    double Z;
    double Time;
    G4ThreeVector Momentum;
    G4ThreeVector Position;
    std::string VolumeName;
    int nClear;
    int eventID;
    int TrackID;

   private:  // tree
    TTree* m_tree;
    bool m_First;
    unsigned int LastStepNumber;

    // Host the Initial conditions TObject
    TTrackInfo* m_TrackInfo;

   private:
    bool m_record_track;
    int m_cut_parent_id;

   public:
    // Attach the TrackInfo object to the tree
    void AttachTrackInfo();
};
#endif