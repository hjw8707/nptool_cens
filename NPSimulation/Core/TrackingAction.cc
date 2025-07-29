#include "TrackingAction.hh"

#include "G4RunManager.hh"
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
    if (m_record_track && aTrack->GetParentID() < m_cut_parent_id) TrackRecording::GetInstance()->RecordTrack(aTrack);
    //////////////////////////////////////////////////////////////////////////////////////////////
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack) {}