#ifndef TrackRecording_h
#define TrackRecording_h

#include "G4Track.hh"
#include "TTrackInfo.h"

class TrackRecording {
   private:
    TrackRecording();
    ~TrackRecording();

   public:
    static TrackRecording* GetInstance();
    void Clear();
    void RecordTrack(const G4Track* aTrack);

   private:
    TTrackInfo* m_TrackInfo;
};

#endif