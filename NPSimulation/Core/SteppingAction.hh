#ifndef SteppingAction_h
#define SteppingAction_h 1
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : January 2021                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  A quite Standard Geant4 SteppingAction class.                            *
 *  Call the Fill method of the output tree.                                 *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// G4 header defining G4 types
#include "globals.hh"

// STL
#include <fstream>
#include <map>
#include <string>

// NPL
#include "TTrackInfo.h"

// G4 header
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4UserSteppingAction.hh"

// Root
#include "TTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class SteppingAction : public G4UserSteppingAction {
   public:
    SteppingAction();
    ~SteppingAction();

   public:
    void UserSteppingAction(const G4Step* step);
    void ResetOpticalDeathStats();
    void WriteOpticalDeathStats(const std::string& filename) const;

   private:
    struct OpticalDeathStats {
        long long count = 0;
        double last_step_length = 0;
        double track_length = 0;
        double max_track_length = 0;
    };

    int m_cut_parent_id;
    // Online text streaming (--online-data-streaming N): write the first N events
    // of each run to online_stream.dat (flushed) so a web viewer can read it live.
    int m_online_n;
    int m_online_run;
    std::ofstream m_online;
    bool m_optical_death_stats;
    std::map<std::string, OpticalDeathStats> m_optical_death_by_key;
    long long m_optical_death_total;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
