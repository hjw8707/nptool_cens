/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: ValerianAlcindor  contact address:
 *valcindor@@ikp.tu-darmstadt.de
 *                                                                           *
 * Creation Date  : September 2021                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                  []                                             *
 *  A quite Standard Geant4 EventAction class.                               *
 *  Call the Fill method of the output tree.                                 *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "SteppingAction.hh"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Run.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4UnitsTable.hh"
#include "NPOptionManager.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction() {
    m_cut_parent_id = NPOptionManager::getInstance()->GetCutParentID();
    m_online_n = NPOptionManager::getInstance()->GetOnlineStream();
    m_online_run = -1;
    m_optical_death_stats = (std::getenv("NPS_OPTICAL_DEATH_STATS") != nullptr);
    m_optical_death_total = 0;
}

SteppingAction::~SteppingAction() { WriteOpticalDeathStats("optical_photon_death_stats.txt"); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* step) {
    G4Track* track = step->GetTrack();
    if (track->GetParentID() > m_cut_parent_id) track->SetTrackStatus(fKillTrackAndSecondaries);

    if (m_optical_death_stats && track->GetDefinition()->GetParticleName() == "opticalphoton") {
        if (track->GetTrackStatus() != fAlive) {
            const G4StepPoint* pre = step->GetPreStepPoint();
            const G4StepPoint* post = step->GetPostStepPoint();
            const G4VProcess* process = post ? post->GetProcessDefinedStep() : nullptr;
            const G4VPhysicalVolume* pre_volume = pre ? pre->GetPhysicalVolume() : nullptr;
            const G4VPhysicalVolume* post_volume = post ? post->GetPhysicalVolume() : nullptr;

            const std::string process_name = process ? process->GetProcessName() : "NoProcess";
            const std::string pre_volume_name = pre_volume ? pre_volume->GetName() : "NoPreVolume";
            const std::string post_volume_name = post_volume ? post_volume->GetName() : "OutOfWorld";
            const std::string key = process_name + " | " + pre_volume_name + " -> " + post_volume_name;
            const double step_length = step->GetStepLength();
            const double track_length = track->GetTrackLength();

            OpticalDeathStats& death_stats = m_optical_death_by_key[key];
            ++death_stats.count;
            death_stats.last_step_length += step_length;
            death_stats.track_length += track_length;
            if (track_length > death_stats.max_track_length) death_stats.max_track_length = track_length;
            ++m_optical_death_total;
        }
    }

    // ---- online text streaming for a live web viewer ----------------------------
    if (m_online_n > 0) {
        G4RunManager* rm = G4RunManager::GetRunManager();
        const G4Event* evt = rm->GetCurrentEvent();
        const G4Run* run = rm->GetCurrentRun();
        int eventID = evt ? evt->GetEventID() : 0;
        int runID = run ? run->GetRunID() : 0;
        // New run (a fresh /run/beamOn) → truncate so the file holds only this run.
        if (runID != m_online_run) {
            m_online_run = runID;
            if (m_online.is_open()) m_online.close();
            m_online.open("online_stream.dat", std::ios::out | std::ios::trunc);
            m_online << "# run " << runID
                     << " ; columns: evt trk parent pdg charge x y z edep_MeV" << std::endl;
        }
        if (eventID < m_online_n && m_online.is_open()) {
            // Geant4 internal units are already mm / MeV, so emit raw values.
            const G4ParticleDefinition* pd = track->GetParticleDefinition();
            auto emit = [&](const G4ThreeVector& q, double edep) {
                m_online << eventID << ' ' << track->GetTrackID() << ' ' << track->GetParentID() << ' '
                         << (pd ? pd->GetPDGEncoding() : 0) << ' ' << (pd ? pd->GetPDGCharge() : 0) << ' '
                         << q.x() << ' ' << q.y() << ' ' << q.z() << ' ' << edep << '\n';
            };
            emit(step->GetPreStepPoint()->GetPosition(), step->GetTotalEnergyDeposit());
            // On the track's final step, also emit the END point — otherwise a particle
            // that dies after a long step (e.g. a photon escaping into the vacuum
            // envelope) would end at its last vertex instead of where it actually went.
            if (track->GetTrackStatus() != fAlive)
                emit(step->GetPostStepPoint()->GetPosition(), 0.0);
            m_online.flush();  // flushed so another process can tail it live
        }
    }
}

void SteppingAction::ResetOpticalDeathStats() {
    m_optical_death_by_key.clear();
    m_optical_death_total = 0;
}

void SteppingAction::WriteOpticalDeathStats(const std::string& filename) const {
    if (!m_optical_death_stats) return;

    std::ofstream out(filename.c_str());
    out << "# Optical photon death statistics\n";
    out << "# total " << m_optical_death_total << "\n";
    out << "\n";
    out << "# Optical photon death by process and volume\n";
    out << "# columns: count fraction mean_last_step_mm mean_track_length_mm max_track_length_mm process | pre_volume -> "
           "post_volume\n";

    for (const auto& entry : m_optical_death_by_key) {
        const OpticalDeathStats& stats = entry.second;
        const double fraction = m_optical_death_total > 0 ? double(stats.count) / double(m_optical_death_total) : 0;
        const double mean_last_step = stats.count > 0 ? stats.last_step_length / double(stats.count) : 0;
        const double mean_track_length = stats.count > 0 ? stats.track_length / double(stats.count) : 0;
        out << stats.count << " " << fraction << " " << mean_last_step << " " << mean_track_length << " "
            << stats.max_track_length << " " << entry.first << "\n";
    }
}
