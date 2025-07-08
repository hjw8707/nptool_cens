#include "TrackingAction.hh"

#include "G4RunManager.hh"
#include "NPFunction.h"
#include "NPOptionManager.h"
#include "RootOutput.h"

TrackingAction::TrackingAction() {
    m_record_track = NPOptionManager::getInstance()->GetRecordTrack();
    m_cut_parent_id = NPOptionManager::getInstance()->GetCutParentID();
    m_tree = RootOutput::getInstance()->GetTree();

    m_First = true;
    ParticleName = "";
    KineticEnergy = -1000;
    Theta = -1000;
    Phi = -1000;
    Mass = -1000;
    Charge = -1000;
    Z = -1000;
    A = -1000;
    Time = -1000;
    Momentum.setX(-1000);
    Momentum.setY(-1000);
    Momentum.setZ(-1000);
    Position.setX(-1000);
    Position.setY(-1000);
    Position.setZ(-1000);
    LastStepNumber = -1000;
    VolumeName = "";
    nClear = 0;
    TrackID = -1000;

    m_TrackInfo = new TTrackInfo();
    AttachTrackInfo();
    if (!RootOutput::getInstance()->GetTree()->FindBranch("TrackInfo"))
        RootOutput::getInstance()->GetTree()->Branch("TrackInfo", "TTrackInfo", &m_TrackInfo);
}

TrackingAction::~TrackingAction() {}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack) {
    if (aTrack->GetParentID() > m_cut_parent_id) {
        m_record_track = false;
    }
    if (m_record_track) {
        TrackRecording(aTrack);
    }
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::AttachTrackInfo() {
    // Reasssigned the branch address
    if (RootOutput::getInstance()->GetTree()->FindBranch("TrackInfo"))
        RootOutput::getInstance()->GetTree()->SetBranchAddress("TrackInfo", &m_TrackInfo);
}

// FIXME Still underdeveloppement
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::TrackRecording(const G4Track* aTrack) {
    if (m_First) m_TrackInfo->Clear();
    m_First = false;

    if (eventID < G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()) {
        m_TrackInfo->Clear();
        nClear++;
    }

    ParticleName = aTrack->GetParticleDefinition()->GetParticleName();
    KineticEnergy = aTrack->GetDynamicParticle()->GetKineticEnergy();
    Mass = aTrack->GetDynamicParticle()->GetMass();
    Charge = aTrack->GetDynamicParticle()->GetCharge();
    Z = aTrack->GetParticleDefinition()->GetAtomicNumber();
    A = aTrack->GetParticleDefinition()->GetAtomicMass();

    Momentum = aTrack->GetMomentum();
    Time = aTrack->GetGlobalTime();
    Theta = Momentum.theta() * 180. / M_PI;
    Phi = Momentum.phi() * 180. / M_PI;
    Position = aTrack->GetPosition();
    VolumeName = aTrack->GetVolume()->GetName();
    TrackID = aTrack->GetTrackID();

    double c_light = 299.792458;  // To go from T.m to MeV/e
    double Brho = sqrt(KineticEnergy * KineticEnergy + 2 * KineticEnergy * Mass) / (c_light * Charge);

    m_TrackInfo->SetKineticEnergy(KineticEnergy);
    m_TrackInfo->SetTheta(Theta);
    m_TrackInfo->SetPhi(Phi);
    m_TrackInfo->SetMass(Mass);
    m_TrackInfo->SetCharge(Charge);
    m_TrackInfo->SetZ(Z);
    m_TrackInfo->SetA(A);
    m_TrackInfo->SetBrho(Brho);
    m_TrackInfo->SetTime(Time);
    TVector3 Mom;
    Mom.SetX(Momentum.x());
    Mom.SetY(Momentum.y());
    Mom.SetZ(Momentum.z());

    m_TrackInfo->SetMomentum(Mom);

    m_TrackInfo->SetMomentumX(Momentum.x());
    m_TrackInfo->SetMomentumY(Momentum.y());
    m_TrackInfo->SetMomentumZ(Momentum.z());

    m_TrackInfo->SetPositionX(Position.x());
    m_TrackInfo->SetPositionY(Position.y());
    m_TrackInfo->SetPositionZ(Position.z());

    m_TrackInfo->SetVolumeName(VolumeName);
    m_TrackInfo->SetIndex(eventID + TrackID);

    if (ParticleName != "e-" && ParticleName != "e+")
        m_TrackInfo->SetParticleName(NPL::ChangeNameFromG4Standard(ParticleName));
    else
        m_TrackInfo->SetParticleName(ParticleName);
    eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
}
