#include "TrackRecording.hh"

#include "NPFunction.h"
#include "RootOutput.h"

TrackRecording::TrackRecording() {
    m_TrackInfo = new TTrackInfo();

    // Reasssigned the branch address
    if (RootOutput::getInstance()->GetTree()->FindBranch("TrackInfo"))
        RootOutput::getInstance()->GetTree()->SetBranchAddress("TrackInfo", &m_TrackInfo);
    else
        RootOutput::getInstance()->GetTree()->Branch("TrackInfo", "TTrackInfo", &m_TrackInfo);
}

TrackRecording::~TrackRecording() { delete m_TrackInfo; }

TrackRecording* TrackRecording::GetInstance() {
    static TrackRecording instance;
    return &instance;
}

void TrackRecording::Clear() { m_TrackInfo->Clear(); }

void TrackRecording::RecordTrack(const G4Track* aTrack) {
    G4String ParticleName = aTrack->GetParticleDefinition()->GetParticleName();
    G4ThreeVector Momentum = aTrack->GetMomentum();
    G4ThreeVector Position = aTrack->GetPosition();

    G4double KineticEnergy = aTrack->GetDynamicParticle()->GetKineticEnergy();
    G4double Mass = aTrack->GetDynamicParticle()->GetMass();
    G4double Charge = aTrack->GetDynamicParticle()->GetCharge();

    G4double c_light = 299.792458;  // To go from T.m to MeV/e
    G4double Brho = sqrt(KineticEnergy * KineticEnergy + 2 * KineticEnergy * Mass) / (c_light * Charge);

    if (ParticleName != "e-" && ParticleName != "e+")
        m_TrackInfo->SetParticleName(NPL::ChangeNameFromG4Standard(ParticleName));
    else
        m_TrackInfo->SetParticleName(ParticleName);
    m_TrackInfo->SetKineticEnergy(KineticEnergy);
    m_TrackInfo->SetTheta(Momentum.theta() * 180. / M_PI);
    m_TrackInfo->SetPhi(Momentum.phi() * 180. / M_PI);
    m_TrackInfo->SetMass(Mass);
    m_TrackInfo->SetCharge(Charge);
    m_TrackInfo->SetZ(aTrack->GetParticleDefinition()->GetAtomicNumber());
    m_TrackInfo->SetA(aTrack->GetParticleDefinition()->GetAtomicMass());
    m_TrackInfo->SetBrho(Brho);
    m_TrackInfo->SetTime(aTrack->GetGlobalTime());
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

    m_TrackInfo->SetVolumeName(aTrack->GetVolume()->GetName());
    m_TrackInfo->SetIndex(aTrack->GetTrackID());
}