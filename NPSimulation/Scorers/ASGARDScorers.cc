/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

#include <algorithm>
#include "ASGARDScorers.hh"
#include "G4UnitsTable.hh"
using namespace std;
using namespace CLHEP;
using namespace ASGARDSCORERS;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_ASGARD::PS_ASGARD(G4String name,G4int ,G4int depth)
  :G4VPrimitiveScorer(name, depth),HCID(-1){
  m_Position = G4ThreeVector(-1000,-1000,-1000);
  m_CloverNumber  = -1;
  m_CrystalNumber = -1;
  m_SegmentNumber = -1;
  m_Index = -1 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_ASGARD::~PS_ASGARD(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_ASGARD::ProcessHits(G4Step* aStep, G4TouchableHistory*){
  // contain Energy, time, segment, crystal, clover
  G4double* Infos = new G4double[10];
  if (aStep->GetTotalEnergyDeposit() > 0){
    Infos[0] = aStep->GetTotalEnergyDeposit(); // energy per step (not per segment yet)
    Infos[1] = aStep->GetPreStepPoint()->GetGlobalTime(); // beta decay time (E+17) is included!

    // get step information at different logical volume levels 0: segment, 1: crystal, 2: clover
    m_SegmentNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
    m_CrystalNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
    m_CloverNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2);
    
    //    G4cout<<m_CloverNumber<<" "<<m_CrystalNumber<<" "<<m_SegmentNumber<<" "<<Infos[0]<<" "<<Infos[1]<<G4endl;
    
    m_Position  = aStep->GetPreStepPoint()->GetPosition();
    // Interaction coordinates (true info, used to fill the InteractionCoordinates branch)
    Infos[2] = m_Position.x();
    Infos[3] = m_Position.y();
    Infos[4] = m_Position.z();
    Infos[5] = m_Position.theta();
    Infos[6] = m_Position.phi();

    // define index for energy summation in segment
    m_Index = m_CloverNumber * 1e2 + m_CrystalNumber * 1e1 + m_SegmentNumber;
    
    // Sum up all hits in the same segment according to m_Index
    map<G4int, G4double**>::iterator it;
    it= EvtMap->GetMap()->find(m_Index);
    if(it!=EvtMap->GetMap()->end()){
      G4double* dummy = *(it->second); // dummy Infos for energy sum
      Infos[0]+=dummy[0];
    }
    EvtMap->set(m_Index, Infos);
  }
  return TRUE;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_ASGARD::Initialize(G4HCofThisEvent* HCE){
  EvtMap = new NPS::HitsMap<G4double*>(GetMultiFunctionalDetector()->GetName(), GetName());
  if (HCID < 0) {
    HCID = GetCollectionID(0);
  }
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_ASGARD::EndOfEvent(G4HCofThisEvent*){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_ASGARD::clear(){
  std::map<G4int, G4double**>::iterator    MapIterator;
  for (MapIterator = EvtMap->GetMap()->begin() ; MapIterator != EvtMap->GetMap()->end() ; MapIterator++){
    delete *(MapIterator->second);
  }

  EvtMap->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_ASGARD::DrawAll(){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_ASGARD::PrintAll(){
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
  G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;
}
