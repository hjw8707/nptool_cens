/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : February 2013                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  File old the scorer to record Hit energy,time and position               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "ATOMXScorers.hh"
#include "G4UnitsTable.hh"
using namespace ATOMXScorers ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
vector<ATOMXData>::iterator ATOMXDataVector::find(const unsigned int& index){
  for(vector<ATOMXData>::iterator it= m_Data.begin()  ; it !=m_Data.end() ; it++){
    if((*it).GetIndex()==index)
      return it;
  }
  return m_Data.end();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_ATOMX::PS_ATOMX(G4String name, int depth)  :G4VPrimitiveScorer(name, depth){
  m_Level = depth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_ATOMX::ProcessHits(G4Step* aStep, G4TouchableHistory*){
  static G4StepPoint* point;
  point = aStep->GetPreStepPoint();
  t_Position = point->GetPosition();
  
  m_DataVector.Set(
          aStep->GetTrack()->GetTrackID(),
          aStep->GetTotalEnergyDeposit(),
          point->GetGlobalTime(),
          point->GetPosition().x(),
          point->GetPosition().y(),
          point->GetPosition().z(),
          aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding(),
          aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy());

  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_ATOMX::Initialize(G4HCofThisEvent*){
  // Clear is called by EventAction
  clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_ATOMX::EndOfEvent(G4HCofThisEvent*){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_ATOMX::clear(){
  m_DataVector.clear();
}
