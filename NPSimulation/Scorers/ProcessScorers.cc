/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: P. MORFOUACE  contact address: pierre.morfouace2@cea.fr  *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Scorer specific to the porcesses occuring in the simulation              *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * The scorer hold the processes name                                        *
 *****************************************************************************/
#include "ProcessScorers.hh"
#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"
#include "RootOutput.h"
using namespace ProcessScorers;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Process::PS_Process(G4String name, G4int depth) : G4VPrimitiveScorer(name, depth) {
  // m_NestingLevel = NestingLevel;
  auto tree = RootOutput::getInstance()->GetTree();
  tree->Branch("ProcScor_TrackID", &t_trackID);
  tree->Branch("ProcScor_ParentID", &t_parentID);
  tree->Branch("ProcScor_Particle", &t_particlename);
  tree->Branch("ProcScor_Volume", &t_volumename);
  tree->Branch("ProcScor_Process", &t_processname);
  tree->Branch("ProcScor_Time", &t_processtime);
  tree->Branch("ProcScor_Efficiency", &t_efficiency);
  tree->Branch("ProcScor_El1_InEl_2", &t_FC_process);
  tree->Branch("ProcScor_Pos_Z", &t_pos_Z);
  tree->Branch("ProcScor_KinEn", &t_kinetic_energy);
  tree->Branch("ProcScor_EnergyDeposit", &t_energy_deposit);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Process::~PS_Process() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_Process::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  // Contain Process Name
  G4String processname;
  double pos_z;
  if (aStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL) {
    processname = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    t_processname.push_back(processname);
    if(processname!="Transportation"){
      pos_z = aStep->GetPostStepPoint()->GetPosition()[2];
      t_pos_Z.push_back(pos_z);
    }
    
    //cout << "Process :" << processname << endl;
    t_processtime.push_back(aStep->GetPostStepPoint()->GetGlobalTime());
  }

  G4int trackID;
  G4int parentID;
  G4String volumename, particlename;
  //G4double step_length;
  G4double track_kineE;
  G4double track_Edep;
  
  if(aStep->GetTrack()->GetNextVolume() != 0){
    trackID = aStep->GetTrack()->GetTrackID();
    parentID = aStep->GetTrack()->GetParentID();
    volumename = aStep->GetTrack()->GetVolume()->GetName();
    particlename = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
    //cout << "Particule :" << particlename << endl;
    //step_length = aStep->GetTrack()->GetStepLength();
    track_kineE = aStep->GetTrack()->GetKineticEnergy();
    track_Edep = aStep->GetTotalEnergyDeposit();

    t_trackID.push_back(trackID);
    t_parentID.push_back(parentID);
    t_particlename.push_back(particlename);
    t_volumename.push_back(volumename);
    t_kinetic_energy.push_back(track_kineE);
    t_energy_deposit.push_back(track_Edep);

    //Looking after the neutrons 
    if(particlename=="neutron"){

      //If it's an original neutron
      if(parentID==0){
        //Checking what reactions it prefers between El and InEl 
        if(processname=="hadElastic")t_FC_process.push_back(5);
        if(processname=="neutronInelastic")t_FC_process.push_back(6);
        //Checking wall efficiency by discriminating whether it has interactied or not 
        if(HasBeenTracked.find(trackID)==HasBeenTracked.end()){//first interaction
          if(processname!="Transportation"){
            // # SPECIAL NEBULA PROJECT #
            if(volumename=="NebulaModule" || volumename=="NebulaVeto"){
              //if the hit happens in Nebula, it counts as efficiency
              HasBeenTracked.insert(trackID);
              t_efficiency.push_back(1);
            }
            else{
              //the purpose of this is for dicrimination, 0 being "didn't interact", 1 "interacted", 5 "interacted elsewhere"
              HasBeenTracked.insert(trackID);
              t_efficiency.push_back(5);
            }
            //Checking if first interaction is El or InEl
            if(processname=="hadElastic")t_FC_process.push_back(1);
            if(processname=="neutronInelastic")t_FC_process.push_back(2);
          }
        }
      }
      //Otherwise it's a neutron coming from secondary reactions
      else{
        if(processname=="hadElastic")t_FC_process.push_back(9);
        if(processname=="neutronInelastic")t_FC_process.push_back(10);
      }
    }//endif neutron


    if(HasBeenTracked.find(trackID)==HasBeenTracked.end()){
      if(particlename=="gamma"){
        HasBeenTracked.insert(trackID);
        //cout << trackID  << " " << track_kineE << endl;
        t_gamma_energy.push_back(track_kineE);
      }
      if(particlename=="proton"){
        HasBeenTracked.insert(trackID);
        if(track_kineE>0.1){
          t_proton_energy.push_back(track_kineE);
          t_proton_time.push_back(aStep->GetPreStepPoint()->GetGlobalTime());
        }
      }
    }//endif first interaction (other than neutron)
  }//endif track unended
  else{
    trackID = aStep->GetTrack()->GetTrackID();
    if(HasBeenTracked.find(trackID)==HasBeenTracked.end()){
      t_efficiency.push_back(0);
    }
  }

  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Process::Initialize(G4HCofThisEvent*) { clear(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Process::EndOfEvent(G4HCofThisEvent*) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Process::clear() {
  t_trackID.clear();
  t_parentID.clear();
  t_particlename.clear();
  t_volumename.clear();
  t_processname.clear();
  t_processtime.clear();
  t_gamma_energy.clear();
  t_proton_energy.clear();
  t_proton_time.clear();
  t_FC_process.clear();
  t_efficiency.clear();
  t_pos_Z.clear();
  t_kinetic_energy.clear();
  t_energy_deposit.clear();

  HasBeenTracked.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Process::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Process::PrintAll() {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

