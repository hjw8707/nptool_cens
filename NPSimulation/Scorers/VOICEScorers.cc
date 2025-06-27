/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Sunghan Bae     contact adresse: shbae2703@gmail.com     *  
 *                                                                           *
 * Creation Date  : Feb. 2024                                                *
 * Last update    : 22/02/24                                                 *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Scorer for VOICE ionization chamber					     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *****************************************************************************/

#include "VOICEScorers.hh"
#include "G4UnitsTable.hh"
#include "G4VPhysicalVolume.hh"
using namespace VOICEScorers;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PS_Electrode::PS_Electrode(G4String name, G4int Level,
        G4double OpenX, G4double OpenY, G4double PCBZ, G4int depth)
    : G4VPrimitiveScorer(name, depth), HCID(-1) {
	
        m_Level = Level;
	m_OpenX = OpenX;
	m_OpenY = OpenY;
	m_PCBZ = PCBZ;
	m_position = G4ThreeVector(-1000,-1000,-1000);
	m_detectorNumber = -1;
	m_index = -1;
	multiness =0;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
PS_Electrode::~PS_Electrode() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_Electrode::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    //Energy, position[3], Time, PIDindex
    G4double* EnergyAndTime = new G4double[12];

    EnergyAndTime[4] = aStep->GetTrack()->GetGlobalTime();

    m_position = aStep->GetTrack()->GetPosition();

    EnergyAndTime[1] = m_position.x();
    EnergyAndTime[2] = m_position.y();
    EnergyAndTime[3] = m_position.z();
    EnergyAndTime[5] = m_position.theta();
    EnergyAndTime[6] = m_position.phi();

    EnergyAndTime[7] = aStep->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
    EnergyAndTime[8] = aStep->GetTrack()->GetParticleDefinition()->GetAtomicMass();
    EnergyAndTime[9] = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

    EnergyAndTime[0] = aStep->GetTotalEnergyDeposit();


    m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10;
    map<long, long>::iterator itmnmap = multinessmap.find(m_index);

    if(itmnmap!= multinessmap.end()){ multiness = itmnmap->second;
    }else {multiness = 0;
	    multinessmap[m_index]=multiness;
    }

    m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10+multiness*10000000;
    
    map<G4int, G4double**>::iterator it = EvtMap->GetMap()->find(m_index);
    if(it!= EvtMap->GetMap()->end()){
	G4double* dummy = *(it->second);
	if(abs(dummy[3]-EnergyAndTime[3])<m_PCBZ){
		EnergyAndTime[0] += dummy[0];	
	}else{
		multiness+=1;
    		m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10;
		multinessmap[m_index]=multiness;
    		m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10+multiness*10000000;
	}
    }
    EvtMap->set(m_index,EnergyAndTime);

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Electrode::Initialize(G4HCofThisEvent* HCE) {
	EvtMap= new NPS::HitsMap<G4double*>(GetMultiFunctionalDetector()->GetName(), GetName());
 	if(HCID < 0) HCID = GetCollectionID(0);	
	HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
	}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Electrode::EndOfEvent(G4HCofThisEvent*) {
	multinessmap.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PS_Electrode::clear() {
	std::map<G4int, G4double**>::iterator MapIterator;
	for (MapIterator = EvtMap->GetMap()->begin();
			MapIterator != EvtMap->GetMap()->end(); MapIterator++){

		delete *(MapIterator->second);
	}
	EvtMap->clear();
	multinessmap.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Electrode::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Electrode::PrintAll() {
	G4cout<<"  MultiFunctionalDet	" <<detector->GetName() << G4endl;
	G4cout<<"  PrimitiveScorer	" <<GetName() << G4endl;
	G4cout<<"  Number of entries	" <<EvtMap->entries() << G4endl;
}
//....oooOO0OOooo.......oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//For end plane

PS_END::PS_END(G4String name, G4int Level,
        G4double ENDthick, G4int depth)
    : G4VPrimitiveScorer(name, depth), HCID(-1) {
        m_Level = Level;
	m_ENDthick = ENDthick;
	m_position = G4ThreeVector(-1000,-1000,-1000);
	m_detectorNumber = -1;
	m_index = -1;
	multiness=0;
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
PS_END::~PS_END() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_END::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    //Energy, position[3], Time, PIDindex
    G4double* EnergyAndTime = new G4double[12];

    EnergyAndTime[4] = aStep->GetTrack()->GetGlobalTime();

    m_position = aStep->GetTrack()->GetPosition();

    EnergyAndTime[1] = m_position.x();
    EnergyAndTime[2] = m_position.y();
    EnergyAndTime[3] = m_position.z();
    EnergyAndTime[5] = m_position.theta();
    EnergyAndTime[6] = m_position.phi();
    
    EnergyAndTime[7] = aStep->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
    EnergyAndTime[8] = aStep->GetTrack()->GetParticleDefinition()->GetAtomicMass();
    EnergyAndTime[9] = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

    EnergyAndTime[0] = aStep->GetTotalEnergyDeposit();

    m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10;
    map<long, long>::iterator itmnmap = multinessmap.find(m_index);

    if(itmnmap!= multinessmap.end()){ multiness = itmnmap->second;
    }else {multiness = 0;
		multinessmap[m_index]=multiness;
    }

    m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10+multiness*10000000;
    
    map<G4int, G4double**>::iterator it = EvtMap->GetMap()->find(m_index);
    if(it!= EvtMap->GetMap()->end()){
	G4double* dummy = *(it->second);
	if(abs(dummy[3]-EnergyAndTime[3])<m_ENDthick){
		EnergyAndTime[0] += dummy[0];	
	}else{
		multiness+=1;
    		m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10;
		multinessmap[m_index]=multiness;
    		m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10+multiness*10000000;
	}
    }
    EvtMap->set(m_index,EnergyAndTime);
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_END::Initialize(G4HCofThisEvent* HCE) {
	EvtMap= new NPS::HitsMap<G4double*>(GetMultiFunctionalDetector()->GetName(), GetName());
 	if(HCID < 0) HCID = GetCollectionID(0);	
	HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
	}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_END::EndOfEvent(G4HCofThisEvent*) {
	multinessmap.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PS_END::clear() {
	std::map<G4int, G4double**>::iterator MapIterator;
	for (MapIterator = EvtMap->GetMap()->begin();
			MapIterator != EvtMap->GetMap()->end(); MapIterator++){

		delete *(MapIterator->second);
	}
	EvtMap->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_END::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_END::PrintAll() {
	G4cout<<"  MultiFunctionalDet	" <<detector->GetName() << G4endl;
	G4cout<<"  PrimitiveScorer	" <<GetName() << G4endl;
	G4cout<<"  Number of entries	" <<EvtMap->entries() << G4endl;
}
//....oooOO0OOooo.......oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//For Gas region

PS_Gas::PS_Gas(G4String name, vector<double> zpos, G4double StepSize, G4int Level, 
        G4int depth)
    : G4VPrimitiveScorer(name, depth), HCID(-1) {
        m_Level = Level;
	m_position = G4ThreeVector(-1000,-1000,-1000);
	m_detectorNumber = -1;
	m_index = -1;
	multiness=0;
	m_zpos = zpos;
    	m_i_stage =0;
	m_StepSize = StepSize;
//	m_checker.clear();
//	cout<<m_zpos.size()<<endl;
//	for(int i=0; i<m_zpos.size(); i++){
//		m_checker.push_back(0);
//	}
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
PS_Gas::~PS_Gas() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_Gas::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    //Energy, position[3], Time, PIDindex
    
    G4double* EnergyAndTime = new G4double[17];

    EnergyAndTime[7] = aStep->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
    EnergyAndTime[8] = aStep->GetTrack()->GetParticleDefinition()->GetAtomicMass();
    EnergyAndTime[9] = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
    m_position = aStep->GetTrack()->GetPosition();

    m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10;
    map<long, long>::iterator itmnmap = multinessmap.find(m_index);
    map<long, int>::iterator itgmmap = gasmulmap.find(m_index);

    if(itmnmap!= multinessmap.end()){ 
	    multiness = itmnmap->second;
    }else {
	    multiness = 0;
	    multinessmap[m_index]=multiness;
    }

    if(itgmmap!= gasmulmap.end()){ 
	    m_i_stage = itgmmap->second;
    }else {
	    m_i_stage = (lower_bound(m_zpos.begin(),m_zpos.end(),m_position.z())-m_zpos.begin());
	    gasmulmap[m_index]=m_i_stage;
    }

    m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10+multiness*10000000+1;
    map<G4int, G4double**>::iterator it = EvtMap->GetMap()->find(m_index);

    if(it==EvtMap->GetMap()->end()){
    
	 if(aStep->IsFirstStepInVolume()){
	 EnergyAndTime[0]= aStep->GetTrack()->GetKineticEnergy();
	 EnergyAndTime[1] = m_position.x();
	 EnergyAndTime[2] = m_position.y();
	 EnergyAndTime[3] = m_position.z();
	 EnergyAndTime[4] = aStep->GetTrack()->GetGlobalTime();
	 EnergyAndTime[5] = m_position.theta();
	 EnergyAndTime[6] = m_position.phi();
	 EnergyAndTime[10] = -1;
	 EnergyAndTime[11] = -1;
	 EnergyAndTime[12] = -1;
	 EnergyAndTime[13] = -1;
	 EnergyAndTime[14] = -1;
	 EnergyAndTime[15] = -1;
	 EnergyAndTime[16] = -1;
   
	 EvtMap->set(m_index,EnergyAndTime);
	// endflag=false;
   	 }else if(m_i_stage>0){
	
     	 m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10+(multiness-1)*10000000+1;
   	 it = EvtMap->GetMap()->find(m_index);
	 G4double* dummy = *(it->second);

	 EnergyAndTime[0] = dummy[10];
	 EnergyAndTime[1] = dummy[11];
	 EnergyAndTime[2] = dummy[12];
	 EnergyAndTime[3] = dummy[13];
	 EnergyAndTime[4] = dummy[14];
	 EnergyAndTime[5] = dummy[15];
	 EnergyAndTime[6] = dummy[16];
	 EnergyAndTime[10] = -1;
	 EnergyAndTime[11] = -1;
	 EnergyAndTime[12] = -1;
	 EnergyAndTime[13] = -1;
	 EnergyAndTime[14] = -1;
	 EnergyAndTime[15] = -1;
	 EnergyAndTime[16] = -1;
   
     	 m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10+(multiness)*10000000+1;
	 EvtMap->set(m_index,EnergyAndTime);
	 }
      }else if(abs(m_zpos[m_i_stage]-m_position.z()) < m_StepSize){

	    G4double* dummy = *(it->second);
//	    if(EnergyAndTime[10]>0) return true;	
	    EnergyAndTime[0] = dummy[0];
	    EnergyAndTime[1] = dummy[1];
	    EnergyAndTime[2] = dummy[2];
	    EnergyAndTime[3] = dummy[3];
	    EnergyAndTime[4] = dummy[4];
	    EnergyAndTime[5] = dummy[5];
	    EnergyAndTime[6] = dummy[6];

	    EnergyAndTime[10] = aStep->GetTrack()->GetKineticEnergy();
	    EnergyAndTime[11] = m_position.x();
	    EnergyAndTime[12] = m_position.y();
	    EnergyAndTime[13] = m_position.z();
	    EnergyAndTime[14] = aStep->GetTrack()->GetGlobalTime();
	    EnergyAndTime[15] = m_position.theta();
	    EnergyAndTime[16] = m_position.phi();

	    EvtMap->set(m_index,EnergyAndTime);
	    multiness+=1;
	    m_index = EnergyAndTime[7]*10000+EnergyAndTime[8]*10;
	    multinessmap[m_index]=multiness;
	    m_i_stage+=1;
	    gasmulmap[m_index]=m_i_stage;

      }

return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Gas::Initialize(G4HCofThisEvent* HCE) {
	EvtMap= new NPS::HitsMap<G4double*>(GetMultiFunctionalDetector()->GetName(), GetName());
	if(HCID < 0) HCID = GetCollectionID(0);	
	HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Gas::EndOfEvent(G4HCofThisEvent*) {
//	cout<<m_i_stage<<" m_i_stage"<<endl;
	multinessmap.clear();
	gasmulmap.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PS_Gas::clear() {
	std::map<G4int, G4double**>::iterator MapIterator;
	for (MapIterator = EvtMap->GetMap()->begin();
			MapIterator != EvtMap->GetMap()->end(); MapIterator++){

		delete *(MapIterator->second);
	}
//	m_i_stage=0;
	EvtMap->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Gas::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Gas::PrintAll() {
	G4cout<<"  MultiFunctionalDet	" <<detector->GetName() << G4endl;
	G4cout<<"  PrimitiveScorer	" <<GetName() << G4endl;
	G4cout<<"  Number of entries	" <<EvtMap->entries() << G4endl;
}
//....oooOO0OOooo.......oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
