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
 *  File old the scorer specific to the Sharc Detector                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This new type of scorer is aim to become the standard for DSSD,SSSD and   *
 * PAD detector (any STARK Detector)                                       *
 *****************************************************************************/
#include "STARKScorers.hh"
#include "G4UnitsTable.hh"

using namespace STARKSCORERS ;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// X6
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_STARK_X6::PS_STARK_X6(G4String name, G4int Level,
			 G4double stripPlaneX, G4double stripPlaneY,
			 G4int nFrontStrip, G4int nBackStrip,
			 G4int depth)
  : G4VPrimitiveScorer(name, depth), HCID(-1) {
  m_stripPlaneX = stripPlaneX;
  m_stripPlaneY = stripPlaneY;
  m_nFrontStrip = nFrontStrip;
  m_nBackStrip  = nBackStrip;
  m_frontStripPitch = m_stripPlaneX / m_nFrontStrip; // X direction for front
  m_backStripPitch  = m_stripPlaneY / m_nBackStrip;  // Y direction for back
  m_Level           = Level;
  
  m_Position = G4ThreeVector(-1000,-1000,-1000);
  m_detectorNumber = -1;
  m_frontStripNumber = -1;
  m_backStripNumber = -1;
  m_Index = -1;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_STARK_X6::~PS_STARK_X6(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_STARK_X6::ProcessHits(G4Step* aStep, G4TouchableHistory*){
  // contain E_DW, E_UP, E_Bot, Time, DetNbr, and StripWidth
  G4double* EnergyAndTime = new G4double[12];

  EnergyAndTime[3] = aStep->GetPreStepPoint()->GetGlobalTime();

  m_detectorNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_Level);
  m_Position  = aStep->GetPreStepPoint()->GetPosition();

  // Interaction coordinates (used to fill the InteractionCoordinates branch)
  EnergyAndTime[7]  = m_Position.x();
  EnergyAndTime[8]  = m_Position.y();
  EnergyAndTime[9]  = m_Position.z();
  EnergyAndTime[10]  = m_Position.theta();
  EnergyAndTime[11] = m_Position.phi();

  m_Position = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(m_Position);

  m_frontStripNumber = (int)((-m_Position.x() + m_stripPlaneX / 2.) / m_frontStripPitch ) + 1  ; // reverse direction for X (from 1)
  m_backStripNumber = (int)((m_Position.y() + m_stripPlaneY / 2.) / m_backStripPitch ) + 1  ; // forward direction for Y (from 1)

  ////////////////////////////////////////////////////////////
  // Not active area: wrong strip numbers -> not hit
  if (m_frontStripNumber < 1 || m_frontStripNumber > m_nFrontStrip ||
      m_backStripNumber < 1  || m_backStripNumber > m_nBackStrip)
    return false;
  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // Energy deposition for resistive strip
  // The front energy is divided in two depending on the position
  // position along the resistive strip
  double P = (m_Position.y())/(0.5*m_stripPlaneY); // P = -1 ~ 1
  ////////////////////////////////////////////////////////////
  
  // Downstream Energy
  EnergyAndTime[0] = aStep->GetTotalEnergyDeposit()*(1+P)*0.5; // downsteam collects MORE energy when P is positive
  // Upstream Energy
  EnergyAndTime[1] = aStep->GetTotalEnergyDeposit()-EnergyAndTime[0];
  // back side Energy
  EnergyAndTime[2] = aStep->GetTotalEnergyDeposit(); // just sum

  //Rare case where particle is close to edge of silicon plan
  EnergyAndTime[4] = m_detectorNumber;
  EnergyAndTime[5] = m_frontStripNumber;
  EnergyAndTime[6] = m_backStripNumber;

  m_Index = m_detectorNumber * 1e3 + m_frontStripNumber * 1e6 + m_backStripNumber * 1e9;
  // Check if the particle has interact before, if yes, add up the energies.
  map<G4int, G4double**>::iterator it = EvtMap->GetMap()->find(m_Index);
  if (it != EvtMap->GetMap()->end()){
    G4double* dummy = *(it->second);
    EnergyAndTime[0] += dummy[0];
    EnergyAndTime[1] += dummy[1];
    EnergyAndTime[2] += dummy[2];}
  EvtMap->set(m_Index, EnergyAndTime);
  return true;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_X6::Initialize(G4HCofThisEvent* HCE){
    EvtMap = new NPS::HitsMap<G4double*>(GetMultiFunctionalDetector()->GetName(), GetName());
    if (HCID < 0) HCID = GetCollectionID(0);
    HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_X6::EndOfEvent(G4HCofThisEvent*){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_X6::clear(){
    std::map<G4int, G4double**>::iterator MapIterator;
    for (MapIterator = EvtMap->GetMap()->begin() ;
	 MapIterator != EvtMap->GetMap()->end() ; MapIterator++){
        delete *(MapIterator->second); }
    EvtMap->clear();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_X6::DrawAll(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_X6::PrintAll(){
    G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
    G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
    G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// BB10
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_STARK_BB10::PS_STARK_BB10(G4String name, G4int Level,
			     G4double geoX, G4double geoY,
			     G4int nStrip, G4int depth)
  : G4VPrimitiveScorer(name, depth), HCID(-1) {
  m_geoX       = geoX;
  m_geoY       = geoY;
  m_nStrip     = nStrip;
  m_stripPitch = m_geoX / m_nStrip; // X direction for front
  m_Level      = Level;
  
  m_Position   = G4ThreeVector(-1000,-1000,-1000);
  m_detectorNumber = -1;
  m_stripNumber = -1;
  m_Index = -1;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_STARK_BB10::~PS_STARK_BB10(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_STARK_BB10::ProcessHits(G4Step* aStep, G4TouchableHistory*){
  // contain E_DW, E_UP, E_Bot, Time, DetNbr, and StripWidth
  G4double* EnergyAndTime = new G4double[10];

  EnergyAndTime[2] = aStep->GetPreStepPoint()->GetGlobalTime();
  m_detectorNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_Level);
  m_Position  = aStep->GetPreStepPoint()->GetPosition();

  // Interaction coordinates (used to fill the InteractionCoordinates branch)
  EnergyAndTime[5]  = m_Position.x();
  EnergyAndTime[6]  = m_Position.y();
  EnergyAndTime[7]  = m_Position.z();
  EnergyAndTime[8]  = m_Position.theta();
  EnergyAndTime[9]  = m_Position.phi();

  m_Position = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(m_Position);

  m_stripNumber = (int)((-m_Position.x() + m_geoX / 2.) / m_stripPitch ) + 1; // reverse direction for X
  ////////////////////////////////////////////////////////////
  // Not active area: wrong strip numbers -> not hit
  if (m_stripNumber < 1 || m_stripNumber > m_nStrip)
    return false;
  if (m_Position.y() < -m_geoY/2. || m_Position.y() > m_geoY/2.) // no strip for Y just check the real position
    return false;
  ////////////////////////////////////////////////////////////
  
  EnergyAndTime[0] = aStep->GetTotalEnergyDeposit(); // Junction side
  EnergyAndTime[1] = aStep->GetTotalEnergyDeposit(); // Ohmic side

  //Rare case where particle is close to edge of silicon plan
  EnergyAndTime[3] = m_detectorNumber;
  EnergyAndTime[4] = m_stripNumber;

  m_Index = m_detectorNumber * 1e3 + m_stripNumber * 1e6;

  // Check if the particle has interact before, if yes, add up the energies.
  map<G4int, G4double**>::iterator it = EvtMap->GetMap()->find(m_Index);
  if (it != EvtMap->GetMap()->end()){
    G4double* dummy = *(it->second);
    EnergyAndTime[0] += dummy[0];
    EnergyAndTime[1] += dummy[1];}
  EvtMap->set(m_Index, EnergyAndTime);
  return true;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_BB10::Initialize(G4HCofThisEvent* HCE){
    EvtMap = new NPS::HitsMap<G4double*>(GetMultiFunctionalDetector()->GetName(), GetName());
    if (HCID < 0) HCID = GetCollectionID(0);
    HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_BB10::EndOfEvent(G4HCofThisEvent*){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_BB10::clear(){
    std::map<G4int, G4double**>::iterator MapIterator;
    for (MapIterator = EvtMap->GetMap()->begin() ;
	 MapIterator != EvtMap->GetMap()->end() ; MapIterator++){
        delete *(MapIterator->second); }
    EvtMap->clear();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_BB10::DrawAll(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_BB10::PrintAll(){
    G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
    G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
    G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// QQQ5
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_STARK_QQQ5::PS_STARK_QQQ5(G4String name, G4int Level,
			     G4double inR, G4double outR, G4double phi0, G4double phi1,
			     G4int NRStrip,G4int NAStrip, G4int depth)
  : G4VPrimitiveScorer(name, depth), HCID(-1) {
  m_inR       = inR;      m_outR      = outR;
  m_phi0      = phi0;     m_phi1      = phi1;
  m_NRStrip   = NRStrip;  m_NAStrip   = NAStrip;
  m_AStripPitch = (m_outR - m_inR) / m_NAStrip;
  m_RStripPitch = (m_phi1 - m_phi0) / m_NRStrip;
  m_Level      = Level;
  m_Position   = G4ThreeVector(-1000,-1000,-1000);
  m_detectorNumber = -1;
  m_radialStripNumber = -1;
  m_annularStripNumber = -1;
  m_Index = -1;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_STARK_QQQ5::~PS_STARK_QQQ5(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_STARK_QQQ5::ProcessHits(G4Step* aStep, G4TouchableHistory*){
  // contain E_DW, E_UP, E_Bot, Time, DetNbr, and StripWidth
  G4double* EnergyAndTime = new G4double[11];

  EnergyAndTime[2] = aStep->GetPreStepPoint()->GetGlobalTime();
  m_detectorNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(m_Level);
  m_Position  = aStep->GetPreStepPoint()->GetPosition();

  // Interaction coordinates (used to fill the InteractionCoordinates branch)
  EnergyAndTime[ 6]  = m_Position.x();
  EnergyAndTime[ 7]  = m_Position.y();
  EnergyAndTime[ 8]  = m_Position.z();
  EnergyAndTime[ 9]  = m_Position.theta();
  EnergyAndTime[10]  = m_Position.phi();

  m_Position = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(m_Position);
  
  m_annularStripNumber = (int)((m_outR - m_Position.rho()) / m_AStripPitch ) + 1 ; // reverse in R(Rho), from 1
  m_radialStripNumber = (int)((m_Position.phi() - m_phi0) / m_RStripPitch) + 1 ; // from 1

  ////////////////////////////////////////////////////////////
  // Not active area: wrong strip numbers -> not hit
  if (m_radialStripNumber < 1  || m_radialStripNumber > m_NRStrip ||
      m_annularStripNumber < 1 || m_annularStripNumber > m_NAStrip)
    return false;
  ////////////////////////////////////////////////////////////

  
  EnergyAndTime[0] = aStep->GetTotalEnergyDeposit(); // Junction side
  EnergyAndTime[1] = aStep->GetTotalEnergyDeposit(); // Ohmic side

  //Rare case where particle is close to edge of silicon plan
  EnergyAndTime[3] = m_detectorNumber;
  EnergyAndTime[4] = m_annularStripNumber;
  EnergyAndTime[5] = m_radialStripNumber;

  m_Index = m_detectorNumber * 1e3 + m_radialStripNumber * 1e6 + m_annularStripNumber * 1e9;

  // Check if the particle has interact before, if yes, add up the energies.
  map<G4int, G4double**>::iterator it = EvtMap->GetMap()->find(m_Index);
  if (it != EvtMap->GetMap()->end()){
    G4double* dummy = *(it->second);
    EnergyAndTime[0] += dummy[0];
    EnergyAndTime[1] += dummy[1];}
  EvtMap->set(m_Index, EnergyAndTime);
  return true;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_QQQ5::Initialize(G4HCofThisEvent* HCE){
    EvtMap = new NPS::HitsMap<G4double*>(GetMultiFunctionalDetector()->GetName(), GetName());
    if (HCID < 0) HCID = GetCollectionID(0);
    HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_QQQ5::EndOfEvent(G4HCofThisEvent*){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_QQQ5::clear(){
    std::map<G4int, G4double**>::iterator MapIterator;
    for (MapIterator = EvtMap->GetMap()->begin() ;
	 MapIterator != EvtMap->GetMap()->end() ; MapIterator++){
        delete *(MapIterator->second); }
    EvtMap->clear();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_QQQ5::DrawAll(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_STARK_QQQ5::PrintAll(){
    G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl ;
    G4cout << " PrimitiveScorer " << GetName() << G4endl               ;
    G4cout << " Number of entries " << EvtMap->entries() << G4endl     ;}
