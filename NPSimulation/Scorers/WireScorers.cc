/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey ANNE  contact address: anne@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  : July 2024                                                *
 * Last update    : January 2025                                             *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Scorer for wires in gas detector                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "WireScorers.hh"
#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"
#include "RootOutput.h"
#include "TMath.h"

#include "TVector3.h"
using namespace WireScorers;
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
unsigned int WireData::CalculateIndex(const vector<unsigned int>& level) {

  unsigned int size = level.size();
  unsigned int result = 0;
  unsigned int multiplier = 1;
  for (unsigned int i = 0; i < size; i++) {
    result += level[i] * multiplier;
    multiplier *= 1000;
  }
  return result;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
vector<WireData>::iterator WireDataVector::find(const unsigned int& index) {
  for (vector<WireData>::iterator it = m_Data.begin(); it != m_Data.end(); it++) {
    if ((*it).GetIndex() == index)
      return it;
  }
  return m_Data.end();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Wire::PS_Wire(G4String name, vector<G4int> NestingLevel, G4int TotalNumberLayer, G4int NumberWireByLayer, G4double DriftSpeed,std::map<unsigned int , double> map_WireAngle, G4int depth) : G4VPrimitiveScorer(name, depth) {
  m_NestingLevel = NestingLevel;
  m_TotalNumberLayer = TotalNumberLayer;
  m_NumberWireByLayer = NumberWireByLayer;
  m_DriftSpeed = DriftSpeed;
  m_WireAngle = map_WireAngle;
  
  auto tree = RootOutput::getInstance()->GetTree();
  // tree->Branch("Layer_Number", &t_LayerNumber);
  // tree->Branch("Wire_Number", &t_WireNumber); 
  // tree->Branch("Wire_Time", &t_Time);
  // tree->Branch("Edge", &t_Edge);

  G4String PosInWireX = name + "_PosInWireX";
  G4String PosInWireY = name + "_PosInWireY";
  G4String PosInWireZ = name + "_PosInWireZ";

  tree->Branch(PosInWireX, &t_PosInWireX);
  tree->Branch(PosInWireY, &t_PosInWireY);
  tree->Branch(PosInWireZ, &t_PosInWireZ);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_Wire::~PS_Wire() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_Wire::ProcessHits(G4Step* aStep, G4TouchableHistory*) {

  // Contain information as many copy number as nested volume
  unsigned int mysize = m_NestingLevel.size();
   
  G4String particlename = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
 
  //Solid information
  static G4VSolid* solid;
  solid = aStep->GetTrack()->GetTouchable()->GetSolid();
  
  G4ThreeVector WirePosition =  aStep->GetTrack()->GetTouchable()->GetTranslation(); //wire's center position
  //cout << solid->GetName() << "  WirePos X = " << WirePosition[0] << "  WirePos Y = " << WirePosition[1] << "  WirePos Z = " << WirePosition[2] << endl;
  const G4RotationMatrix* WireRotation = aStep->GetTrack()->GetTouchable()->GetRotation();

  // we only keep information about charged particles
  if (particlename != "neutron") {
    
    //Detector number and Wire copy number
    t_DetectorNumber = aStep->GetTrack()->GetTouchable()->GetCopyNumber(1); //Wires are contained in a box
    int CopyNbr;
    CopyNbr  =  aStep->GetTrack()->GetTouchable()->GetCopyNumber(0);
    t_WireNumber = (CopyNbr-1)%m_NumberWireByLayer; //-1 because copynbr is the ID contained in the xml file and ID begins at 1.
    t_LayerNumber = (floor((CopyNbr-1)/m_NumberWireByLayer));
 
    //cout <<  "solid name = " << solid->GetName() <<  "  particle name = " << particlename << " Copy number = " << CopyNbr << " layer number = " << t_LayerNumber << " wire number =  " << t_WireNumber << endl;

    t_Level.clear();

    G4ThreeVector PrePosWire =  aStep->GetPreStepPoint()->GetPosition(); //lab position when particle enters the cylinder 
    G4ThreeVector PostPosWire =  aStep->GetPostStepPoint()->GetPosition(); //lab position when particle exists the cylinder
   
    G4ThreeVector InOutDir = PostPosWire - PrePosWire; //particle's direction vector 
    G4ThreeVector PointC = PrePosWire + 0.5*InOutDir; //closest point to the cylinder center on the particle's path

    //cout << "  presposwire =  " << PrePosWire << "  postposwire =  " << PostPosWire << "  PointC =  "<< PointC << endl;
    
    t_PosInWireX.push_back(PointC[0]); //leading egde
    t_PosInWireX.push_back(-10000); // trailing edge
    
    t_PosInWireY.push_back(PointC[1]);
    t_PosInWireY.push_back(-10000);
    
    t_PosInWireZ.push_back(PointC[2]);
    t_PosInWireZ.push_back(-10000);

    G4ThreeVector DirectionWireZ; // wire direction vector z' in the lab  x,y,z vectors basis 
    DirectionWireZ = WireRotation->rowZ();

    G4ThreeVector PointA = WirePosition; // first point on the symmetry axis of the cylinder
    G4ThreeVector PointB = WirePosition + DirectionWireZ;  // second point on the symmetry axis of the cylinder

    //cout << " CenterPointWire =  " << PointA << "  CenterPointWire + shift =  " << PointB << endl;
    
    G4ThreeVector AB = PointB - PointA;
    G4ThreeVector AC = PointC - PointA;

    G4ThreeVector CrossABAC = AB.cross(AC);
    G4double MagCrossABAC = CrossABAC.mag();
    G4double MagAB = AB.mag();
  
    G4double CD = MagCrossABAC/MagAB; // DriftLength
    t_DrifLength = CD;

    
    t_Time = (t_DrifLength/m) / (m_DriftSpeed*299792458);///!!!NEEDS TO BE FIXED BEFORE USE!!!///
    
    for (unsigned int i = 0; i < mysize; i++) {
      t_Level.push_back(aStep->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber(m_NestingLevel[i]));
    }
    
    m_Data.Set(t_Time*1e9 , t_Level, t_LayerNumber, t_WireNumber,1, t_DrifLength, t_DetectorNumber); // leading egde
    m_Data.Set((t_Time+(t_Time/10))*1e9 , t_Level, t_LayerNumber, t_WireNumber,0, t_DrifLength+(t_DrifLength*100), t_DetectorNumber);// trailing edge

  }
 
  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Wire::Initialize(G4HCofThisEvent*) { clear();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Wire::EndOfEvent(G4HCofThisEvent*) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Wire::clear() {
  m_Data.clear();
  t_Level.clear();

  t_PosInWireX.clear();
  t_PosInWireY.clear();
  t_PosInWireZ.clear();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Wire::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_Wire::PrintAll() {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
