/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: hjw8707@gmail.com                        *
 *                                                                           *
 * Creation Date  : 7월 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Plunger simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <cmath>
#include <limits>
#include <sstream>
// G4 Geometry object
#include "G4Box.hh"
#include "G4IntersectionSolid.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4UserLimits.hh"

// G4 sensitive
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"

// G4 various object
#include "G4Colour.hh"
#include "G4FastSimulationManager.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4VFastSimulationModel.hh"
#include "G4VisAttributes.hh"

// NPTool header
#include "BeamReaction.hh"
#include "CalorimeterScorers.hh"
#include "Decay.hh"
#include "InteractionScorers.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"
#include "NPSHitsMap.hh"
#include "Plunger.hh"
#include "RootOutput.h"

// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Plunger_NS {
// Energy and time Resolution
const double EnergyThreshold = 0.1 * MeV;
const double ResoTime = 4.5 * ns;
const double ResoEnergy = 1.0 * MeV;
const double Radius = 50 * mm;
const double Width = 100 * mm;
const double Thickness = 300 * mm;
const string Material = "BC400";
}  // namespace Plunger_NS
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Plunger Specific Method
Plunger::Plunger()
    : m_TargetFound(false),
      m_StopperFound(false),
      m_ChamberFound(false),
      m_TargetR(0),
      m_TargetThickness(0),
      m_TargetMaterial(""),
      m_StopperR(0),
      m_StopperThickness(0),
      m_StopperMaterial(""),
      m_ChamberR(0),
      m_ChamberThickness(0),
      m_ChamberMaterial(""),
      m_ChamberPipeR(0),
      m_ChamberPipeZ0(0),
      m_ChamberPipeZ1(0) {
    m_Event = new TPlungerData();
    m_TargetLV = nullptr;
    m_StopperLV = nullptr;
    m_ChamberLV = nullptr;
    m_ReactionRegionLV = nullptr;

    // RGB Color + Transparency
    m_VisTarget = new G4VisAttributes(G4Colour(0, 1, 0));
    m_VisStopper = new G4VisAttributes(G4Colour(0, 0, 1));
    // 챔버를 회색(0.5, 0.5, 0.5)으로 하고, 투명도를 0.3으로 설정하여 안이 비치게 함
    m_VisChamber = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.3));
}

Plunger::~Plunger() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Plunger::BuildTarget() {
    if (!m_TargetLV) {
        G4Tubs* targetCylinder =
            new G4Tubs("Plunger_TargetCylinder", 0, m_TargetR, m_TargetThickness * 0.5, 0. * deg, 360. * deg);

        G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(m_TargetMaterial);
        m_TargetLV = new G4LogicalVolume(targetCylinder, DetectorMaterial, "logic_Plunger_TargetCylinder", 0, 0, 0);
        m_TargetLV->SetVisAttributes(m_VisTarget);
    }
    return m_TargetLV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Plunger::BuildStopper() {
    if (!m_StopperLV) {
        G4Tubs* stopperCylinder =
            new G4Tubs("Plunger_StopperCylinder", 0, m_StopperR, m_StopperThickness * 0.5, 0. * deg, 360. * deg);
        G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(m_StopperMaterial);
        m_StopperLV = new G4LogicalVolume(stopperCylinder, DetectorMaterial, "logic_Plunger_StopperCylinder", 0, 0, 0);
        m_StopperLV->SetVisAttributes(m_VisStopper);
    }
    return m_StopperLV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Plunger::BuildChamber() {
    if (!m_ChamberLV) {
        // 챔버는 구형 껍데기에 z축 방향으로 앞뒤로 파이프가 연결된 형태로 생성합니다.

        // 구형 챔버 껍데기 생성
        G4double chamberInnerR = m_ChamberR - m_ChamberThickness;
        G4double chamberOuterR = m_ChamberR;
        G4Sphere* sphereChamber = new G4Sphere("Plunger_Chamber_Sphere", chamberInnerR, chamberOuterR, 0. * deg,
                                               360. * deg, 0. * deg, 180. * deg);

        // 파이프가 연결되는 앞뒤(z+/-) 방향에 파이프 외경만큼의 구멍을 구 껍데기에 뚫어줍니다.
        // 구멍을 뚫기 위해 파이프 외경과 동일한 반지름의 실린더를 구 껍데기와 Boolean Subtraction 연산

        // 파이프 내경(실제 구멍 반지름)
        G4double pipeHoleR = m_ChamberPipeR;
        G4double pipeHoleHalfLength = m_ChamberR * 2;  // 챔버 반지름보다는 충분히 길게

        // 구멍용 실린더
        G4Tubs* pipeHole =
            new G4Tubs("Plunger_Chamber_PipeHole", 0, pipeHoleR, pipeHoleHalfLength, 0. * deg, 360. * deg);

        // 구 껍데기에서 구멍 빼기
        G4SubtractionSolid* sphereChamberWithHole = new G4SubtractionSolid(
            "Plunger_Chamber_SphereWithHole", sphereChamber, pipeHole, 0, G4ThreeVector(0, 0, 0));

        // 파이프(실린더) 생성
        G4double pipeInnerR = m_ChamberPipeR - m_ChamberThickness;
        G4double pipeOuterR = m_ChamberPipeR;
        G4double pipeHalfLength = (m_ChamberPipeZ1 - m_ChamberPipeZ0) * 0.5;
        G4double pipePosZ = (m_ChamberPipeZ1 + m_ChamberPipeZ0) * 0.5;

        G4Tubs* pipeForSubtraction = new G4Tubs("Plunger_Chamber_PipeForSubtraction", pipeInnerR, pipeOuterR,
                                                pipeHalfLength, 0. * deg, 360. * deg);

        G4Sphere* sphereForSubtraction = new G4Sphere("Plunger_Chamber_SphereForSubtraction", 0, m_ChamberR + 0.1 * mm,
                                                      0. * deg, 360. * deg, 0. * deg, 180. * deg);  // 0.1mm 오버랩 방지

        G4SubtractionSolid* pipeChamber = new G4SubtractionSolid(
            "Plunger_Chamber_PipeChamber", pipeForSubtraction, sphereForSubtraction, 0,
            G4ThreeVector(0, 0, -pipePosZ));  // 파이프 중심을 구 껍데기 중심에 맞춤 (구는 고정, 파이프만 움직임)

        // 파이프를 구 껍데기와 합치기 (Boolean 연산)
        G4UnionSolid* chamberWithPipe = new G4UnionSolid("Plunger_Chamber_WithPipe", sphereChamberWithHole, pipeChamber,
                                                         0, G4ThreeVector(0, 0, pipePosZ));

        // 논리 볼륨 생성
        G4Material* chamberMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(m_ChamberMaterial);
        m_ChamberLV = new G4LogicalVolume(chamberWithPipe, chamberMaterial, "logic_Plunger_Chamber", 0, 0, 0);

        // 시각화 속성 지정
        m_ChamberLV->SetVisAttributes(m_VisChamber);
    }
    return m_ChamberLV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Plunger::ReadConfiguration(NPL::InputParser parser) {
    vector<NPL::InputBlock*> blocks;
    ////////////////////////////////////////////////////
    // Plunger Target
    ////////////////////////////////////////////////////
    blocks = parser.GetAllBlocksWithTokenAndValue("Plunger", "Target");
    if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl;

    if (blocks.size() > 0) {
        if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << endl << "////  Plunger Target " << endl;
        m_TargetFound = true;
        m_TargetR = blocks[0]->GetDouble("R", "mm");
        m_TargetThickness = blocks[0]->GetDouble("Thickness", "mm");
        m_TargetPosZ = blocks[0]->GetDouble("Z", "mm");
        m_TargetMaterial = blocks[0]->GetString("Material");
    }
    ////////////////////////////////////////////////////
    // Plunger Stopper
    ////////////////////////////////////////////////////
    blocks = parser.GetAllBlocksWithTokenAndValue("Plunger", "Stopper");
    if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl;

    if (blocks.size() > 0) {
        if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << endl << "////  Plunger Stopper " << endl;
        m_StopperFound = true;
        m_StopperR = blocks[0]->GetDouble("R", "mm");
        m_StopperThickness = blocks[0]->GetDouble("Thickness", "mm");
        m_StopperPosZ = blocks[0]->GetDouble("Z", "mm");
        m_StopperMaterial = blocks[0]->GetString("Material");
    }
    ////////////////////////////////////////////////////
    // Plunger Chamber
    ////////////////////////////////////////////////////
    blocks = parser.GetAllBlocksWithTokenAndValue("Plunger", "Chamber");
    if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl;

    if (blocks.size() > 0) {
        if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << endl << "////  Plunger Chamber " << endl;
        m_ChamberFound = true;
        m_ChamberR = blocks[0]->GetDouble("R", "mm");
        m_ChamberThickness = blocks[0]->GetDouble("Thickness", "mm");
        m_ChamberMaterial = blocks[0]->GetString("Material");
        m_ChamberPipeR = blocks[0]->GetDouble("PipeR", "mm");
        m_ChamberPipeZ0 = blocks[0]->GetDouble("PipeZ0", "mm");
        m_ChamberPipeZ1 = blocks[0]->GetDouble("PipeZ1", "mm");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Plunger::SetReactionRegion(G4LogicalVolume* world) {
    ////////////////////////////////////////////////////
    // Reaction Region (only in the target)
    ////////////////////////////////////////////////////
    G4LogicalVolume* reactionRegionLV = m_TargetLV;
    if (!m_ReactionRegion) {
        m_ReactionRegion = new G4Region("NPSimulationReactionProcess");
        m_ReactionRegion->AddRootLogicalVolume(reactionRegionLV);
        m_ReactionRegion->SetUserLimits(new G4UserLimits(1. * mm));
    }
    new NPS::BeamReaction("BeamReaction", m_ReactionRegion);
    // new NPS::Decay("Decay", m_ReactionRegion);

    ////////////////////////////////////////////////////
    // Decay Region (in the whole chamber)
    ////////////////////////////////////////////////////
    // G4Sphere* decayRegionSphere =
    //     new G4Sphere("Plunger_DecayRegion_Sphere", 0, m_ChamberR, 0. * deg, 360. * deg, 0. * deg, 180. * deg);
    // G4Material* decayRegionMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    // G4LogicalVolume* decayRegionLV =
    //     new G4LogicalVolume(decayRegionSphere, decayRegionMaterial, "logic_Plunger_DecayRegion", 0, 0, 0);
    // decayRegionLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    // new G4PVPlacement(G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, 0)), decayRegionLV, "PlungerDecayRegion",
    //                   world, false, 0);

    // if (!m_DecayRegion) {
    //     m_DecayRegion = new G4Region("NPSimulationDecayProcess");
    //     m_DecayRegion->AddRootLogicalVolume(decayRegionLV);
    //     m_DecayRegion->SetUserLimits(new G4UserLimits(1. * mm));
    // }
    // new NPS::Decay("Decay", m_DecayRegion);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Plunger::ConstructDetector(G4LogicalVolume* world) {
    if (m_TargetFound) {
        new G4PVPlacement(G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, m_TargetPosZ)), BuildTarget(),
                          "PlungerTarget", world, false, 0);
    }
    if (m_StopperFound) {
        new G4PVPlacement(G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, m_StopperPosZ)), BuildStopper(),
                          "PlungerStopper", world, false, 0);
    }
    if (m_ChamberFound) {
        new G4PVPlacement(G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, 0)), BuildChamber(), "PlungerChamber",
                          world, false, 0);
    }
    SetReactionRegion(world);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Plunger::InitializeRootOutput() {
    RootOutput* pAnalysis = RootOutput::getInstance();
    TTree* pTree = pAnalysis->GetTree();
    if (!pTree->FindBranch("Plunger")) {
        pTree->Branch("Plunger", "TPlungerData", &m_Event);
    }
    pTree->SetBranchAddress("Plunger", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Plunger::ReadSensitive(const G4Event*) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void Plunger::InitializeScorers() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Plunger::Construct() { return (NPS::VDetector*)new Plunger(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_nps_Plunger {
   public:
    proxy_nps_Plunger() {
        NPS::DetectorFactory::getInstance()->AddToken("Plunger", "Plunger");
        NPS::DetectorFactory::getInstance()->AddDetector("Plunger", Plunger::Construct);
    }
};

proxy_nps_Plunger p_nps_Plunger;
}
