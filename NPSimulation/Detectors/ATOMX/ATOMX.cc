/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: jungwoo  contact address: phyjics@gmail.com                        *
 *                                                                           *
 * Creation Date  : May 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  ATOMX simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#ifndef jw_cout
#include <string>
#define jw_cout std::cout<<"\033[0;32m"<<Form("+%d %s # \033[0m",__LINE__,std::string(__FILE__).c_str())
#endif

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ProductionCuts.hh"

// NPTool header
#include "ATOMX.hh"
#include "ATOMXScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

#include "STARK.hh"

using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ATOMX Specific Method
ATOMX::ATOMX(){
    m_Event = new TATOMXData();
    m_ATOMXScorer = 0;
    m_StepLimitRegion = NULL;
    m_Pos = G4ThreeVector(0,0,0);
    m_ChamberSize = G4ThreeVector(480, 450, 480);
    m_GasSize = G4ThreeVector(450, 400, 450);
    m_PadSize = G4ThreeVector(400, 2  , 400);
    m_MMSSize = G4ThreeVector(400, 20 , 400);
    m_ReactionBoxSize = 20.;
    m_Mylar_Rmax = 3.5;
    m_Mylar_Thickness = 7;
    m_GasMaterial.clear();
    m_GasMaterial.push_back("He");
    m_GasMaterial.push_back("CO2");
    m_GasFraction.clear();
    m_GasFraction.push_back(97);
    m_GasFraction.push_back(3);
    m_Temperature = 295;
    m_Pressure = 0.1;
    m_ReactionZ1 = 0;
    m_ReactionZ2 = 0;

}

ATOMX::~ATOMX(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* ATOMX::BuildDetector()
{
    bool limitReactionZ = (m_ReactionZWidth>0);
    if (limitReactionZ)
        cout << "Limiting sensitive region : " << m_ReactionZ1 << " -> " << m_ReactionZ2 << endl;

    G4Material* Cu = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
    G4Material* Al = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4Material* Mylar = MaterialManager::getInstance()->GetMaterialFromLibrary("Mylar");

    double chamberThicknessZ = 0.5 * (m_ChamberSize.z() - m_GasSize.z());
    double throughZ = -0.5*m_ChamberSize.z() + 0.5*chamberThicknessZ;
    double windowZ  = -0.5*chamberThicknessZ + 0.5*m_Mylar_Thickness;
    double sensitiveZ = -0.5*m_GasSize.z() + 0.5*(m_ReactionZ1+m_ReactionZ2);
    double padY = 0.5*m_GasSize.y();

    G4Box*  solidChamber  = new G4Box ("solid_Box", m_ChamberSize.x() * 0.5, m_ChamberSize.y() * 0.5, m_ChamberSize.z() * 0.5);
    G4Box*  solidGas      = new G4Box ("solid_Gas", m_GasSize.x() * 0.5, m_GasSize.y() * 0.5, m_GasSize.z() * 0.5);
    G4Box*  solidPad      = new G4Box ("solid_Pad", m_PadSize.x() * 0.5, m_PadSize.y() * 0.5, m_PadSize.z() * 0.5);
    G4Box*  solidMMS      = new G4Box ("solid_MMS", m_MMSSize.x() * 0.5, m_MMSSize.y() * 0.5, m_MMSSize.z() * 0.5);
    G4Tubs* solidThrough1 = new G4Tubs("solid_Th1", 0, m_Mylar_Rmax, chamberThicknessZ * 0.5, 0 * deg, 360 * deg);
    G4Tubs* solidWindow   = new G4Tubs("solid_Win", 0, m_Mylar_Rmax, m_Mylar_Thickness * 0.5, 0 * deg, 360 * deg);
    G4Box*  solidSensitive = nullptr;
    if (limitReactionZ)
        solidSensitive = new G4Box("solid_Sen", m_ReactionBoxSize * 0.5, m_ReactionBoxSize * 0.5, m_ReactionZWidth * 0.5);

    unsigned const int NumberOfGasMix = m_GasMaterial.size();

    double density = 0;
    double density_sum = 0;
    vector<G4Material*> GasComponent;
    vector<double> FractionMass;

    for (unsigned int i = 0; i < NumberOfGasMix; i++) {
        GasComponent.push_back(
                MaterialManager::getInstance()->GetGasFromLibrary(m_GasMaterial[i], m_Pressure, m_Temperature));
    }
    for (unsigned int i = 0; i < NumberOfGasMix; i++) {
        density += ((double)m_GasFraction[i] / 100) * GasComponent[i]->GetDensity();
        density_sum += GasComponent[i]->GetDensity();
    }

    for (unsigned int i = 0; i < NumberOfGasMix; i++) {
        FractionMass.push_back(GasComponent[i]->GetDensity() / density_sum);
    }

    G4Material* GasMaterial = new G4Material("GasMix", density, NumberOfGasMix, kStateGas, m_Temperature, m_Pressure);
    G4Material* DriftGasMaterial =
        new G4Material("DriftGasMix", density, NumberOfGasMix, kStateGas, m_Temperature, m_Pressure);

    for (unsigned int i = 0; i < NumberOfGasMix; i++) {
        GasMaterial->AddMaterial(GasComponent[i], FractionMass[i]);
        DriftGasMaterial->AddMaterial(GasComponent[i], FractionMass[i]);
    }
    for (unsigned int i = 0; i < NumberOfGasMix; i++)
        cout << i << " " << GasComponent[i] << " " << FractionMass[i] << endl;

    m_logicChamber = new G4LogicalVolume(solidChamber, GasMaterial, "logic_Chamber", 0, 0, 0);
    m_logicGas = new G4LogicalVolume(solidGas, DriftGasMaterial, "ATOMX_LV_Gas", 0, 0, 0);
    G4LogicalVolume* logicThrough1 = new G4LogicalVolume(solidThrough1, DriftGasMaterial, "logic_Gas", 0, 0, 0);
    if (limitReactionZ)
        m_sensitive = new G4LogicalVolume(solidSensitive, DriftGasMaterial, "logic_SensitiveGas", 0, 0, 0);
    G4LogicalVolume* logicPad = new G4LogicalVolume(solidPad, Cu, "logic_Pad", 0, 0, 0);
    G4LogicalVolume* logicMMS = new G4LogicalVolume(solidMMS, Al, "logic_MMS", 0, 0, 0);
    G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow, Mylar, "logic_Win", 0, 0, 0);

    int copyNo = fCopyNo;
    G4RotationMatrix* Rot = new G4RotationMatrix();
    new G4PVPlacement(G4Transform3D(*Rot, G4ThreeVector(0, 0, 0)), m_logicGas, "ATOMX_Gas", m_logicChamber, false, ++copyNo);
    new G4PVPlacement(G4Transform3D(*Rot, G4ThreeVector(0, 0, throughZ)), logicThrough1, "ATOMX_Through1", m_logicChamber, false, ++copyNo);
    if (limitReactionZ)
    new G4PVPlacement(G4Transform3D(*Rot, G4ThreeVector(0, 0, sensitiveZ)), m_sensitive, "ATOMX_SensitiveGas", m_logicGas, false, ++copyNo);
    new G4PVPlacement(G4Transform3D(*Rot, G4ThreeVector(0, padY, 0)), logicPad, "ATOMX_Pad", m_logicGas, false, ++copyNo);
    new G4PVPlacement(G4Transform3D(*Rot, G4ThreeVector(0, 0, windowZ)), logicWindow, "ATOMX_Window", logicThrough1, false, ++copyNo);

    m_logicGas -> SetSensitiveDetector(m_ATOMXScorer);
    if (limitReactionZ)
        m_sensitive -> SetSensitiveDetector(m_ATOMXScorer);

    G4VisAttributes* visChamber = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.3));
    G4VisAttributes* visGas = new G4VisAttributes(G4Colour(0, 0.5, 0.5, 0.3));
    G4VisAttributes* visWindow = new G4VisAttributes(G4Colour(1, 0, 0, 0.25));
    G4VisAttributes* visPad = new G4VisAttributes(G4Colour(255, 223, 50, 0.8));
    G4VisAttributes* visMMS = new G4VisAttributes(G4Colour(100, 100, 100, 0.4));
    G4VisAttributes* visSensitive = new G4VisAttributes(G4Colour(0, 0.3, 0.3, 0.3));
    visPad -> SetForceWireframe(true);

    m_logicChamber -> SetVisAttributes(visChamber);
    m_logicGas     -> SetVisAttributes(visGas);
    logicThrough1  -> SetVisAttributes(visGas);
    logicWindow    -> SetVisAttributes(visWindow);
    logicPad       -> SetVisAttributes(visPad);
    logicMMS       -> SetVisAttributes(visMMS);
    if (limitReactionZ)
        m_sensitive -> SetVisAttributes(visSensitive);

    return m_logicChamber;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void ATOMX::ReadConfiguration(NPL::InputParser parser){

    vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ATOMX");
    if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl; 

    vector<string> tokens = { "TPC_Pos" };

    for(unsigned int i = 0 ; i < blocks.size() ; i++){
        if(blocks[i]->HasTokenList(tokens))
        {
            if(NPOptionManager::getInstance()->GetVerboseLevel())
                cout << endl << "////  ATOMX " << i+1 <<  endl;
            vector<string> GasName;
            vector<int> GasFraction;
            vector<int> vReactionZ;
            if (blocks[i]->HasToken("TPC_Pos"             )) m_Pos = NPS::ConvertVector(blocks[i]->GetTVector3("TPC_Pos","mm"));
            if (blocks[i]->HasToken("TPC_Chamber_Size"    )) m_ChamberSize = NPS::ConvertVector(blocks[i]->GetTVector3("TPC_Chamber_Size","mm"));
            if (blocks[i]->HasToken("TPC_Gas_Size"        )) m_GasSize = NPS::ConvertVector(blocks[i]->GetTVector3("TPC_Gas_Size","mm"));
            if (blocks[i]->HasToken("TPC_Pad_Size"        )) m_PadSize = NPS::ConvertVector(blocks[i]->GetTVector3("TPC_Pad_Size","mm"));
            if (blocks[i]->HasToken("TPC_MMS_Size"        )) m_MMSSize = NPS::ConvertVector(blocks[i]->GetTVector3("TPC_MMS_Size","mm"));
            if (blocks[i]->HasToken("TPC_Mylar_Rmax"      )) m_Mylar_Rmax = blocks[i]->GetDouble("TPC_Mylar_Rmax","cm");
            if (blocks[i]->HasToken("TPC_Mylar_Thickness" )) m_Mylar_Thickness = blocks[i]->GetDouble("TPC_Mylar_Thickness","micrometer");
            if (blocks[i]->HasToken("Gas_EProductionCut"  )) m_EProductionCut = blocks[i]->GetDouble("Gas_EProductionCut", "mm");
            if (blocks[i]->HasToken("Gas_StepLimit"       )) m_StepLimit = blocks[i]->GetDouble("Gas_StepLimit", "mm");
            if (blocks[i]->HasToken("Gas_Material"        )) GasName = blocks[i]->GetVectorString("Gas_Material");
            if (blocks[i]->HasToken("Gas_Fraction"        )) GasFraction = blocks[i]->GetVectorInt("Gas_Fraction");
            if (blocks[i]->HasToken("Gas_Temperature"     )) m_Temperature = blocks[i]->GetDouble("Gas_Temperature", "kelvin");
            if (blocks[i]->HasToken("Gas_Pressure"        )) m_Pressure = blocks[i]->GetDouble("Gas_Pressure", "bar");
            if (blocks[i]->HasToken("Gas_ReactionBox_Size")) m_ReactionBoxSize = blocks[i]->GetDouble("Gas_ReactionBox_Size","mm");
            if (GasFraction.size()==GasName.size()&&GasName.size()>0) {
                m_GasMaterial.clear();
                m_GasFraction.clear();
                for (unsigned int j = 0; j < GasName.size(); j++) {
                    m_GasMaterial.push_back(GasName[j]);
                    m_GasFraction.push_back(GasFraction[j]);
                }
            }
            vReactionZ = blocks[i]->GetVectorInt("Gas_ReactionZ");
            if (vReactionZ.size()==2) {
                m_ReactionZ1 = vReactionZ[0];
                m_ReactionZ2 = vReactionZ[1];
            }
            if (m_ReactionZ2==m_ReactionZ1) continue;
            else {
                if (m_ReactionZ2<m_ReactionZ1) {
                    m_ReactionZ1 = vReactionZ[1];
                    m_ReactionZ2 = vReactionZ[0];
                }
                m_ReactionZWidth = m_ReactionZ2 - m_ReactionZ1;
            }
        }
        else{
            cout << "ERROR: check your input file formatting (ATOMX::ReadConfiguration)" << endl;
            exit(1);
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void ATOMX::ConstructDetector(G4LogicalVolume* world)
{
    G4ThreeVector Det_pos = m_Pos;
    G4RotationMatrix* Rot = new G4RotationMatrix();
    BuildDetector();
    new G4PVPlacement(G4Transform3D(*Rot,Det_pos), m_logicChamber, "ATOMX", world, false, fCopyNo);

    G4ProductionCuts* ecut = new G4ProductionCuts();
    ecut->SetProductionCut(m_EProductionCut * mm, "e-");

    G4Region *region_reaction = new G4Region("NPSimulationProcess");
    region_reaction->SetProductionCuts(ecut);
    region_reaction->SetUserLimits(new G4UserLimits(m_StepLimit * mm));

    G4Region* region_others = new G4Region("RegionCut");
    region_others->SetProductionCuts(ecut);

    if (m_ReactionZWidth>0) {
        cout << "Reaction sensitive" << endl;
        region_reaction->AddRootLogicalVolume(m_sensitive);
        G4UserLimits* limits = new G4UserLimits(m_StepLimit * mm);
        m_logicGas -> SetUserLimits(limits);
    }
    else
        region_reaction->AddRootLogicalVolume(m_logicGas);

    //region_others->AddRootLogicalVolume(m_logicGas);
    region_others->AddRootLogicalVolume(m_logicChamber);

    G4FastSimulationManager* mng = region_reaction->GetFastSimulationManager();
    unsigned int size = m_ReactionModel.size();
    for (unsigned int i = 0; i < size; i++) {
        mng->RemoveFastSimulationModel(m_ReactionModel[i]);
    }
    m_ReactionModel.clear();
    G4VFastSimulationModel* fsm;
    fsm = new NPS::BeamReaction("BeamReaction", region_reaction);
    m_ReactionModel.push_back(fsm);
    fsm = new NPS::Decay("Decay", region_reaction);
    m_ReactionModel.push_back(fsm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void ATOMX::InitializeRootOutput()
{
    RootOutput *pAnalysis = RootOutput::getInstance();
    TTree *pTree = pAnalysis->GetTree();
    if(!pTree->FindBranch("ATOMX")){
        pTree->Branch("ATOMX", "TATOMXData", &m_Event) ;
    }
    pTree->SetBranchAddress("ATOMX", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void ATOMX::ReadSensitive(const G4Event*)
{
    m_Event->Clear();

    ATOMXScorers::PS_ATOMX* interaction = (ATOMXScorers::PS_ATOMX*) m_ATOMXScorer->GetPrimitive(0);

    unsigned int size = interaction -> GetMult(); 
    for (unsigned int i=0; i<size; i++)
    {
        TVector3 Position(interaction->GetPositionX(i), interaction->GetPositionY(i), interaction->GetPositionZ(i));
        m_Event->fEnergyLoss.push_back(interaction -> GetEnergy(i));
        m_Event->fTime.push_back(interaction -> GetTime(i));
        m_Event->fPosition.push_back(Position);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void ATOMX::InitializeScorers()
{
    bool already_exist = false; 
    m_ATOMXScorer = CheckScorer("ATOMXScorer",already_exist) ;

    if(already_exist) 
        return ;

    G4VPrimitiveScorer* Interaction= new ATOMXScorers::PS_ATOMX("ATOMX", 0) ;
    m_ATOMXScorer->RegisterPrimitive(Interaction);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_ATOMXScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* ATOMX::Construct(){
    return  (NPS::VDetector*) new ATOMX();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
    class proxy_nps_ATOMX{
        public:
            proxy_nps_ATOMX(){
                NPS::DetectorFactory::getInstance()->AddToken("ATOMX","ATOMX");
                NPS::DetectorFactory::getInstance()->AddDetector("ATOMX",ATOMX::Construct);
            }
    };

    proxy_nps_ATOMX p_nps_ATOMX;
}
