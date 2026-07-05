#include "CSSU.hh"
#include <cstdlib>

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <tuple>

#include "CalorimeterScorers.hh"
#include "CLHEP/Random/RandGauss.h"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "InteractionScorers.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSFunction.hh"
#include "NPSDetectorFactory.hh"
#include "NPSHitsMap.hh"
#include "PhotoDiodeScorers.hh"
#include "RootOutput.h"

using namespace std;
using namespace CLHEP;

namespace CSSU_NS {
const double ResoTime = 1.0 * ns;
const double ResoEnergy = 1.0 * keV;
const double EnergyThreshold = 1.0 * keV;
const double DefaultScintYield = 10000. / MeV;
const double DefaultAttenuationLength = 2.0 * m;
const double DefaultScintillatorRIndex = 1.58;
}  // namespace CSSU_NS

CSSU::CSSU() {
  m_Event = new TCSSUData();
  m_CSSUScorer = 0;
  m_PMTScorer = 0;

  m_VisMother = new G4VisAttributes(G4Colour(1, 1, 1, 0.05));
  m_VisMother->SetVisibility(false);
  m_VisScintillator = new G4VisAttributes(G4Colour(0.0, 0.25, 1.0, 0.35));
  m_VisPMT = new G4VisAttributes(G4Colour(0.1, 0.1, 0.1, 0.8));
  m_VisGasMother = new G4VisAttributes(G4Colour(0.1, 0.5, 1.0, 0.12));
  m_VisGasMother->SetForceWireframe(true);
  m_VisGas = new G4VisAttributes(G4Colour(0.2, 0.8, 1.0, 0.45));
  m_VisGas->SetForceSolid(true);
}

CSSU::~CSSU() {}

void CSSU::AddScintillatorBar(G4ThreeVector pos, double length, double width, double thickness, double pmtLength,
                              double pmtDiameter, double scintYield, double attenuationLength, double scintillatorRIndex) {
  m_Bars.push_back({pos, length, width, thickness, pmtLength, pmtDiameter, scintYield, attenuationLength,
                    scintillatorRIndex});
}

void CSSU::AddGasBox(G4ThreeVector pos, double sizeX, double sizeY, double sizeZ, double pressureTorr,
                     double stepLimit, int segmentsZ, double segmentGapZ, double energyThreshold, double resoEnergy) {
  m_GasBoxes.push_back(
      {pos, sizeX, sizeY, sizeZ, pressureTorr, stepLimit, std::max(1, segmentsZ), std::max(0.0, segmentGapZ),
       energyThreshold, resoEnergy});
}

void CSSU::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("CSSU");
  if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << "//// " << blocks.size() << " CSSU blocks found" << endl;

  vector<string> posToken = {"POS", "Type", "Length", "Width", "Thickness"};
  vector<string> xyzToken = {"X", "Y", "Z", "Type", "Length", "Width", "Thickness"};

  for (unsigned int i = 0; i < blocks.size(); ++i) {
    G4ThreeVector pos;
    if (blocks[i]->HasTokenList(posToken)) {
      pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS", "mm"));
    }
    else if (blocks[i]->HasTokenList(xyzToken)) {
      pos = G4ThreeVector(blocks[i]->GetDouble("X", "mm"), blocks[i]->GetDouble("Y", "mm"),
                          blocks[i]->GetDouble("Z", "mm"));
    }
    else {
      cout << "ERROR: CSSU block requires POS or X/Y/Z plus Type, Length, Width and Thickness" << endl;
      exit(1);
    }

    string type = blocks[i]->GetString("Type");
    if (type == "ScintillatorBar") {
      double length = blocks[i]->GetDouble("Length", "mm");
      double width = blocks[i]->GetDouble("Width", "mm");
      double thickness = blocks[i]->GetDouble("Thickness", "mm");
      double pmtLength = 1.0 * mm;
      if (blocks[i]->HasToken("PMTLength")) pmtLength = blocks[i]->GetDouble("PMTLength", "mm");
      else if (blocks[i]->HasToken("PMTThickness")) pmtLength = blocks[i]->GetDouble("PMTThickness", "mm");

      double pmtDiameter = std::max(width, thickness);
      if (blocks[i]->HasToken("PMTDiameter")) pmtDiameter = blocks[i]->GetDouble("PMTDiameter", "mm");
      else if (blocks[i]->HasToken("PMTFace")) pmtDiameter = blocks[i]->GetDouble("PMTFace", "mm");
      double scintYield = CSSU_NS::DefaultScintYield;
      if (blocks[i]->HasToken("ScintillationYield")) {
        scintYield = std::stod(blocks[i]->GetValue("ScintillationYield")) / MeV;
      }
      double attenuationLength = blocks[i]->HasToken("AttenuationLength")
                                     ? blocks[i]->GetDouble("AttenuationLength", "mm")
                                     : CSSU_NS::DefaultAttenuationLength;
      double scintillatorRIndex = blocks[i]->HasToken("ScintillatorRIndex")
                                      ? std::stod(blocks[i]->GetValue("ScintillatorRIndex"))
                                      : CSSU_NS::DefaultScintillatorRIndex;

      AddScintillatorBar(pos, length, width, thickness, pmtLength, pmtDiameter, scintYield, attenuationLength,
                         scintillatorRIndex);
    }
    else if (type == "GasBox") {
      double sizeX = blocks[i]->HasToken("SizeX") ? blocks[i]->GetDouble("SizeX", "mm") : blocks[i]->GetDouble("Length", "mm");
      double sizeY = blocks[i]->HasToken("SizeY") ? blocks[i]->GetDouble("SizeY", "mm") : blocks[i]->GetDouble("Width", "mm");
      double sizeZ = blocks[i]->HasToken("SizeZ") ? blocks[i]->GetDouble("SizeZ", "mm") : blocks[i]->GetDouble("Thickness", "mm");
      double pressureTorr = blocks[i]->HasToken("Pressure") ? std::stod(blocks[i]->GetValue("Pressure")) : 200.0;
      if (blocks[i]->HasToken("PressureTorr")) pressureTorr = std::stod(blocks[i]->GetValue("PressureTorr"));
      double stepLimit = blocks[i]->HasToken("StepLimit") ? blocks[i]->GetDouble("StepLimit", "mm") : 10.0 * mm;
      int segmentsZ = 1;
      if (blocks[i]->HasToken("SegmentsZ")) segmentsZ = std::stoi(blocks[i]->GetValue("SegmentsZ"));
      else if (blocks[i]->HasToken("ActiveSegmentsZ")) segmentsZ = std::stoi(blocks[i]->GetValue("ActiveSegmentsZ"));
      else if (blocks[i]->HasToken("NbSlices")) segmentsZ = std::stoi(blocks[i]->GetValue("NbSlices"));
      double segmentGapZ = 0;
      if (blocks[i]->HasToken("SegmentGapZ")) segmentGapZ = blocks[i]->GetDouble("SegmentGapZ", "mm");
      else if (blocks[i]->HasToken("GapZ")) segmentGapZ = blocks[i]->GetDouble("GapZ", "mm");
      double energyThreshold = blocks[i]->HasToken("EnergyThreshold")
                                   ? blocks[i]->GetDouble("EnergyThreshold", "MeV")
                                   : CSSU_NS::EnergyThreshold;
      double resoEnergy = blocks[i]->HasToken("ResoEnergy")
                              ? blocks[i]->GetDouble("ResoEnergy", "keV")
                              : CSSU_NS::ResoEnergy;
      AddGasBox(pos, sizeX, sizeY, sizeZ, pressureTorr, stepLimit, segmentsZ, segmentGapZ, energyThreshold,
                resoEnergy);
    }
    else {
      cout << "ERROR: CSSU Type " << type << " is not implemented yet" << endl;
      exit(1);
    }
  }
}

void CSSU::ConstructDetector(G4LogicalVolume* world) {
  for (unsigned int i = 0; i < m_Bars.size(); ++i) {
    const int det = i + 1;
    const ScintillatorBar& bar = m_Bars[i];
    G4LogicalVolume* mother = BuildMotherVolume(bar, det);
    G4LogicalVolume* scint = BuildScintillatorBar(bar, det);
    G4LogicalVolume* pmt = BuildPMT(bar, det);

    new G4PVPlacement(0, bar.pos, mother, "Mother", world, false, det);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), scint, "ScintillatorBar", mother, false, det);

    G4RotationMatrix* pmtRotation = new G4RotationMatrix();
    pmtRotation->rotateY(90.0 * deg);

    const double pmtX = 0.5 * bar.length + 0.5 * bar.pmtLength;
    new G4PVPlacement(pmtRotation, G4ThreeVector(-pmtX, 0, 0), pmt, "PMT_Left", mother, false, 2 * det - 1);
    new G4PVPlacement(pmtRotation, G4ThreeVector(+pmtX, 0, 0), pmt, "PMT_Right", mother, false, 2 * det);
  }

  int gasOutputDetNumber = m_Bars.size() + 1;
  int gasSliceCopyNumber = m_Bars.size() + m_GasBoxes.size() + 1;
  for (unsigned int i = 0; i < m_GasBoxes.size(); ++i) {
    const GasBox& gas = m_GasBoxes[i];
    const int outputDet = gasOutputDetNumber++;
    m_GasResponseByOutputDetector[outputDet] = {gas.energyThreshold, gas.resoEnergy};

    const int segmentsZ = std::max(1, gas.segmentsZ);
    const double totalGapZ = gas.segmentGapZ * (segmentsZ - 1);
    if (totalGapZ >= gas.sizeZ) {
      cout << "ERROR: CSSU GasBox SegmentGapZ leaves no active gas thickness" << endl;
      exit(1);
    }
    const double sliceZ = (gas.sizeZ - totalGapZ) / segmentsZ;
    G4Material* gasMaterial = BuildCF4Gas(gas, outputDet);
    G4LogicalVolume* gasMother = BuildGasMother(gas, outputDet, gasMaterial);
    new G4PVPlacement(0, gas.pos, gasMother, "CSSU_GasBox", world, false, outputDet);

    for (int s = 0; s < segmentsZ; ++s) {
      const int det = gasSliceCopyNumber++;
      m_GasSliceToOutputDetector[det] = outputDet;
      const double localZ = -0.5 * gas.sizeZ + 0.5 * sliceZ + s * (sliceZ + gas.segmentGapZ);
      G4LogicalVolume* gasLogic = BuildGasSlice(gas, det, sliceZ, gasMaterial);
      new G4PVPlacement(0, G4ThreeVector(0, 0, localZ), gasLogic, "CSSU_GasSlice", gasMother, false, det);
    }
  }
}

G4Material* CSSU::BuildOpticalAir() {
  static G4Material* material = 0;
  if (material) return material;

  material = new G4Material("OpticalAir", 1.290 * mg / cm3, 2);
  material->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("N"), 7);
  material->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("O"), 3);

  G4double energy[] = {1.5 * eV, 4.0 * eV};
  G4double rindex[] = {1.0, 1.0};
  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", energy, rindex, 2);
  material->SetMaterialPropertiesTable(mpt);
  return material;
}

G4Material* CSSU::BuildScintillatorMaterial(const ScintillatorBar& bar, int detectorNumber) {
  ostringstream name;
  name << "BC400_Scintillator_" << detectorNumber;
  G4Material* material = new G4Material(name.str(), 1.032 * g / cm3, 2);
  material->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("H"), 10);
  material->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("C"), 9);

  G4double energy[] = {2.0 * eV, 2.5 * eV, 3.0 * eV};
  G4double scint[] = {0.4, 1.0, 0.4};
  G4double rindex[] = {bar.scintillatorRIndex, bar.scintillatorRIndex, bar.scintillatorRIndex};
  G4double absorption[] = {bar.attenuationLength, bar.attenuationLength, bar.attenuationLength};

  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("SCINTILLATIONCOMPONENT1", energy, scint, 3);
  mpt->AddProperty("RINDEX", energy, rindex, 3);
  mpt->AddProperty("ABSLENGTH", energy, absorption, 3);
  mpt->AddConstProperty("SCINTILLATIONYIELD", bar.scintYield);
  mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
  mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.4 * ns);
  mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  material->SetMaterialPropertiesTable(mpt);
  return material;
}

// n=1 vacuum: optical photons propagate through it (RINDEX defined) but, being
// vacuum, it does NOT slow the incoming alpha — so the mother envelope can be large
// (photons fly a visible distance before dying at the world boundary) without
// stopping the beam in air. Tunable via NPS_CSSU_ENVELOPE_MM (default 50 mm).
G4Material* CSSU::BuildOpticalVacuum() {
  static G4Material* material = 0;
  if (material) return material;
  G4Material* vac = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  material = new G4Material("CSSU_OpticalVacuum", vac->GetDensity(), 1, vac->GetState(),
                            vac->GetTemperature(), vac->GetPressure());
  material->AddMaterial(vac, 1.0);
  G4double energy[] = {1.5 * eV, 4.0 * eV};
  G4double rindex[] = {1.0, 1.0};
  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", energy, rindex, 2);
  material->SetMaterialPropertiesTable(mpt);
  return material;
}

G4LogicalVolume* CSSU::BuildMotherVolume(const ScintillatorBar& bar, int detectorNumber) {
  ostringstream name;
  name << "Mother_" << detectorNumber;
  double margin = 50.0 * mm;   // vacuum envelope — safe to be large (see BuildOpticalVacuum)
  if (const char* e = std::getenv("NPS_CSSU_ENVELOPE_MM")) {
    double v = atof(e);
    if (v > 0) margin = v * mm;
  }
  G4Box* solid = new G4Box(name.str(), 0.5 * (bar.length + 2.0 * bar.pmtLength + 2.0 * margin),
                           0.5 * (std::max(bar.width, bar.pmtDiameter) + 2.0 * margin),
                           0.5 * (std::max(bar.thickness, bar.pmtDiameter) + 2.0 * margin));
  G4LogicalVolume* logic = new G4LogicalVolume(solid, BuildOpticalVacuum(), name.str(), 0, 0, 0);
  logic->SetVisAttributes(m_VisMother);
  return logic;
}

G4LogicalVolume* CSSU::BuildScintillatorBar(const ScintillatorBar& bar, int detectorNumber) {
  ostringstream name;
  name << "ScintillatorBar_" << detectorNumber;
  G4Box* solid = new G4Box(name.str(), 0.5 * bar.length, 0.5 * bar.width, 0.5 * bar.thickness);
  G4LogicalVolume* logic = new G4LogicalVolume(solid, BuildScintillatorMaterial(bar, detectorNumber), name.str(), 0, 0, 0);
  logic->SetVisAttributes(m_VisScintillator);
  logic->SetSensitiveDetector(m_CSSUScorer);
  return logic;
}

G4LogicalVolume* CSSU::BuildPMT(const ScintillatorBar& bar, int detectorNumber) {
  ostringstream name;
  name << "Photocathode_" << detectorNumber;
  G4Tubs* solid = new G4Tubs(name.str(), 0, 0.5 * bar.pmtDiameter, 0.5 * bar.pmtLength, 0, 360 * deg);
  G4Material* material = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  G4LogicalVolume* logic = new G4LogicalVolume(solid, material, name.str(), 0, 0, 0);
  logic->SetVisAttributes(m_VisPMT);
  logic->SetSensitiveDetector(m_PMTScorer);
  return logic;
}

G4Material* CSSU::BuildCF4Gas(const GasBox& gas, int detectorNumber) {
  ostringstream name;
  name << "CSSU_CF4_" << detectorNumber << "_" << gas.pressureTorr << "Torr";

  const double densityAt200Torr = 0.947 * mg / cm3;
  const double pressureAt200Torr = 26664.5 * hep_pascal;
  const double scale = gas.pressureTorr / 200.0;
  G4Material* material =
      new G4Material(name.str(), densityAt200Torr * scale, 2, kStateGas, 293.15 * kelvin, pressureAt200Torr * scale);

  G4NistManager* nist = G4NistManager::Instance();
  material->AddElement(nist->FindOrBuildElement("C"), 1);
  material->AddElement(nist->FindOrBuildElement("F"), 4);
  return material;
}

G4LogicalVolume* CSSU::BuildGasMother(const GasBox& gas, int detectorNumber, G4Material* material) {
  ostringstream name;
  name << "CSSU_GasMother_" << detectorNumber;
  G4Box* solid = new G4Box(name.str(), 0.5 * gas.sizeX, 0.5 * gas.sizeY, 0.5 * gas.sizeZ);
  G4LogicalVolume* logic = new G4LogicalVolume(solid, material, name.str(), 0, 0, 0);
  if (gas.stepLimit > 0) logic->SetUserLimits(new G4UserLimits(gas.stepLimit));
  logic->SetVisAttributes(m_VisGasMother);
  return logic;
}

G4LogicalVolume* CSSU::BuildGasSlice(const GasBox& gas, int detectorNumber, double sliceZ, G4Material* material) {
  ostringstream name;
  name << "CSSU_GasSlice_" << detectorNumber;
  G4Box* solid = new G4Box(name.str(), 0.5 * gas.sizeX, 0.5 * gas.sizeY, 0.5 * sliceZ);
  G4LogicalVolume* logic = new G4LogicalVolume(solid, material, name.str(), 0, 0, 0);
  if (gas.stepLimit > 0) logic->SetUserLimits(new G4UserLimits(gas.stepLimit));
  logic->SetVisAttributes(m_VisGas);
  logic->SetSensitiveDetector(m_CSSUScorer);
  return logic;
}

void CSSU::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree* pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("CSSU")) {
    pTree->Branch("CSSU", "TCSSUData", &m_Event);
  }
  pTree->SetBranchAddress("CSSU", &m_Event);
}

void CSSU::ReadSensitive(const G4Event* event) {
  m_Event->Clear();

  CalorimeterScorers::PS_Calorimeter* calo = (CalorimeterScorers::PS_Calorimeter*)m_CSSUScorer->GetPrimitive(0);
  unsigned int caloMult = calo->GetMult();
  map<int, pair<double, double> > gasEnergyTime;
  for (unsigned int i = 0; i < caloMult; ++i) {
    int det = calo->GetLevel(i)[0];
    auto gasSlice = m_GasSliceToOutputDetector.find(det);
    if (gasSlice != m_GasSliceToOutputDetector.end()) {
      const int outputDet = gasSlice->second;
      auto& summed = gasEnergyTime[outputDet];
      summed.first += calo->GetEnergy(i);
      if (summed.second == 0 || calo->GetTime(i) < summed.second) summed.second = calo->GetTime(i);
    }
    else {
      double energy = RandGauss::shoot(calo->GetEnergy(i), CSSU_NS::ResoEnergy);
      if (energy > CSSU_NS::EnergyThreshold) {
        double time = RandGauss::shoot(calo->GetTime(i), CSSU_NS::ResoTime);
        m_Event->SetEnergy(det, energy, time);
      }
    }
  }

  for (const auto& gas : gasEnergyTime) {
    const int det = gas.first;
    const auto response = m_GasResponseByOutputDetector.find(det);
    const double threshold = response == m_GasResponseByOutputDetector.end() ? CSSU_NS::EnergyThreshold
                                                                             : response->second.energyThreshold;
    const double reso = response == m_GasResponseByOutputDetector.end() ? CSSU_NS::ResoEnergy : response->second.resoEnergy;
    double energy = RandGauss::shoot(gas.second.first, reso);
    if (energy > threshold) {
      double time = RandGauss::shoot(gas.second.second, CSSU_NS::ResoTime);
      m_Event->SetEnergy(det, energy, time);
    }
  }

  G4int collectionID = G4SDManager::GetSDMpointer()->GetCollectionID("CSSUPMTScorer/PhotoDiode");
  NPS::HitsMap<G4double*>* hitMap = (NPS::HitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(collectionID));
  map<pair<int, int>, pair<int, double> > pmtData;

  for (const auto& hit : *hitMap->GetMap()) {
    G4double* info = *(hit.second);
    if (info[8] <= 0) continue;

    int copy = (int)info[7];
    int det = (copy + 1) / 2;
    int pmt = (copy % 2) ? 1 : 2;
    pair<int, int> key(det, pmt);

    auto it = pmtData.find(key);
    if (it == pmtData.end()) {
      pmtData[key] = make_pair((int)info[8], info[1]);
    }
    else {
      it->second.first += (int)info[8];
      if (info[1] < it->second.second) it->second.second = info[1];
    }
  }

  for (const auto& pmt : pmtData) {
    m_Event->SetPMT(pmt.first.first, pmt.first.second, pmt.second.first, pmt.second.second);
  }
  hitMap->clear();
}

void CSSU::InitializeScorers() {
  bool cssuAlreadyExist = false;
  m_CSSUScorer = CheckScorer("CSSUScorer", cssuAlreadyExist);
  if (!cssuAlreadyExist) {
    vector<int> level;
    level.push_back(0);
    G4VPrimitiveScorer* calo = new CalorimeterScorers::PS_Calorimeter("Calorimeter", level, 0);
    G4VPrimitiveScorer* interaction = new InteractionScorers::PS_Interactions("Interaction", ms_InterCoord, 0);
    m_CSSUScorer->RegisterPrimitive(calo);
    m_CSSUScorer->RegisterPrimitive(interaction);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_CSSUScorer);
  }

  bool pmtAlreadyExist = false;
  m_PMTScorer = CheckScorer("CSSUPMTScorer", pmtAlreadyExist);
  if (!pmtAlreadyExist) {
    G4VPrimitiveScorer* pmt = new PHOTODIODESCORERS::PS_PhotoDiode_Rectangle("PhotoDiode", 0, 1, 1, 1, 1);
    m_PMTScorer->RegisterPrimitive(pmt);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_PMTScorer);
  }
}

NPS::VDetector* CSSU::Construct() { return (NPS::VDetector*)new CSSU(); }

extern "C" {
class proxy_nps_CSSU {
 public:
  proxy_nps_CSSU() {
    NPS::DetectorFactory::getInstance()->AddToken("CSSU", "CSSU");
    NPS::DetectorFactory::getInstance()->AddDetector("CSSU", CSSU::Construct);
  }
};

proxy_nps_CSSU p_nps_CSSU;
}
