/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : January 2009                                             *
 * Last update    : 26/08/2010  (M. Labiche)                                 *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This Class manage the virtual detector and call their method.            *
 *  Every time a new detector is added to the project, a few line in the     *
 *  ReadConfigurationFile method are needed in order to detect the associate *
 *  token.                                                                   *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    15/01/10: Added Chamber (marc.labiche@stfc.ac.uk)                      *
 *    04/12/09: Added PARIS detector (marc.labiche@stfc.ac.uk)               *
 *****************************************************************************/

#include "DetectorConstruction.hh"

// G4
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polyhedron.hh"   // exact triangulation of Boolean/other solids on export
#include "G4Point3D.hh"
#include "G4BooleanSolid.hh" // fall back to a constituent when a Boolean can't tessellate
#include "G4RotationMatrix.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"

#ifdef NPS_GDML
#include "G4GDMLParser.hh"
#endif

// STL
#include <cstdlib>
#include <fstream>
#include <set>
#include <string>
// NPL
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "RootOutput.h"

// NPS
#include "DetectorMessenger.hh"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction() {
  world_log = NULL;
  world_phys = NULL;
  world_box = NULL;
  m_Target = NULL;
  m_Chamber = NULL;
  m_Messenger = new DetectorMessenger(this);
  m_ReadSensitivePtr = &NPS::VDetector::ReadSensitive;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct() { return ReadConfigurationFile(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddDetector(NPS::VDetector* NewDetector) {
  // Add new detector to vector
  m_Detectors.push_back(NewDetector);

  // Initialize Scorer
  NewDetector->InitializeScorers();

  // Construct new detector
  NewDetector->ConstructDetector(world_log);

  // Add Detector to TTree
  NewDetector->InitializeRootOutput();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::ReadConfigurationFile() {
  // Construct the World volume
  G4Material* Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");

  //------------------------------world volume
  G4double world_x = 100.0 * m;
  G4double world_y = 100.0 * m;
  G4double world_z = 100.0 * m;

  if (!world_box)
    world_box = new G4Box("world_box", world_x, world_y, world_z);

  if (!world_log)
    world_log = new G4LogicalVolume(world_box, Vacuum, "world_log", 0, 0, 0);

  if (!world_phys)
    world_phys = new G4PVPlacement(0, G4ThreeVector(), world_log, "world", 0, false, 0);

  G4VisAttributes* VisAtt = new G4VisAttributes(G4VisAttributes::GetInvisible());
  world_log->SetVisAttributes(VisAtt);

  //------------------------------------------------------------------

  //------------------------------------------------------------------

  std::string Path = NPOptionManager::getInstance()->GetDetectorFile();

  ////////General Reading needs////////
  std::string LineBuffer;
  std::string DataBuffer;
  static bool already = false;
  std::set<std::string> check;
  int VerboseLevel = 0;
  if (!already) {
    VerboseLevel = NPOptionManager::getInstance()->GetVerboseLevel();
    already = 1;
  }
  else
    VerboseLevel = 0;

  if (VerboseLevel)
    cout << endl << "\033[1;36m//// Reading detector file  " << Path << endl << endl;

  if (already)
    cout << endl << "\033[1;36m//// Changing detector file to  " << Path << endl << endl;

  // Access the DetectorFactory and ask it to load the Class List
  std::string classlist = getenv("NPTOOL");
  classlist += "/NPLib/ClassList.txt";
  NPS::DetectorFactory* theFactory = NPS::DetectorFactory::getInstance();
  theFactory->ReadClassList(classlist);

  NPL::InputParser parser(Path);
  ////////////////////////////////////////////
  /////////// Search for Target  /////////////
  ////////////////////////////////////////////
  std::vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Target");
  if (blocks.size() == 1) {
    m_Target = new Target();
    m_Target->ReadConfiguration(parser);
    AddDetector(m_Target);
  }
  else {
    blocks = parser.GetAllBlocksWithToken("CryogenicTarget");
    if (blocks.size() == 1) {
      m_Target = new Target();
      m_Target->ReadConfiguration(parser);
      AddDetector(m_Target);
    }
    else {
      cout << "WARNING: No target found in detector input file! Use with caution" << endl;
    }
  }
  ////////////////////////////////////////////
  /////////// Search for Chamber /////////////
  ////////////////////////////////////////////
  blocks.clear();
  blocks = parser.GetAllBlocksWithToken("Chamber");
  if (blocks.size() == 1) {
    m_Chamber = new Chamber();
    m_Chamber->ReadConfiguration(parser);
    AddDetector(m_Chamber);
  }
  ////////////////////////////////////////////
  /////////// Search for Detectors ///////////
  ////////////////////////////////////////////
  // Get the list of main token
  std::vector<std::string> token = parser.GetAllBlocksToken();
  // Look for detectors among them
  for (unsigned int i = 0; i < token.size(); i++) {
    NPS::VDetector* detector = theFactory->Construct(token[i]);
    if (detector != NULL && check.find(token[i]) == check.end()) {
      cout << "/////////////////////////////////////////" << endl;
      cout << "//// Adding Detector " << token[i] << endl;
      detector->ReadConfiguration(parser);
      cout << "/////////////////////////////////////////" << endl;
      // Add array to the VDetector Vector
      AddDetector(detector);
      check.insert(token[i]);
    }
    else if (detector != NULL)
      delete detector;
  }

  cout << "\033[0m";

  // Create the Material sample for DEDX tables
  MaterialManager::getInstance()->CreateSampleVolumes(world_log);
  return world_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ReadAllSensitive(const G4Event* event) {
  unsigned int mysize = m_Detectors.size();
  for (unsigned short i = 0; i < mysize; i++) {
    (m_Detectors[i]->*m_ReadSensitivePtr)(event);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* DetectorConstruction::GetWorldLogic() { return world_log; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ClearGeometry() {

  unsigned int mySize = m_Detectors.size();

  for (unsigned short i = 0; i < mySize; i++) {
    delete m_Detectors[i];
  }

  m_Detectors.clear();
  // Delete previous Sensitive Detector (Scorer)
  delete G4SDManager::GetSDMpointer();

  G4GeometryManager::GetInstance()->OpenGeometry();

  // Cleaning stores by hand (everything but the world)
  // Cleaning the PVP
  G4PhysicalVolumeStore* PVP = G4PhysicalVolumeStore::GetInstance();
  while (PVP->size() > 1) {
    delete PVP->at(PVP->size() - 1);
    PVP->pop_back();
  }
  // Cleaning the Logical Volume
  G4LogicalVolumeStore* LV = G4LogicalVolumeStore::GetInstance();
  while (LV->size() > 1) {
    delete LV->at(LV->size() - 1);
    LV->pop_back();
  }
  // Cleaning the Solid store
  G4SolidStore* SL = G4SolidStore::GetInstance();
  while (SL->size() > 1) {
    delete SL->at(SL->size() - 1);
    SL->pop_back();
  }

  G4Region* reg = G4RegionStore::GetInstance()->GetRegion("NPSimulationProcess", false);
  if (reg)
    G4RegionStore::GetInstance()->DeRegister(reg);

  delete reg;

  m_Target = 0;
  m_Chamber = 0;
  world_log->ClearDaughters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::RedefineGeometry(std::string file) {
  ClearGeometry();

  if (file != "")
    NPOptionManager::getInstance()->SetDetectorFile(file);
  Construct();

  if (m_Target)
    m_Target->SetReactionRegion();

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  G4RunManager::GetRunManager()->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ExportGeometry(string file) {
#ifdef NPS_GDML
  G4GDMLParser parser;
  parser.Write(file.c_str(), world_log);
#else
  // To silenced the compiler
  file = "";
  G4cout << "You need to compile Geant4 with GDML support to use this command" << G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Walk the placed physical-volume tree and dump each solid with its GLOBAL
// transform (rotation matrix + translation). Unlike GDML this handles nested
// assemblies (they are already imprinted into real placements here) and carries
// the exact rotation Geant4 uses at tracking time — no Euler-convention drift.
static void DumpPV(G4VPhysicalVolume* pv, const G4RotationMatrix& Rp,
                   const G4ThreeVector& Tp, std::ofstream& out, bool& first) {
  G4RotationMatrix Robj = pv->GetObjectRotationValue();  // local -> mother
  G4ThreeVector    Tobj = pv->GetObjectTranslation();
  G4RotationMatrix Rg = Rp * Robj;                       // local -> world
  G4ThreeVector    Tg = Rp * Tobj + Tp;

  G4LogicalVolume* lv = pv->GetLogicalVolume();
  G4VSolid* solid = lv ? lv->GetSolid() : nullptr;
  std::string type; double a = 0, b = 0, c = 0;
  double rmin = 0, phi0 = 0, dphi = 0;   // extra Tubs params (inner radius + phi wedge)
  std::vector<double> verts;             // "mesh" type: triangles (local coords), 9 per tri
  if (auto* box = dynamic_cast<G4Box*>(solid)) {
    type = "box"; a = box->GetXHalfLength() * 2; b = box->GetYHalfLength() * 2; c = box->GetZHalfLength() * 2;
  } else if (auto* tub = dynamic_cast<G4Tubs*>(solid)) {
    type = "cylinder"; a = tub->GetOuterRadius(); c = tub->GetZHalfLength() * 2;
    // carry inner radius + phi start/sweep so a partial Tubs (e.g. a plain annular
    // sector) renders as the real fan shape, not a full disc.
    rmin = tub->GetInnerRadius(); phi0 = tub->GetStartPhiAngle(); dphi = tub->GetDeltaPhiAngle();
  } else if (auto* sph = dynamic_cast<G4Sphere*>(solid)) {
    type = "sphere"; a = sph->GetOuterRadius();
  } else if (solid) {
    // Any other solid — a Boolean (Subtraction/Union) or Trd/Cons/Polycone/… — has
    // no single primitive shape. Tessellate the REAL solid via its polyhedron so the
    // exact shape (incl. boolean cuts, e.g. the STARK QQQ5 wafer's rectangular
    // notches) is exported as a triangle mesh in the solid's local frame.
    // Some Boolean combos can't build a polyhedron (G4 prints "No G4Polyhedron for
    // Boolean component", e.g. the Plunger sphere-shell + pipe union). When that
    // happens, fall back to the first constituent so at least the base shape shows.
    G4VSolid* tess = solid;
    G4Polyhedron* poly = tess->GetPolyhedron();
    for (int guard = 0; !poly && guard < 8; ++guard) {
      auto* boolean = dynamic_cast<G4BooleanSolid*>(tess);
      if (!boolean) break;
      tess = boolean->GetConstituentSolid(0);
      poly = tess->GetPolyhedron();
    }
    if (poly) {
      const G4int nf = poly->GetNoFacets();
      auto push = [&](const G4Point3D& p) { verts.push_back(p.x()); verts.push_back(p.y()); verts.push_back(p.z()); };
      for (G4int f = 0; f < nf; ++f) {
        G4int nnode = 0; G4Point3D nd[4];
        poly->GetNextFacet(nnode, nd);
        if (nnode < 3) continue;
        push(nd[0]); push(nd[1]); push(nd[2]);              // triangle
        if (nnode == 4) { push(nd[0]); push(nd[2]); push(nd[3]); }  // split quad
      }
      if (!verts.empty()) type = "mesh";
    }
  }

  double cr = 0.6, cg = 0.7, cb = 0.8;
  const G4VisAttributes* vis = lv ? lv->GetVisAttributes() : nullptr;
  if (vis) { G4Colour col = vis->GetColour(); cr = col.GetRed(); cg = col.GetGreen(); cb = col.GetBlue(); }
  // skip volumes the geometry marks invisible (e.g. MaterialManager's zero-size
  // "sample_log" material-registration dummies parked 10 m away)
  bool visible = !(vis && !vis->IsVisible());

  if (!type.empty() && visible) {
    if (!first) out << ",\n";
    first = false;
    out << "{\"name\":\"" << lv->GetName() << "\",\"type\":\"" << type << "\",\"copy\":" << pv->GetCopyNo();
    if (type == "mesh") {     // exact triangulated shape (local coords), 9 numbers/triangle
      out << ",\"verts\":[";
      for (size_t k = 0; k < verts.size(); ++k) { if (k) out << ","; out << verts[k]; }
      out << "]";
    } else {
      out << ",\"dim\":[" << a << "," << b << "," << c << "]";
      if (type == "cylinder")   // inner radius + phi wedge (radians) for annular sectors
        out << ",\"rmin\":" << rmin << ",\"phi0\":" << phi0 << ",\"dphi\":" << dphi;
    }
    out << ",\"pos\":[" << Tg.x() << "," << Tg.y() << "," << Tg.z() << "]"
        << ",\"rot\":[" << Rg.xx() << "," << Rg.xy() << "," << Rg.xz() << ","
        << Rg.yx() << "," << Rg.yy() << "," << Rg.yz() << ","
        << Rg.zx() << "," << Rg.zy() << "," << Rg.zz() << "]"
        << ",\"color\":[" << cr << "," << cg << "," << cb << "]}";
  }
  for (int i = 0; i < lv->GetNoDaughters(); ++i) DumpPV(lv->GetDaughter(i), Rg, Tg, out, first);
}

void DetectorConstruction::ExportGeometryJSON(string file) {
  std::ofstream out(file.c_str());
  if (!out) { G4cout << "export_geometry: cannot open " << file << G4endl; return; }
  // Finer polyhedra so tessellated arcs (annular sectors etc.) read smooth.
  G4Polyhedron::SetNumberOfRotationSteps(48);
  out << "[\n";
  bool first = true;
  if (world_phys) {
    G4RotationMatrix Rid;         // identity
    G4ThreeVector Tid(0, 0, 0);
    G4LogicalVolume* wl = world_phys->GetLogicalVolume();
    for (int i = 0; i < wl->GetNoDaughters(); ++i)   // skip the world box itself
      DumpPV(wl->GetDaughter(i), Rid, Tid, out, first);
  }
  out << "\n]\n";
  out.close();
  G4cout << "export_geometry: wrote " << file << G4endl;
}

