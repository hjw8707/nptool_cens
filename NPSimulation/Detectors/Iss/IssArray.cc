/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Author: M. Labiche                     address: marc.labiche@stfc.ac.uk    *
 *                                                                           *
 * Creation Date  : July 2019                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This file describe the Iss charge particle Detector                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * Iss is a modular array made of DSSSD (Telescope). Each                    *
 *  Telescope is made of a single Stage:                                     *
 *  - A 300um Silicium, double-sided strip                                   *
 *  - possibility to add a second layer in future                            *
 *****************************************************************************/
#include "IssArray.hh"

// Geant4
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4MaterialTable.hh"
#include "G4PVDivision.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "Randomize.hh"

// For the magnetic field
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4TransportationManager.hh"
#include "MyMagneticField.hh"

// NPS
#include "CalorimeterScorers.hh"
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"
// NPL
#include "NPCore.h"
// ROOT
#include "RootOutput.h"

// CLHEP
#include "CLHEP/Random/RandGauss.h"

// STL
#include <cmath>
#include <set>
#include <sstream>
#include <string>
using namespace std;
using namespace CLHEP;
using namespace ISS;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// IssArray Specific Method
IssArray::IssArray() {
  m_Event = new TIssData();
  InitializeMaterial();
  // if using mask: m_StripScorer = 0;
  m_BOXScorer = 0;
}

IssArray::~IssArray() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void IssArray::AddTelescope(G4ThreeVector X1_Y1, G4ThreeVector X128_Y1, G4ThreeVector X1_Y128, G4ThreeVector X128_Y128,
                            bool wSi) {
  m_DefinitionType.push_back(true);

  m_X1_Y1.push_back(X1_Y1);
  m_X128_Y1.push_back(X128_Y1);
  m_X1_Y128.push_back(X1_Y128);
  m_X128_Y128.push_back(X128_Y128);
  m_wSi.push_back(wSi);

  m_R.push_back(0);
  m_Theta.push_back(0);
  m_Phi.push_back(0);
  m_beta_u.push_back(0);
  m_beta_v.push_back(0);
  m_beta_w.push_back(0);
}

void IssArray::AddTelescope(G4double R, G4double Theta, G4double Phi, G4double beta_u, G4double beta_v, G4double beta_w,
                            bool wSi) {
  G4ThreeVector empty = G4ThreeVector(0, 0, 0);

  m_DefinitionType.push_back(false);

  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_beta_u.push_back(beta_u);
  m_beta_v.push_back(beta_v);
  m_beta_w.push_back(beta_w);
  m_wSi.push_back(wSi);

  m_X1_Y1.push_back(empty);
  m_X128_Y1.push_back(empty);
  m_X1_Y128.push_back(empty);
  m_X128_Y128.push_back(empty);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void IssArray::VolumeMaker(G4int TelescopeNumber, G4double DetTgt, G4ThreeVector MMpos, G4RotationMatrix* MMrot,
                           bool wSi, G4LogicalVolume* world) {

  G4double NbrTelescopes = TelescopeNumber;
  G4String DetectorNumber;
  std::ostringstream Number;
  Number << NbrTelescopes;
  DetectorNumber = Number.str();

  ////////////////////////////////////////////////////////////////
  ////////////// Starting Volume Definition //////////////////////
  ////////////////////////////////////////////////////////////////

  //
  // Mechanical structure to be included only once with the 1st helios module
  //

  const char* ToMyGDML = NULL;
  G4String MyGDMLPath;

  if (getenv("MyGDML")) {
    // const char* ToMyGDML= getenv("MyGDML"); // "/mnt/hgfs/Documents/NPhelisol/gdml/";
    MyGDMLPath = G4String(getenv("MyGDML"));
  }

  if (TelescopeNumber == 1) {
    //
    // Aluminium support plate (gdml)
    //

    if (MyGDMLPath) {

      cout << "Construction using gdml files" << endl;

      m_gdmlparser.Read(MyGDMLPath + "/ISSNudePlate.gdml");
      G4LogicalVolume* m_LogicalAlPlate = m_gdmlparser.GetVolume("StructPlate");

      m_gdmlparser.Read(MyGDMLPath + "/ISSTarg.gdml");
      G4LogicalVolume* m_LogicalSupTarg = m_gdmlparser.GetVolume("StructTarg");

      m_gdmlparser.Read(MyGDMLPath + "/ISSDetStrut.gdml");
      G4LogicalVolume* m_LogicalSupDet = m_gdmlparser.GetVolume("StructDet");

      /*      m_gdmlparser.Read("/mnt/hgfs/Documents/NPhelisol/gdml/ISSNudePlate.gdml");
            G4LogicalVolume* m_LogicalAlPlate= m_gdmlparser.GetVolume("StructPlate");

            m_gdmlparser.Read("/mnt/hgfs/Documents/NPhelisol/gdml/ISSTarg.gdml");
            G4LogicalVolume* m_LogicalSupTarg= m_gdmlparser.GetVolume("StructTarg");

            m_gdmlparser.Read("/mnt/hgfs/Documents/NPhelisol/gdml/ISSDetStrut.gdml");
            G4LogicalVolume* m_LogicalSupDet= m_gdmlparser.GetVolume("StructDet");
      */
      m_LogicalAlPlate->SetVisAttributes(G4VisAttributes::GetInvisible()); // set world of Al Plate invisible
      m_LogicalSupTarg->SetVisAttributes(G4VisAttributes::GetInvisible()); // set world of Al Plate invisible
      m_LogicalSupDet->SetVisAttributes(G4VisAttributes::GetInvisible());  // set world of Al Plate invisible

      // G4VisAttributes *TargVol = new G4VisAttributes(G4Colour(1.0, .5, .5));
      // TargVol->SetForceSolid(false);
      // m_LogicalSupTarg->SetVisAttributes(TargVol);

      // G4VisAttributes *SupDetVol = new G4VisAttributes(G4Colour(1.0, .5, .5));
      // SupDetVol->SetForceWireframe(false);
      // m_LogicalSupDet->SetVisAttributes(SupDetVol);

      G4RotationMatrix* RotY_DetSup = new G4RotationMatrix();
      RotY_DetSup->rotateY(180 * deg);

      cout << "MMpos0: " << MMpos(0) << endl;
      cout << "MMpos1: " << MMpos(1) << endl;
      cout << "MMpos2: " << MMpos(2) << endl;

      if (DetTgt < 0) { // ie: detector at backward angles
        new G4PVPlacement(0, G4ThreeVector(0., -353, 0.), m_LogicalAlPlate, "AlumPlate", world, false, 0);
        new G4PVPlacement(0, G4ThreeVector(13.5, -200. + (27.830127 - 5.66), 87.910256 - 5. + 4.), m_LogicalSupTarg,
                          "SupTarg", world, false, 0);
        new G4PVPlacement(0, G4ThreeVector(0., 0., -6.5 + (DetTgt + 100)), m_LogicalSupDet, "SupDet", world, false, 0);
      }
      else // ie: detector at forward angle
      {
        new G4PVPlacement(RotY_DetSup, G4ThreeVector(0., -370., 0.), m_LogicalAlPlate, "AlumPlate", world, false, 0);
        new G4PVPlacement(RotY_DetSup, G4ThreeVector(-13.5, -200. + (27.830127 - 5.66), -87.910256 + 5. - 4.),
                          m_LogicalSupTarg, "SupTarg", world, false, 0);

        new G4PVPlacement(RotY_DetSup, G4ThreeVector(0., 0., 6.5 + (DetTgt - 100)), m_LogicalSupDet, "SupDet", world,
                          false, 1); // for detecmm
      }

      /*
             new G4PVPlacement(0, G4ThreeVector(0., -370., 0.), m_LogicalAlPlate,"AlumPlate",  world, false, 0 );
             new G4PVPlacement(0, G4ThreeVector(13.5, -200.+(27.830127-5.66), 87.910256-5.+4.),
         m_LogicalSupTarg,"SupTarg",  world, false, 0 ); // for backward angles
             //new G4PVPlacement(RotY_DetSup, G4ThreeVector(0., -370., 0.), m_LogicalAlPlate,"AlumPlate",  world, false,
         0 ); // for forward angles
             //new G4PVPlacement(RotY_DetSup, G4ThreeVector(-13.5, -200.+(27.830127-5.66), -87.910256+5.-4.),
         m_LogicalSupTarg,"SupTarg",  world, false, 0 ); // for forward angles

             //new G4PVPlacement(0, G4ThreeVector(0., 0., -156.5), m_LogicalSupDet,"SupDet",  world, false, 0 ); // for
         detector at -250mm
             //new G4PVPlacement(0, G4ThreeVector(0., 0., -6.5), m_LogicalSupDet,"SupDet",  world, false, 0 ); // for
         detector at -100mm new G4PVPlacement(0, G4ThreeVector(0., 0., -26.5), m_LogicalSupDet,"SupDet",  world, false,
         0 ); // for detector at -120mm
             //new G4PVPlacement(0, G4ThreeVector(0., 0., -106.5), m_LogicalSupDet,"SupDet",  world, false, 0 ); // for
         detector at -200mm
             //new G4PVPlacement(RotY_DetSup, G4ThreeVector(0., 0., 156.5), m_LogicalSupDet,"SupDet",  world, false, 1
         ); // for detector at 156.5mm
             //new G4PVPlacement(RotY_DetSup, G4ThreeVector(0., 0., 36.5), m_LogicalSupDet,"SupDet",  world, false, 1 );
         // for detector at 130.mm

      */
    }
    else {

      //
      //  Detector support structure (rods)
      //

      cout << "Construction without gdml files" << endl;
      cout << "to add gdml geometry files, make sure the environment variable MyGDML in nptool.sh or .csh point to "
              "your gdml directory "
           << endl;

      // Add the Aluminium rod

      G4double RodTgt;

      if (DetTgt < 0)
        RodTgt = -10 - DetTgt;
      if (DetTgt > 0)
        RodTgt = -10 + DetTgt;

      const G4double Zplan[2] = {(RodTgt), (RodTgt + 530)};
      // const G4double Zplan[2]={140*mm, 670*mm}; // if first silicons are at 150mm to target (helisol.detector)
      // const G4double Zplan[2]={90*mm, 620*mm}; // if first silicon at 100mm to target (helisol_2.detector)
      // const G4double Zplan[2]={40*mm, 570*mm}; // if first silicon at 50mm to target (helisol_3.detector)
      const G4double Rinner[2] = {15 * mm, 15 * mm};
      const G4double Router[2] = {24 * mm, 24 * mm};
      G4Polyhedra* Al_rod = new G4Polyhedra("Al_rod", 0., 360 * deg, 6, 2, Zplan, Rinner, Router);

      G4LogicalVolume* Al_rod_log = new G4LogicalVolume(Al_rod, m_MaterialAluminium, "Al_rod", 0, 0, 0);

      G4RotationMatrix* RotZ30Y180 = new G4RotationMatrix();
      RotZ30Y180->rotateY(180 * deg);
      RotZ30Y180->rotateZ(30 * deg);
      G4RotationMatrix* RotZ30 = new G4RotationMatrix();
      RotZ30->rotateZ(30 * deg);

      if (DetTgt < 0) {
        new G4PVPlacement(RotZ30Y180, G4ThreeVector(0., 0., 0.), Al_rod_log, "Al_rod", world, false, 1); // backward
      }
      else {
        new G4PVPlacement(RotZ30, G4ThreeVector(0., 0., 0.), Al_rod_log, "Al_rod", world, false, 1); // forward
      }

      // new G4PVPlacement(RotZ30, G4ThreeVector(0.,0., 100., Al_rod_log, "Al_rod", world, false, 0); // forward
      // new G4PVPlacement(RotZ30Y180, G4ThreeVector(0.,0., -150.), Al_rod_log, "Al_rod", world, false, 1);  // backward

      G4VisAttributes* VisAtt1 = new G4VisAttributes(G4Colour(1, 0.8, 0));
      Al_rod_log->SetVisAttributes(VisAtt1);
    }

    //
    // Add the Aluminium chamber
    //
    G4double Al_chamber_rmin = 50. * cm;
    G4double Al_chamber_rmax = 50.5 * cm;
    G4double Al_chamber_z = 150.0 * cm;

    // G4Tubs* Al_chamber_tub
    //   = new G4Tubs("Al_chamber_tub", Al_chamber_rmin, Al_chamber_rmax, Al_chamber_z, 0.*deg, 180*deg);
    G4Tubs* Al_chamber_tub =
        new G4Tubs("Al_chamber_tub", Al_chamber_rmin, Al_chamber_rmax, Al_chamber_z, 0. * deg, 240 * deg);

    G4LogicalVolume* Al_chamber_log = new G4LogicalVolume(Al_chamber_tub, m_MaterialAluminium, "Al_chamber", 0, 0, 0);

    G4RotationMatrix* RotZ = new G4RotationMatrix();
    RotZ->rotateZ(-90 * deg);

    new G4PVPlacement(RotZ, G4ThreeVector(0., 0., 0.), Al_chamber_log, "Al_chamber", world, false, 0);

    G4VisAttributes* VisAtt2 = new G4VisAttributes(G4Colour(0, 0.8, 0.5));
    // VisAtt2->SetForceWireframe("True");
    Al_chamber_log->SetVisAttributes(VisAtt2);

    //
    // Add the Aluminium shield after target to stop most shallow angle particles
    //
    G4double Al_tgt_rmin = 0.03 * cm;
    G4double Al_tgt_rmax = 0.13 * cm;
    G4double Al_tgt_z = .05 * cm;

    G4Tubs* Al_tgt_shield = new G4Tubs("Al_tgt_shield", Al_tgt_rmin, Al_tgt_rmax, Al_tgt_z, 0. * deg, 360 * deg);

    G4LogicalVolume* Al_tgt_log = new G4LogicalVolume(Al_tgt_shield, m_MaterialAluminium, "Al_tgt_shield", 0, 0, 0);

    // new G4PVPlacement(RotZ, G4ThreeVector(0.,0.,-0.2*cm), Al_tgt_log, "Al_tgt_shield", world, false, 0);

  } // endof of TelescopNumber==1

  ////////////////////////////////////////////////////////////////
  ///////////////// First & single Stage Construction /////////////////////
  ////////////////////////////////////////////////////////////////

  // Mother volume:
  G4Box* solidMM = new G4Box("ISSDet" + DetectorNumber, 0.5 * FaceLength, 0.5 * FaceWidth, 0.5 * Length);
  G4LogicalVolume* logicMM = new G4LogicalVolume(solidMM, m_MaterialVacuum, "ISSDet" + DetectorNumber, 0, 0, 0);
  G4String Name = "ISSDet" + DetectorNumber;

  new G4PVPlacement(G4Transform3D(*MMrot, MMpos), logicMM, Name, world, false, TelescopeNumber);

  if (m_non_sensitive_part_visiualisation) {
    G4VisAttributes* FrameVisAtt = new G4VisAttributes(G4Colour(0.80, 0.80, 0.80));
    FrameVisAtt->SetForceWireframe(true);
    logicMM->SetVisAttributes(FrameVisAtt);
  }
  else
    logicMM->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4ThreeVector positionVacBox = G4ThreeVector(0, 0, VacBox_PosZ);

  G4Box* solidVacBox = new G4Box("solidVacBox", 0.5 * SiliconFaceLength, 0.5 * SiliconFaceWidth, 0.5 * VacBoxThickness);
  G4LogicalVolume* logicVacBox = new G4LogicalVolume(solidVacBox, m_MaterialVacuum, "logicVacBox", 0, 0, 0);

  new G4PVPlacement(0, positionVacBox, logicVacBox, Name + "_VacBox", logicMM, false, TelescopeNumber);

  logicVacBox->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VisAttributes* SiliconVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));

  // Silicon:
  if (wSi) {

    // Aluminium layer

    G4ThreeVector positionAluStripFront = G4ThreeVector(0, 0, AluStripFront_PosZ);
    G4ThreeVector positionAluStripBack = G4ThreeVector(0, 0, AluStripBack_PosZ);

    G4Box* solidAluStrip =
        new G4Box("AluBox", 0.5 * SiliconFaceLength, 0.5 * SiliconFaceWidth, 0.5 * AluStripThickness);
    G4LogicalVolume* logicAluStrip = new G4LogicalVolume(solidAluStrip, m_MaterialAluminium, "logicAluStrip", 0, 0, 0);

    new G4PVPlacement(0, positionAluStripFront, logicAluStrip, Name + "_AluStripFront", logicMM, false,
                      TelescopeNumber);
    new G4PVPlacement(0, positionAluStripBack, logicAluStrip, Name + "_AluStripBack", logicMM, false, TelescopeNumber);

    logicAluStrip->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Silicon detector itself

    G4ThreeVector positionSilicon = G4ThreeVector(0, 0, Silicon_PosZ);

    G4Box* solidSilicon =
        new G4Box("solidSilicon", 0.5 * SiliconFaceLength, 0.5 * SiliconFaceWidth, 0.5 * SiliconThickness);
    G4LogicalVolume* logicSilicon = new G4LogicalVolume(solidSilicon, m_MaterialSilicon, "logicSilicon", 0, 0, 0);

    new G4PVPlacement(0, positionSilicon, logicSilicon, Name + "_Silicon", logicMM, false, TelescopeNumber);

    /// Set Silicon strip sensible
    logicSilicon->SetSensitiveDetector(m_BOXScorer);
    // if using mask: logicSilicon->SetSensitiveDetector(m_StripScorer);

    /// Visualisation of Silicon Strip
    logicSilicon->SetVisAttributes(SiliconVisAtt);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void IssArray::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ISSDet");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detector found" << endl;
  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (NPOptionManager::getInstance()->GetVerboseLevel())
      cout << endl << "//// Iss det " << i + 1 << endl;
    // Cartesian Case
    vector<string> cart = {"X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128", "SI"};
    // Spherical Case
    vector<string> sphe = {"R", "THETA", "PHI", "BETA", "SI"};

    if (blocks[i]->HasToken("MField")) {
      double Bz = blocks[i]->GetDouble("MField", "T");
      MyMagneticField* myField = new MyMagneticField(G4ThreeVector(0., 0., Bz));
      G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
      fieldMgr->SetDetectorField(myField);
      fieldMgr->CreateChordFinder(myField);
    }

    if (blocks[i]->HasTokenList(cart)) {
      G4ThreeVector A = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y1", "mm"));
      G4ThreeVector B = NPS::ConvertVector(blocks[i]->GetTVector3("X128_Y1", "mm"));
      G4ThreeVector C = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y128", "mm"));
      G4ThreeVector D = NPS::ConvertVector(blocks[i]->GetTVector3("X128_Y128", "mm"));
      int SI = blocks[i]->GetInt("SI");

      AddTelescope(A, B, C, D, SI == 1);
    }

    else if (blocks[i]->HasTokenList(sphe)) {

      double Theta = blocks[i]->GetDouble("THETA", "deg");
      double Phi = blocks[i]->GetDouble("PHI", "deg");
      double R = blocks[i]->GetDouble("R", "mm");
      vector<double> beta = blocks[i]->GetVectorDouble("BETA", "deg");
      int SI = blocks[i]->GetInt("SI");
      AddTelescope(R, Theta, Phi, beta[0], beta[1], beta[2], SI == 1);
    }

    else {
      cout << "WARNING: Missing token for ISSDet blocks, check your input "
              "file"
           << endl;
      exit(1);
    }

    if (blocks[i]->GetString("VIS") == "all")
      m_non_sensitive_part_visiualisation = true;
  }
}

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void IssArray::ConstructDetector(G4LogicalVolume* world) {
  G4RotationMatrix* MMrot = NULL;
  G4ThreeVector MMpos = G4ThreeVector(0, 0, 0);
  G4ThreeVector MMu = G4ThreeVector(0, 0, 0);
  G4ThreeVector MMv = G4ThreeVector(0, 0, 0);
  G4ThreeVector MMw = G4ThreeVector(0, 0, 0);
  G4ThreeVector MMCenter = G4ThreeVector(0, 0, 0);
  G4double DetTgt = 0;
  bool Si = true;

  G4int NumberOfTelescope = m_DefinitionType.size();

  for (G4int i = 0; i < NumberOfTelescope; i++) {

    if (i == 0) { // first telescope found
      MMu = m_X1_Y1[i];
      DetTgt = MMu(2);
    }

    // By Point
    if (m_DefinitionType[i]) {
      // (u,v,w) unitary vector associated to telescope referencial
      // (u,v) // to silicon plan
      // w perpendicular to (u,v) plan and pointing CsI
      MMu = m_X128_Y1[i] - m_X1_Y1[i];
      MMu = MMu.unit();

      MMv = m_X1_Y128[i] - m_X1_Y1[i];
      MMv = MMv.unit();

      MMw = MMv.cross(MMu);
      // if (MMw.z() > 0)MMw = MMv.cross(MMu)  ;
      MMw = MMw.unit();

      MMCenter = (m_X1_Y1[i] + m_X1_Y128[i] + m_X128_Y1[i] + m_X128_Y128[i]) / 4;

      // Passage Matrix from Lab Referential to Telescope Referential
      MMrot = new G4RotationMatrix(MMv, MMu, MMw);
      MMpos = MMw * Length * 0.5 + MMCenter;
    }

    // By Angle
    else {
      G4double Theta = m_Theta[i];
      G4double Phi = m_Phi[i];

      // (u,v,w) unitary vector associated to telescope referencial
      // (u,v) // to silicon plan
      // w perpendicular to (u,v) plan and pointing ThirdStage
      // Phi is angle between X axis and projection in (X,Y) plan
      // Theta is angle between  position vector and z axis
      G4double wX = m_R[i] * sin(Theta / rad) * cos(Phi / rad);
      G4double wY = m_R[i] * sin(Theta / rad) * sin(Phi / rad);
      G4double wZ = m_R[i] * cos(Theta / rad);
      MMw = G4ThreeVector(wX, wY, wZ);

      // vector corresponding to the center of the module
      G4ThreeVector CT = MMw;

      // vector parallel to one axis of silicon plane
      G4double ii = cos(Theta / rad) * cos(Phi / rad);
      G4double jj = cos(Theta / rad) * sin(Phi / rad);
      G4double kk = -sin(Theta / rad);
      G4ThreeVector Y = G4ThreeVector(ii, jj, kk);

      MMw = MMw.unit();
      MMu = MMw.cross(Y);
      MMv = MMw.cross(MMu);
      MMv = MMv.unit();
      MMu = MMu.unit();

      // Passage Matrix from Lab Referential to Telescope Referential
      // Iss
      MMrot = new G4RotationMatrix(MMu, MMv, MMw);
      // Telescope is rotate of Beta angle around MMv axis.
      MMrot->rotate(m_beta_u[i], MMu);
      MMrot->rotate(m_beta_v[i], MMv);
      MMrot->rotate(m_beta_w[i], MMw);
      // translation to place Telescope
      MMpos = MMw * Length * 0.5 + CT;
    }

    Si = m_wSi[i];

    VolumeMaker(i + 1, DetTgt, MMpos, MMrot, Si, world);
  }

  delete MMrot;
}

// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method

void IssArray::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree* pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("ISS")) {
    pTree->Branch("ISS", "TIssData", &m_Event);
  }
  pTree->SetBranchAddress("ISS", &m_Event);
}

// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void IssArray::ReadSensitive(const G4Event*) {
  m_Event->Clear();

  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////// Used to Read Event Map of detector
  /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////

  /////////////////////
  // Read the Scorer associate to the Silicon Strip in case of PS_Images (See NPLib/Detectors/Iss/ressources)
  ///////////
  // BOX
  DSSDScorers::PS_Rectangle* BOXScorer = (DSSDScorers::PS_Rectangle*)m_BOXScorer->GetPrimitive(0);

  // Loop on the BOX map
  unsigned int sizeFront = BOXScorer->GetLengthMult();

  for (unsigned int i = 0; i < sizeFront; i++) {

    double Energy = BOXScorer->GetEnergyLength(i);

    if (Energy > ThresholdSi) {
      double Time = BOXScorer->GetTimeLength(i);
      int DetNbr = BOXScorer->GetDetectorLength(i);
      int StripFront = BOXScorer->GetStripLength(i);

      m_Event->SetFront(DetNbr, StripFront, RandGauss::shoot(Energy, ResoStrip), RandGauss::shoot(Time, ResoTime),
                        RandGauss::shoot(Time, ResoTime));
    }
  }

  unsigned int sizeBack = BOXScorer->GetWidthMult();
  for (unsigned int i = 0; i < sizeBack; i++) {

    double Energy = BOXScorer->GetEnergyWidth(i);

    if (Energy > ThresholdSi) {
      double Time = BOXScorer->GetTimeWidth(i);
      int DetNbr = BOXScorer->GetDetectorWidth(i);
      int StripBack = BOXScorer->GetStripWidth(i);

      m_Event->SetBack(DetNbr, Silicon_Back_NumberOfStrip - StripBack + 1, RandGauss::shoot(Energy, ResoStrip),
                       RandGauss::shoot(Time, ResoTime), RandGauss::shoot(Time, ResoTime));
    }
  }
  // clear map for next event
  BOXScorer->clear();

  /////////////////////
  // Read the Scorer associate to the Silicon Strip in case of using mask ( See NPLib/Detectors/Iss/ressources)
  /* DSSDScorers::PS_Images* SiScorer
     = (DSSDScorers::PS_Images*)m_StripScorer->GetPrimitive(0);

   bool     SiScoredHit; // flag true if first stage scores a hit above threshold
   set<int> trig; // list of telescope that got a Si trigger
   unsigned int sizeFront = SiScorer->GetFrontMult();
   unsigned int sizeBack  = SiScorer->GetBackMult();

   // Check for double match Strip :
   // rare case where a particle hit a strip and then an interstrip
   // since the map idex is build on pixel value, we end up with the same strip
   // fired twice, which is impossible in reality.
   std::map< unsigned int, std::pair<double,double> > mapFront;
   std::map< unsigned int, std::pair<double,double> >::iterator it;

   for (unsigned int i = 0; i < sizeFront; i++) {
     double energy      = SiScorer->GetEnergyFront(i);
     int    detectorNbr = SiScorer->GetDetectorFront(i);
     double time        = SiScorer->GetTimeFront(i);
     // Pixel value at interaction point
     unsigned int a, r, g, b;
     //  pixel
     SiScorer->GetARGBFront(i, a, r, g, b);
     if (r == 0) {
       mapFront[b+detectorNbr*1e6].first+=energy;
       mapFront[b+detectorNbr*1e6].second=time;
     }

     else { // Interstrip X, keep maximum shared energy
       double rand = G4UniformRand();
       if (rand > 0.5) {
         double energyX = rand * energy;
           mapFront[b+detectorNbr*1e6].first+=energyX;
           mapFront[b+detectorNbr*1e6].second=time;
         }

       else {
           double energyX = (1 - rand) * energy;
           mapFront[g+detectorNbr*1e6].first+=energyX;
           mapFront[g+detectorNbr*1e6].second=time;
         }
       }
     }

   for(it=mapFront.begin();it!=mapFront.end();it++){
     double energyX = RandGauss::shoot(it->second.first, ResoStrip);
     double timeX = TimeOffset - RandGauss::shoot(it->second.second, ResoTimeMust);
     unsigned int strip = it->first-1000000*(it->first/1000000);
     unsigned int det   = it->first/1000000;
     if (energyX > ThresholdSi) {
       trig.insert(det);
       SiScoredHit = true;
       m_Event->SetStripXE(det, strip ,
           NPL::EnergyToADC(energyX, 0, 63, 8192, 16384));
       m_Event->SetStripXT(det, strip ,
           NPL::EnergyToADC(timeX, 0, 1000, 8192, 16384));
     }
   }

   // Check for double match Strip :
   // rare case where a particle hit a strip and then an interstrip
   // since the map idex is build on pixel value, we end up with the same strip
   // fired twice, which is impossible in reality.
   std::map< unsigned int, std::pair<double,double> > mapBack;

   for (unsigned int i = 0; i < sizeBack; i++) {
     double energy      = SiScorer->GetEnergyBack(i);
     int    detectorNbr = SiScorer->GetDetectorBack(i);
     double time        = SiScorer->GetTimeBack(i);

       // Pixel value at interaction point
       unsigned int a, r, g, b;
       //  pixel
       SiScorer->GetARGBBack(i, a, r, g, b);
       if (r == 0) {
           mapBack[b+detectorNbr*1e6].first+=energy;
           mapBack[b+detectorNbr*1e6].second=time;
       }
       else { // Interstrip Y, keep both strip with shared energy
         double rand     = G4UniformRand();
         double energyY1 = rand * energy;
           mapBack[b+detectorNbr*1e6].first+=energyY1;
           mapBack[b+detectorNbr*1e6].second=time;

         double energyY2 = (1 - rand) * energy;
         mapBack[g+detectorNbr*1e6].first+=energyY2;
         mapBack[g+detectorNbr*1e6].second=time;
         }
       }

   for(it=mapBack.begin();it!=mapBack.end();it++){
     double energyY = RandGauss::shoot(it->second.first, ResoStrip);
     double timeX = TimeOffset - RandGauss::shoot(it->second.second, ResoTimeMust);
     unsigned int strip = it->first-1000000*(it->first/1000000);
     unsigned int det   = it->first/1000000;
     if (energyY > ThresholdSi) {
       trig.insert(det);
       SiScoredHit = true;
       m_Event->SetStripYE(det, strip ,
           NPL::EnergyToADC(energyY, 0, 63, 8192, 0));
       m_Event->SetStripYT(det, strip ,
           NPL::EnergyToADC(timeX, 0, 1000, 8192, 16384));
     }
   }
 */

  // Look for 2nd and 3rd stage only if 1st stage is hit
  /*  if (SiScoredHit) {
    }
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void IssArray::InitializeScorers() {
  //	Silicon Associate Scorer

  bool already_exist = false;
  // m_StripScorer      = CheckScorer("ISS_StripScorer", already_exist);
  m_BOXScorer = CheckScorer("ISS_BOXScorer", already_exist);

  // if the scorer were created previously nothing else need to be made
  if (already_exist) {
    cout << "SCORER already exist" << endl;
    return;
  }

  string nptool = getenv("NPTOOL");

  G4VPrimitiveScorer* BOXScorer = new DSSDScorers::PS_Rectangle(
      "IssBOX", 0, SiliconFaceLength, SiliconFaceWidth, Silicon_Front_NumberOfStrip, Silicon_Back_NumberOfStrip);

  /* is using mask:
   G4VPrimitiveScorer* SiScorer = new DSSDScorers::PS_Images(
        "SiScorer", nptool + "/NPLib/Detectors/Iss/ressources/maskFront.png",
        nptool + "/NPLib/Detectors/Iss/ressources/maskBack.png", 127./12800, 19.9/1100, 0,
        0, 0xffff0000, 0);
  */

  G4VPrimitiveScorer* InterScorer = new InteractionScorers::PS_Interactions("SiScorer", ms_InterCoord, 0);

  // and register it to the multifunctionnal detector
  m_BOXScorer->RegisterPrimitive(BOXScorer);
  m_BOXScorer->RegisterPrimitive(InterScorer);
  // if using mask:
  //  m_StripScorer->RegisterPrimitive(SiScorer);
  //  m_StripScorer->RegisterPrimitive(InterScorer);

  //	Add All Scorer to the Global Scorer Manager
  G4SDManager::GetSDMpointer()->AddNewDetector(m_BOXScorer);
  // G4SDManager::GetSDMpointer()->AddNewDetector(m_StripScorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void IssArray::InitializeMaterial() {

  m_MaterialSilicon = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  m_MaterialAluminium = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  m_MaterialIron = MaterialManager::getInstance()->GetMaterialFromLibrary("Fe");
  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4RotationMatrix* Rotation(double tetaX, double tetaY, double tetaZ) {
  double PI = 3.141592653589793238;
  double radX = tetaX * PI / 180.;
  double radY = tetaY * PI / 180.;
  double radZ = tetaZ * PI / 180.;

  G4ThreeVector col1 = G4ThreeVector(cos(radZ) * cos(radY), -sin(radZ) * cos(radY), -sin(radY));
  G4ThreeVector col2 = G4ThreeVector(sin(radZ) * cos(radX) - cos(radZ) * sin(radY) * sin(radX),
                                     cos(radZ) * cos(radX) + sin(radZ) * sin(radY) * sin(radX), -cos(radY) * sin(radX));
  G4ThreeVector col3 = G4ThreeVector(sin(radZ) * sin(radX) + cos(radZ) * sin(radY) * sin(radX),
                                     cos(radZ) * sin(radX) - sin(radZ) * sin(radY) * cos(radX), cos(radY) * cos(radX));

  return (new G4RotationMatrix(col1, col2, col3));
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* IssArray::Construct() { return (NPS::VDetector*)new IssArray(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_nps_iss {
 public:
  proxy_nps_iss() {
    NPS::DetectorFactory::getInstance()->AddToken("ISSDet", "ISS");
    NPS::DetectorFactory::getInstance()->AddDetector("ISSDet", IssArray::Construct);
  }
};

proxy_nps_iss p_nps_iss;
}
