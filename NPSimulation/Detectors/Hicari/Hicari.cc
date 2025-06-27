/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Freddy Flavigny  contact: flavigny@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  : april 2022                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class allows to build the HiCARI Ge array geometry                  *
 *  it is largely based on HICARI/GRETINA/AGATA simulation packages          *
 *  Since it is a bit generic, it could be used for other geometries/arrays  *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment: WORK IN PROGRESS, Lots of things still missing                   *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <cmath>
#include <limits>
#include <sstream>
// G4 Geometry object
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Para.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"

#include "G4IntersectionSolid.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

// G4 sensitive
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"

// G4 various object
#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"

// NPTool header
#include "Hicari.hh"
//#include "CalorimeterScorers.hh"
#include "GeScorers.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"
#include "NPSHitsMap.hh"
#include "RootOutput.h"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Hicari_NS {
  // Energy and time Resolution
  // const double EnergyThreshold = 10 * keV;
  // const double ResoTime = 4.5*ns ;  //not used
  // const double ResoEnergy = 2. * keV;
} // namespace Hicari_NS
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Hicari Specific Method
Hicari::Hicari() {
  m_Event = new THicariData();
  m_HicariScorer = 0;

  // InitializeMaterials();

  G4String path = "./";
  if (path.find("./", 0) != string::npos) {
    G4int position = path.find("./", 0);
    if (position == 0)
      path.erase(position, 2);
  }
  iniPath = path;

  cryostatStatus = false;
  // Slot 0 Position (starting point for Euler angle rotations)
  cryostatPos0.setX(0.);
  cryostatPos0.setY(0.);
  cryostatPos0.setZ(406.);
  // approximate polygonal feed-through space with cylinder
  G4double extend = 93. * mm;
  cryostatPos0.setMag(cryostatPos0.mag() - extend);
  cryostatZplanes[0] = 0.;
  cryostatZplanes[1] = extend;
  cryostatZplanes[2] = 145. * mm + extend;
  cryostatZplanes[3] = 155. * mm + extend;
  cryostatZplanes[4] = 217. * mm + extend;
  cryostatZplanes[5] = 217. * mm + extend;
  cryostatZplanes[6] = 232. * mm + extend;
  cryostatRinner[0] = 0.;
  cryostatRinner[1] = 0.;
  cryostatRinner[2] = 0.;
  cryostatRinner[3] = 0.;
  cryostatRinner[4] = 0.;
  cryostatRinner[5] = 0.;
  cryostatRinner[6] = 0.;
  cryostatRouter[0] = 100. * mm;
  cryostatRouter[1] = 130. * mm;
  cryostatRouter[2] = 130. * mm;
  cryostatRouter[3] = 140. * mm;
  cryostatRouter[4] = 140. * mm;
  cryostatRouter[5] = 165. * mm;
  cryostatRouter[6] = 165. * mm;

  readOut = false;

  matCryst = NULL;
  matWalls = NULL;
  matBackWalls = NULL;
  matHole = NULL;
  matCryo = NULL;
  matCrystName = "Germanium";
  matWallsName = "Al";
  // matBackWallsName   = "BackWallMaterial";
  matBackWallsName = "Al";
  // matHoleName        = "G4_Galactic";
  matHoleName = "Vacuum";
  matCryoName = "Al";

  nEuler = 0;
  eulerFile = iniPath + "aeuler";
  nPgons = 0;
  solidFile = iniPath + "asolid";
  nWalls = 0;
  wallsFile = iniPath + "awalls";
  nClAng = 0;
  clustFile = iniPath + "aclust";
  sliceFile = iniPath + "aslice";

  nDets = 0;
  iCMin = 0;
  iGMin = 0;
  maxSec = 1;
  maxSli = 1;

  arrayRmin = 0.;
  arrayRmax = 0.;
  thetaShift = 0.;
  phiShift = 0.;
  thetaPrisma = 0.;
  posShift = G4ThreeVector();

  useCylinder = true;
  usePassive = true;
  drawReadOut = false;
  makeCapsule = false;
  maxSolids = 0;
  totSegments = 0;
  stepFactor = 1;
  stepHasChanged = false;
  printVolumes = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Hicari::~Hicari() {
  clust.clear();
  pgons.clear();
  euler.clear();
  walls.clear();
  nSegments.clear();
  tSegments.clear();
  pgSegLl.clear();
  pgSegLu.clear();
  pgSegRl.clear();
  pgSegRu.clear();
  segVolume.clear();
  segCenter.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int Hicari::InitializeMaterials() {
  /*
    m_Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    m_Aluminum = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    m_Copper = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");

    m_BGO = new G4Material("BGO", 7.13*g/cm3, 3, kStateSolid);  //BGO does not exist in nptool !!
    m_BGO->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("Bi"),4);
    m_BGO->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("Ge"),3);
    m_BGO->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("O"),12);

    m_CsI = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");
  */

  // search the material by its name
  // G4Material* ptMaterial = G4Material::GetMaterial(matCrystName);
  // G4Material* ptMaterial = G4Material::GetMaterial("");
  G4Material* ptMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(matCrystName);
  if (ptMaterial) {
    matCryst = ptMaterial;
    G4cout << "\n----> The crystals material is " << matCryst->GetName() << G4endl;
  }
  else {
    G4cout << " Could not find the material " << matCrystName << G4endl;
    G4cout << " Could not build the array! " << G4endl;
    return 1;
  }

  // search the material by its name
  // ptMaterial = G4Material::GetMaterial(matWallsName);
  ptMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(matWallsName);
  if (ptMaterial) {
    matWalls = ptMaterial;
    G4cout << "\n----> The wall material is " << matWalls->GetName() << G4endl;
  }
  else {
    G4cout << " Could not find the material " << matWallsName << G4endl;
    G4cout << " Could not build the walls! " << G4endl;
  }

  // search the material by its name
  // ptMaterial = G4Material::GetMaterial(matBackWallsName);
  ptMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(matBackWallsName);
  if (ptMaterial) {
    matBackWalls = ptMaterial;
    G4cout << "\n----> The back wall material is " << matBackWalls->GetName() << G4endl;
  }
  else {
    G4cout << " Could not find the material " << matBackWallsName << G4endl;
    G4cout << " Could not build the walls behind the crystals! " << G4endl;
  }

  // ptMaterial = G4Material::GetMaterial(matHoleName);
  ptMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(matHoleName);
  if (ptMaterial) {
    matHole = ptMaterial;
    G4cout << "\n----> The hole material is " << matHole->GetName() << G4endl;
  }
  else {
    G4cout << " Could not find the material " << matHoleName << G4endl;
    G4cout << " Could not build the capsules! " << G4endl;
  }

  // ptMaterial = G4Material::GetMaterial(matCryoName);
  ptMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(matCryoName);
  if (ptMaterial) {
    matCryo = ptMaterial;
    G4cout << "\n----> The cryostat material is " << matCryo->GetName() << G4endl;
  }
  else {
    G4cout << " Could not find the material " << matCryoName << G4endl;
    G4cout << " Could not build the cryostats! " << G4endl;
  }

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// void Hicari::BuildClover(int i_clo, G4LogicalVolume* world)
//{

//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// void Hicari::AddDetector(double X, double Y, double Z, double ThetaX, double ThetaY, double ThetaZ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Hicari::ReadConfiguration(NPL::InputParser parser) {

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Hicari");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " block found " << endl;

  vector<string> filenames = {"asolid", "aclust", "aeuler", "awalls", "aslice"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(filenames)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Hicari " << i + 1 << endl;
      solidFile = blocks[i]->GetString("asolid");
      clustFile = blocks[i]->GetString("aclust");
      eulerFile = blocks[i]->GetString("aeuler");
      wallsFile = blocks[i]->GetString("awalls");
      sliceFile = blocks[i]->GetString("aslice");
      // AddDetector(X,Y,Z,ThetaX, ThetaY, ThetaZ);
    }
    else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  ReadSolidFile();
  ReadClustFile();
  ReadWallsFile();
  ReadEulerFile();
  ReadSliceFile();
}

/////////////////////////////////////////////////////////////
///////////////// methods to read the files
/////////////////////////////////////////////////////////////

void Hicari::ReadSolidFile() {
  FILE* fp;
  char line[256];
  G4int lline, i1, i2, i3, nvdots, opgon;
  float x, y, z, X, Y, Z;

  nPgons = 0;
  nDets = 0;
  nClus = 0;

  if ((fp = fopen(solidFile, "r")) == NULL) {
    G4cout << "\nError opening data file " << solidFile << G4endl;
    exit(-1);
  }

  G4cout << "\nReading description of crystals from file " << solidFile << " ..." << G4endl;

  pgons.clear();

  nvdots = 0;
  opgon = -1;
  maxPgons = -1;
  CpolyhPoints* pPg = NULL;

  while (fgets(line, 255, fp) != NULL) {
    lline = strlen(line);
    if (lline < 2)
      continue;
    if (line[0] == '#')
      continue;
    if (sscanf(line, "%d %d %d %f %f %f %f %f %f", &i1, &i2, &i3, &x, &y, &z, &X, &Y, &Z) != 9) {
      nPgons++;
      break;
    }
    if (opgon != i1) { // first-pass initializaton for each solid
      nPgons++;
      opgon = i1;
      pgons.push_back(CpolyhPoints());
      pPg = &pgons.back();
      pPg->whichGe = i1;
      if (i1 > maxPgons)
        maxPgons = i1;
      pPg->npoints = 2 * i2;
      pPg->tubX = -1. * mm;
      pPg->tubY = -1. * mm;
      pPg->tubZ = -1. * mm;
      pPg->tubr = -1. * mm;
      pPg->tubR = -1. * mm;
      pPg->tubL = -1. * mm;
      pPg->capSpace = -1. * mm;
      pPg->capThick = -1. * mm;
      pPg->passThick1 = -1. * mm;
      pPg->passThick2 = -1. * mm;
      pPg->colx = 0.;
      pPg->coly = 0.;
      pPg->colz = 0.;
      pPg->vertex.resize(pPg->npoints);
      pPg->cylinderMakesSense = true;
      pPg->makeCapsule = true;
      pPg->isPlanar = false;
      pPg->segSize_x = -1. * mm;
      pPg->segSize_y = -1. * mm;
      pPg->maxSize_x = -1000. * m;
      pPg->maxSize_y = -1000. * m;
      pPg->minSize_x = 1000. * m;
      pPg->minSize_y = 1000. * m;
      pPg->guardThick[0] = -1. * mm;
      pPg->guardThick[1] = -1. * mm;
      pPg->guardThick[2] = -1. * mm;
      pPg->guardThick[3] = -1. * mm;
      pPg->nSeg_x = 1;
      pPg->nSeg_y = 1;
    }
    if (i2 == 0 && i3 == 0) {
      pPg->tubr = ((G4double)x) * mm;
      pPg->tubR = ((G4double)y) * mm;
      pPg->tubL = ((G4double)z) * mm;
      pPg->tubX = ((G4double)X) * mm;
      pPg->tubY = ((G4double)Y) * mm;
      pPg->tubZ = ((G4double)Z) * mm;
    }
    else if (i2 == 0 && i3 == 1) {
      pPg->thick = ((G4double)x) * mm;
      pPg->passThick1 = ((G4double)y) * mm;
      pPg->passThick2 = ((G4double)z) * mm;
      pPg->capSpace = ((G4double)X) * mm;
      pPg->capThick = ((G4double)Y) * mm;
    }
    else if (i2 == 0 && i3 == 2) {
      pPg->colx = ((G4double)x);
      pPg->coly = ((G4double)y);
      pPg->colz = ((G4double)z);
    }
    else if (i2 == 0 && i3 == 3) { // planar
      pPg->cylinderMakesSense = false;
      pPg->makeCapsule = false;
      pPg->isPlanar = true;
      pPg->nSeg_x = (G4int)x;
      pPg->nSeg_y = (G4int)y;
      pPg->guardThick[0] = ((G4double)z) * mm;
      pPg->guardThick[1] = ((G4double)X) * mm;
      pPg->guardThick[2] = ((G4double)Y) * mm;
      pPg->guardThick[3] = ((G4double)Z) * mm;
    }
    else {
      pPg->vertex[i3] = G4Point3D(((G4double)x), ((G4double)y), ((G4double)z)) * mm;
      pPg->vertex[i3 + i2] = G4Point3D(((G4double)X), ((G4double)Y), ((G4double)Z)) * mm;
      nvdots += 2;
    }
  }

  fclose(fp);
  G4cout << nPgons << " polyhedra for a total of " << nvdots << " vertex points read." << G4endl;

  G4int npt, npg, nn;
  G4double tolerance = 0.5 * mm;
  /*
  #ifdef GRETA
  #ifdef GRETA_DEBUG
    G4bool isOpen = true;
    if( (fp=fopen("asolidG", "w"))==NULL )
      isOpen = false;
  #endif
  #endif
  */
  for (npg = 0; npg < nPgons; npg++) {
    pPg = &pgons[npg];
    npt = pPg->npoints;
    if (!npt)
      continue;

    // calculates z of the two faces of the original polyhedron
    pPg->centerFace1 = G4Point3D();
    pPg->centerFace2 = G4Point3D();
    for (nn = 0; nn < npt / 2; nn++) {
      pPg->centerFace1 += pPg->vertex[nn];
      pPg->centerFace2 += pPg->vertex[nn + npt / 2];
    }
    pPg->centerFace1 /= npt / 2;
    pPg->centerFace2 /= npt / 2;

    pPg->zFace1 = pPg->centerFace1.z();
    pPg->zFace2 = pPg->centerFace2.z();
    pPg->zCenter = 0.5 * (pPg->zFace1 + pPg->zFace2);
    /*
    #ifdef GRETA
    #ifdef GRETA_DEBUG
    // prints out a modified version of the original file, having the first face at zero z coordinate
        if( isOpen ) {
          fprintf( fp, "# solid %d\n", npg);
          for( G4int iii=0; iii<npt/2; iii++ )
            fprintf( fp, "%4.1d %4.1d %4.1d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n",
              pPg->whichGe, npt/2, iii,
              pPg->vertex[iii].x()/mm,
              pPg->vertex[iii].y()/mm,
              pPg->vertex[iii].z()/mm-pPg->zFace1/mm,
              pPg->vertex[iii+npt/2].x()/mm,
              pPg->vertex[iii+npt/2].y()/mm,
              pPg->vertex[iii+npt/2].z()/mm-pPg->zFace1/mm );
          fprintf( fp, "%4.1d %4.1d %4.1d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n",
              pPg->whichGe, 0, 0,
              pPg->tubr/mm, pPg->tubR/mm, pPg->tubL/mm, pPg->tubX/mm, pPg->tubY/mm, pPg->tubZ/mm );
          fprintf( fp, "%4.1d %4.1d %4.1d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n",
              pPg->whichGe, 0, 1,
              pPg->thick/mm, pPg->passThick1/mm, pPg->passThick2/mm, pPg->capSpace/mm, pPg->capThick/mm,
    (pPg->capSpace+pPg->capThick)/mm ); fprintf( fp, "%4.1d %4.1d %4.1d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf
    %12.6lf\n", pPg->whichGe, 0, 2, pPg->colx, pPg->coly, pPg->colz, 0., 0., 0. ); fprintf( fp, "# end solid %d\n#\n",
    npg);
        }
    #endif
    #endif
    */
    // calculates the minimum radius of a cylinder surrounding the polyhedron
    G4double minR2 = 0.;
    G4double theR2;
    for (nn = 0; nn < npt; nn++) {
      theR2 = pow(pPg->vertex[nn].x(), 2.) + pow(pPg->vertex[nn].y(), 2.);
      if (theR2 > minR2)
        minR2 = theR2;
    }
    pPg->minR = sqrt(minR2) + 2. * tolerance; // safety margin!

    // check: to avoid tolerance problems, increase the cylinder length
    if (fabs((pPg->zFace2 - pPg->zFace1) - pPg->tubL) < tolerance)
      pPg->tubL += 2. * tolerance;

    // check the validity of the cylinder
    if (pPg->cylinderMakesSense && (pPg->tubR < 0.)) {
      pPg->cylinderMakesSense = false;
      G4cout << " Warning! Cylinder will not be built in solid " << pPg->whichGe << G4endl;
    }
    if (pPg->cylinderMakesSense && (pPg->tubL < 0.)) {
      pPg->cylinderMakesSense = false;
      G4cout << " Warning! Cylinder will not be built in solid " << pPg->whichGe << G4endl;
    }
    if (pPg->cylinderMakesSense && (pPg->tubr < 0.)) {
      G4cout << " Warning! Setting inner cylinder radius to zero in solid " << pPg->whichGe << G4endl;
      pPg->tubr = 0.;
    }
    if (pPg->cylinderMakesSense && (pPg->tubr > pPg->tubR)) {
      pPg->cylinderMakesSense = false;
      G4cout << " Warning! Cylinder will not be built in solid " << pPg->whichGe << G4endl;
    }
    // if the cylinder does not exceed the polyhedron, keep the cylinder coordinates!
    if (pPg->zCenter - pPg->tubL / 2. > pPg->zFace1)
      pPg->zFace1 = pPg->zCenter - pPg->tubL / 2.;
    if (pPg->zCenter + pPg->tubL / 2. < pPg->zFace2)
      pPg->zFace2 = pPg->zCenter + pPg->tubL / 2.;
    // additional check: if crystal is not long enough, no hole!!!
    if ((pPg->zFace1 + pPg->thick) > pPg->zFace2) {
      G4cout << " Warning! Setting inner cylinder radius to zero in solid " << pPg->whichGe << G4endl;
      pPg->tubr = 0.;
    }

    // passive areas
    // at the back of the detector
    if (pPg->passThick1 < 0.) {
      pPg->passThick1 = 0.;
      G4cout << " Warning! Passive layer will not be built in solid " << pPg->whichGe << G4endl;
    }
    else if (pPg->passThick1 > pPg->tubL) {
      pPg->passThick1 = pPg->tubL;
      G4cout << " Warning! Setting passive layer thickness to " << pPg->tubL / mm << " mm in solid " << pPg->whichGe
             << G4endl;
    }

    // around the coaxial hole
    if (pPg->cylinderMakesSense) {
      if (pPg->passThick2 < 0.) {
        pPg->passThick2 = 0.;
        G4cout << " Warning! Passive layer will not be built in solid " << pPg->whichGe << G4endl;
      }
      else if (pPg->passThick2 > (pPg->tubR - pPg->tubr) || pPg->passThick2 > pPg->thick) {
        pPg->passThick2 = min(pPg->tubR - pPg->tubr, pPg->thick);
        G4cout << " Warning! Setting passive layer thickness to " << (pPg->tubR - pPg->tubr) / mm << " mm in solid "
               << pPg->whichGe << G4endl;
      }
    }

    if (pPg->makeCapsule && (pPg->capSpace <= 0.)) {
      pPg->makeCapsule = false;
      G4cout << " Warning! Capsule will not be built for solid " << pPg->whichGe << G4endl;
    }
    if (pPg->makeCapsule && (pPg->capThick <= 0.)) {
      pPg->makeCapsule = false;
      G4cout << " Warning! Capsule will not be built for solid " << pPg->whichGe << G4endl;
    }

    // planar detectors
    if (pPg->isPlanar) {
      // number of segments
      if (pPg->nSeg_x <= 0) {
        G4cout << " Warning! Invalid segment number for solid " << pPg->whichGe << ", set to 1 instead." << G4endl;
        pPg->nSeg_x = 1;
      }
      if (pPg->nSeg_y <= 0) {
        G4cout << " Warning! Invalid segment number for solid " << pPg->whichGe << ", set to 1 instead." << G4endl;
        pPg->nSeg_y = 1;
      }

      // max, min coordinates
      for (nn = 0; nn < pPg->npoints; nn++) {
        if (pPg->vertex[nn].x() > pPg->maxSize_x)
          pPg->maxSize_x = pPg->vertex[nn].x();
        if (pPg->vertex[nn].y() > pPg->maxSize_y)
          pPg->maxSize_y = pPg->vertex[nn].y();
        if (pPg->vertex[nn].x() < pPg->minSize_x)
          pPg->minSize_x = pPg->vertex[nn].x();
        if (pPg->vertex[nn].y() < pPg->minSize_y)
          pPg->minSize_y = pPg->vertex[nn].y();
      }

      // segment size
      pPg->segSize_x = (pPg->maxSize_x - pPg->minSize_x) / pPg->nSeg_x;
      pPg->segSize_y = (pPg->maxSize_y - pPg->minSize_y) / pPg->nSeg_y;
    }
  }
  /*
  #ifdef GRETA
  #ifdef GRETA_DEBUG
    fclose(fp);
  #endif
  #endif
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Hicari::ReadWallsFile() {
  FILE* fp;
  char line[256];
  G4int lline, i1, i2, i3, i4, i5, i6, nvdots, opgon;
  float x, y, z, X, Y, Z;
  CpolyhPoints* pPg = NULL;

  nWalls = 0;
  nWlTot = 0;
  if (!wallsFile)
    return;

  if ((fp = fopen(wallsFile, "r")) == NULL) {
    G4cout << "\nError opening data file " << wallsFile << G4endl;
    G4cout << "No walls included." << G4endl;
    return;
  }

  G4cout << "\nReading description of walls from file " << wallsFile << " ..." << G4endl;

  walls.clear();

  nvdots = 0;
  opgon = -1;
  /*
  #ifdef GRETA
  #ifdef GRETA_DEBUG
    // writes out a modified copy of the original file
    G4bool isOpen = true;
    FILE *fp1;
    if( (fp1=fopen("awallsG","w"))==NULL )
      isOpen = false;
  #endif
  #endif
  */
  while (fgets(line, 255, fp) != NULL) {
    lline = strlen(line);
    if (lline < 1)
      continue;
    if (line[0] == '#')
      continue;
    if (sscanf(line, "%d %d %d %d %d %d %f %f %f %f %f %f", &i1, &i2, &i3, &i4, &i5, &i6, &x, &y, &z, &X, &Y, &Z) != 12)
      break;
    if (opgon != i3) {
      nWalls++;
      opgon = i3;
      walls.push_back(CpolyhPoints());
      pPg = &walls.back();
      pPg->whichGe = i1;
      pPg->whichCrystal = i2;
      pPg->whichWall = i3;
      pPg->npoints = 2 * i5;
      pPg->vertex.resize(pPg->npoints);
    }
    pPg->vertex[i6] = G4Point3D(((G4double)x), ((G4double)y), ((G4double)z)) * mm;
    pPg->vertex[i6 + i5] = G4Point3D(((G4double)X), ((G4double)Y), ((G4double)Z)) * mm;
    nvdots += 2;
    /*
    #ifdef GRETA
    #ifdef GRETA_DEBUG
        if( isOpen ) {
          CclusterAngles *pCa = NULL;
          CeulerAngles   *pEa = NULL;
          for( G4int iii=0; iii<nClAng; iii++ ) {
            pCa = &clust[iii];
            if( pCa->whichClus == pPg->whichGe ) break;
          }
          for( G4int jjj=0; jjj<pCa->nsolids; jjj++ ) {
            pEa = &(pCa->solids[jjj]);
            if( pEa->numPhys == pPg->whichCrystal ) break;
          }
    //      fprintf( fp1, "%4.1d %4.1d %4.1d %4.1d %4.1d %4.1d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n",
          fprintf( fp1, "%4.1d %4.1d %4.1d %4.1d %4.1d %4.1d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
            pPg->whichGe, pPg->whichCrystal, pPg->whichWall, i4, pPg->npoints/2, i6,
      (*(pEa->pTransf) * G4Point3D( x * mm, y * mm,  z * mm )).x()/mm,
      (*(pEa->pTransf) * G4Point3D( x * mm, y * mm,  z * mm )).y()/mm,
      (*(pEa->pTransf) * G4Point3D( x * mm, y * mm,  z * mm )).z()/mm,
      (*(pEa->pTransf) * G4Point3D( X * mm, Y * mm,  Z * mm )).x()/mm,
      (*(pEa->pTransf) * G4Point3D( X * mm, Y * mm,  Z * mm )).y()/mm,
      (*(pEa->pTransf) * G4Point3D( X * mm, Y * mm,  Z * mm )).z()/mm );
    //	  ((pEa->rotMat) * G4Point3D( ((G4double)x), ((G4double)y),  ((G4double)z) )).x(),
    //	  ((pEa->rotMat) * G4Point3D( ((G4double)x), ((G4double)y),  ((G4double)z) )).y(),
    //	  ((pEa->rotMat) * G4Point3D( ((G4double)x), ((G4double)y),  ((G4double)z) )).z(),
    //	  ((pEa->rotMat) * G4Point3D( ((G4double)X), ((G4double)Y),  ((G4double)Z) )).x(),
    //	  ((pEa->rotMat) * G4Point3D( ((G4double)X), ((G4double)Y),  ((G4double)Z) )).y(),
    //	  ((pEa->rotMat) * G4Point3D( ((G4double)X), ((G4double)Y),  ((G4double)Z) )).z() );
        }
    #endif
    #endif
    */
  }

  fclose(fp);
  /*
  #ifdef GRETA
  #ifdef GRETA_DEBUG
    fclose(fp1);
  #endif
  #endif
  */
  G4cout << nWalls << " polyhedra for a total of " << nvdots << " vertex points read." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Hicari::ReadClustFile() {
  FILE* fp;
  char line[256];
  G4int lline, i1, i2, i3, nsolids, oclust;
  float psi, theta, phi;
  float x, y, z;

  nClAng = 0;
  maxSolids = 0;

  if ((fp = fopen(clustFile, "r")) == NULL) {
    G4cout << "\nError opening data file " << clustFile << G4endl;
    exit(-1);
  }

  G4RotationMatrix rm;

  G4cout << "\nReading description of clusters from file " << clustFile << " ..." << G4endl;

  clust.clear();

  nsolids = 0;
  oclust = -1;
  CclusterAngles* pPg = NULL;
  CeulerAngles* pEa = NULL;
  /*
  #ifdef GRETA
  #ifdef GRETA_DEBUG
    // writes out a modified copy of the original file
    G4bool isOpen = true;
    FILE *fp1;
    if( (fp1=fopen("aclustG","w"))==NULL )
      isOpen = false;
  #endif
  #endif
  */
  while (fgets(line, 255, fp) != NULL) {
    lline = strlen(line);
    if (lline < 2)
      continue;
    if (line[0] == '#')
      continue;
    if (sscanf(line, "%d %d %d %f %f %f %f %f %f", &i1, &i2, &i3, &psi, &theta, &phi, &x, &y, &z) != 9) {
      nClAng++;
      break;
    }
    if (oclust != i1) {
      nClAng++;
      oclust = i1;
      clust.push_back(CclusterAngles());
      pPg = &clust.back();
      pPg->whichClus = i1;
      pPg->nsolids = 0;
      pPg->solids.clear();
      pPg->nwalls = 0;
      pPg->pAssV = new G4AssemblyVolume();
    }
    pPg->solids.push_back(CeulerAngles());
    pEa = &pPg->solids.back();
    pEa->whichGe = i2;
    pEa->numPhys = i3;

    pEa->ps = ((G4double)psi) * deg;
    pEa->th = ((G4double)theta) * deg;
    pEa->ph = ((G4double)phi) * deg;

    pEa->rotMat.set(0, 0, 0);
    pEa->rotMat.rotateZ(((G4double)psi) * deg);
    pEa->rotMat.rotateY(((G4double)theta) * deg);
    pEa->rotMat.rotateZ(((G4double)phi) * deg);

    pEa->trasl = G4ThreeVector(((G4double)x), ((G4double)y), ((G4double)z)) * mm;

    pEa->pTransf = new G4Transform3D(pEa->rotMat, pEa->trasl);

    pPg->nsolids++;
    nsolids++;
    /*
    #ifdef GRETA
    #ifdef GRETA_DEBUG
        if( isOpen ) {
          CpolyhPoints *pPs = NULL;
          for( G4int iii=0; iii<nPgons; iii++ ) {
            pPs = &pgons[iii];
            if( pPs->whichGe == pEa->whichGe ) break;
          }
    //      fprintf( fp1, "%4.1d %4.1d %4.1d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n",
          fprintf( fp1, "%4.1d %4.1d %4.1d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
            pPg->whichClus, pEa->whichGe, pEa->numPhys, pEa->ps/deg, pEa->th/deg, pEa->ph/deg,
            (*(pEa->pTransf) * G4Point3D( 0., 0., pPs->zFace1)).x()/mm,
            (*(pEa->pTransf) * G4Point3D( 0., 0., pPs->zFace1)).y()/mm,
            (*(pEa->pTransf) * G4Point3D( 0., 0., pPs->zFace1)).z()/mm );
        }
    #endif
    #endif
    */
  }

  fclose(fp);
  /*
  #ifdef GRETA
  #ifdef GRETA_DEBUG
    fclose(fp1);
  #endif
  #endif
    */
  for (G4int ii = 0; ii < nClAng; ii++) {
    pPg = &clust[ii];
    if (pPg->nsolids > maxSolids)
      maxSolids = pPg->nsolids;
  }

  G4cout << " Read " << nClAng << " cluster description for a total of " << nsolids << " individual solids." << G4endl;
}

void Hicari::ReadEulerFile() {
  FILE* fp;
  char line[256];
  G4int lline, i1, i2;
  float psi, theta, phi, x, y, z;

  if ((fp = fopen(eulerFile, "r")) == NULL) {
    G4cout << "\nError opening data file " << eulerFile << G4endl;
    exit(-1);
  }

  euler.clear();

  G4cout << "\nReading Euler angles from file " << eulerFile << " ..." << G4endl;
  nEuler = 0;

  G4RotationMatrix rm;
  CeulerAngles* pEa = NULL;

  while (fgets(line, 255, fp) != NULL) {
    lline = strlen(line);
    if (lline < 2)
      continue;
    if (line[0] == '#')
      continue;
    if (sscanf(line, "%d %d %f %f %f %f %f %f", &i1, &i2, &psi, &theta, &phi, &x, &y, &z) != 8)
      break;
    euler.push_back(CeulerAngles());
    pEa = &euler[nEuler];

    pEa->numPhys = i1;
    pEa->whichGe = i2;

    pEa->rotMat.set(0, 0, 0);
    pEa->rotMat.rotateZ(((G4double)psi) * deg);
    pEa->rotMat.rotateY(((G4double)theta) * deg);
    pEa->rotMat.rotateZ(((G4double)phi) * deg);

    pEa->ps = ((G4double)psi) * deg;
    pEa->th = ((G4double)theta) * deg;
    pEa->ph = ((G4double)phi) * deg;

    pEa->trasl = G4ThreeVector(((G4double)x), ((G4double)y), ((G4double)z)) * mm;

    pEa->pTransf = new G4Transform3D(pEa->rotMat, pEa->trasl);

    nEuler++;
  }

  fclose(fp);
  G4cout << nEuler << " Euler angles read." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// This method reads the file with the slice planes
void Hicari::ReadSliceFile() {
  FILE* fp;
  char line[256];
  //  G4int     ns, npts, sameSlice, nSlices; // REMOVE?
  G4int nso, sameSlice, nSlices;
  float zz1, ZZ1;
  G4double zz, ZZ;
  CpolyhPoints* pPg;

  maxSec = 1; // When data file is missing only 1 sector is considered (non-segmented)
  maxSli = 1; // When data file is missing only 1 slice is considered

  if ((fp = fopen(sliceFile, "r")) == NULL) {
    G4cout << "\nError opening data file " << sliceFile << G4endl;
    G4cout << " Warning! No segmentation will be considered." << G4endl;
    readOut = false;
    // When data file is missing only 1 slice is considered
    return;
  }

  G4cout << " Reading slice planes from file " << sliceFile << G4endl;
  readOut = true;

  // Initializes arrays assuming one slice
  // zSliceI[i] is the z coordinate of the i-th plane
  // zSliceI[0] = 1st face,..., zSliceI[nSlice] = 2nd face
  // so that zSliceI.size() = nslice+1
  G4int nPg;
  for (nPg = 0; nPg < nPgons; nPg++) {
    pPg = &pgons[nPg];
    //    npts = pPg->npoints; // REMOVE?
    pPg->zSliceI.clear();
    pPg->zSliceO.clear();
    pPg->zSliceI.push_back(pPg->zFace1);
    pPg->zSliceO.push_back(pPg->zFace1);
    pPg->zSliceI.push_back(pPg->zFace2);
    pPg->zSliceO.push_back(pPg->zFace2);
    pPg->nslice = 1;
  }

  G4cout << "\nReading slicing planes from data file " << sliceFile << G4endl;

  // first line gives the format
  sameSlice = -1;
  while (fgets(line, 255, fp)) {
    if (line[0] == '#')
      continue;
    sscanf(line, "%d", &sameSlice);
    break;
  }
  if (sameSlice < 0) {
    G4cout << "\nError reading slice type in file " << sliceFile << G4endl;
    G4cout << " Warning! Considering only one slice per solid." << G4endl;
    return;
  }

  if (sameSlice == 0) {
    // second line gives number of slices
    nSlices = -1;
    while (fgets(line, 255, fp)) {
      if (line[0] == '#')
        continue;
      sscanf(line, "%d", &nSlices);
      if (nSlices < 1)
        break;
      if (nSlices > maxSli)
        maxSli = nSlices;
      // fills the arrays
      for (nPg = 0; nPg < nPgons; nPg++) {
        pPg = &pgons[nPg];
        zz = pPg->zFace1;
        G4double dz = (pPg->zFace2 - pPg->zFace1) / nSlices;
        pPg->zSliceI.clear();
        pPg->zSliceO.clear();
        for (G4int kk = 0; kk < nSlices; kk++) {
          pPg->zSliceI.push_back(zz + kk * dz);
          pPg->zSliceO.push_back(zz + kk * dz);
        }
        pPg->zSliceI.push_back(pPg->zFace2);
        pPg->zSliceO.push_back(pPg->zFace2);
        pPg->nslice = nSlices;
        if (pPg->npoints / 2 > maxSec)
          maxSec = pPg->npoints / 2;
      }
      break;
    }
    if (nSlices < 1) {
      G4cout << "\nError reading number of slices in file " << sliceFile << G4endl;
      G4cout << " Warning! Considering only one slice per solid." << G4endl;
    }
    fclose(fp);
    return;
  }

  // variable slices
  while (fgets(line, 255, fp)) {
    if (line[0] == '#')
      continue;
    //    if( sscanf(line,"%d %lf %lf", &ns, &zz, &ZZ) == 2 )
    if (sscanf(line, "%d %f %f", &nso, &zz1, &ZZ1) == 2)
      ZZ1 = zz1;
    zz = (G4double)zz1;
    ZZ = (G4double)ZZ1;
    if (nso < 0 || nso > maxPgons) { // needed to avoid problems in case the minimum index in pgons is not zero
      G4cout << " Warning! Solid " << nso << " out of range: ignoring  slice  " << zz << " -- " << ZZ << G4endl;
      continue;
    }
    for (nPg = 0; nPg < nPgons; nPg++) {
      pPg = &pgons[nPg];
      if (pPg->whichGe != nso)
        continue; // looks for the right solid
      G4double z2 = pPg->zSliceI[pPg->nslice - 1];
      G4double Z2 = pPg->zSliceO[pPg->nslice - 1];
      z2 += zz * mm;
      Z2 += ZZ * mm;
      if (z2 > pPg->zFace2 || Z2 > pPg->zFace2) {
        G4cout << " Warning! Slice " << zz << " -- " << ZZ << " of solid " << nso << " exceeds length of solid "
               << "( " << z2 << " -- " << Z2 << "   > " << pPg->zFace2 - pPg->zFace1 << " )" << G4endl;
        continue;
      }
      pPg->zSliceI.back() = z2;
      pPg->zSliceO.back() = Z2;
      pPg->zSliceI.push_back(pPg->zFace2);
      pPg->zSliceO.push_back(pPg->zFace2);
      pPg->nslice++;
      if (pPg->npoints / 2 > maxSec)
        maxSec = pPg->npoints / 2;
      if (pPg->nslice > maxSli)
        maxSli = pPg->nslice;
      break;
    }
  }
  fclose(fp);
  return;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Hicari::ConstructDetector(G4LogicalVolume* world) {

  if (InitializeMaterials())
    return;

  ConstructGeCrystals();

  if (nWalls)
    ConstructTheWalls();

  ConstructTheCapsules();
  /*
    G4int depth;
    if(makeCapsule)
      depth = 2;
    else
      depth = 1;
  */
  // theDetector->GetGammaSD()->SetDepth(depth);

  ConstructTheClusters();
  PlaceTheClusters(world);

  if (readOut) {
    // delete old structures
    nSegments.clear();

    pgSegLl.clear();
    pgSegLu.clear();
    pgSegRl.clear();
    pgSegRu.clear();

    segVolume.clear();
    segCenter.clear();

    nSegments.resize(nPgons);
    tSegments.resize(nPgons);
    totSegments = 0;

    for (G4int ii = 0; ii < nPgons; ii++) {
      G4int nn = CalculateSegments(ii);
      nSegments[ii] = nn;
      tSegments[ii] = totSegments;
      totSegments += nn;
    }
    ConstructSegments();
  }

  WriteCrystalAngles("./crystalangles.dat");
  WriteCrystalPositions("./crystalpositions.dat");
  WriteSegmentPositions("./segmentpositions.dat");
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////// methods to construct and place the actual volumes
/////////////////////////////////////////////////////////////////////////////////////
void Hicari::ConstructGeCrystals() {

  // G4RunManager* runManager = G4RunManager::GetRunManager();

  // DetectorConstruction* theDetector  = (DetectorConstruction*) runManager->GetUserDetectorConstruction();

  if (!matCryst) {
    G4cout << G4endl << "----> Missing material, cannot build the crystals!" << G4endl;
    return;
  }

  char sName[50];

  // identity matrix
  G4RotationMatrix rm;
  rm.set(0, 0, 0);

  G4int ngen, nPg, nGe, nPh = 0, nPt;

  // data to construct the cylinders and the passive parts
  // they are generated in their final position to avoid problems in placement
  G4double* InnRadGe;
  G4double* OutRadGe;
  G4double* zSliceGe;

  G4cout << G4endl << "Generating crystals ... " << G4endl;
  ngen = 0;
  for (nPg = 0; nPg < nPgons; nPg++) {
    CpolyhPoints* pPg = &pgons[nPg];
    nGe = pPg->whichGe;
    if (nGe < 0) {
      G4cout << "ConstructGeCrystals : crystal " << nPg << " skipped because nGe ( " << nGe << " ) out of range "
             << G4endl;
      continue;
    }
    nPt = pPg->npoints;
    if (nPt < 6) {
      G4cout << "ConstructGeCrystals : crystal " << nPg << " skipped because of too few ( " << nPt << " ) points "
             << G4endl;
      continue;
    }

    sprintf(sName, "gePoly%2.2d", nGe);
    pPg->pPoly = new CConvexPolyhedron(G4String(sName), pPg->vertex);

    cout << nGe << ", colx " << pPg->colx << ", coly " << pPg->coly << ", colz " << pPg->colz << endl;

    pPg->pDetVA = new G4VisAttributes(G4Color(pPg->colx, pPg->coly, pPg->colz));

    G4double zFace1 = pPg->zFace1;

    pPg->pDetL1 = NULL;
    pPg->pDetL2 = NULL;

    if (useCylinder && pPg->cylinderMakesSense) {

      // no coaxial hole!
      if (pPg->tubr == 0) {
        zSliceGe = new G4double[2];
        zSliceGe[0] = pPg->zCenter - pPg->tubL / 2.;
        zSliceGe[1] = pPg->zCenter + pPg->tubL / 2.;

        InnRadGe = new G4double[2];
        InnRadGe[0] = pPg->tubr;
        InnRadGe[1] = pPg->tubr;

        OutRadGe = new G4double[2];
        OutRadGe[0] = pPg->tubR;
        OutRadGe[1] = pPg->tubR;

        sprintf(sName, "geTubs%2.2d", nGe);
        pPg->pCoax = new G4Polycone(G4String(sName), 0. * deg, 360. * deg, 2, zSliceGe, InnRadGe, OutRadGe);
        sprintf(sName, "geCapsPolyTubs%2.2d", nGe);
        pPg->pCaps =
            new G4IntersectionSolid(G4String(sName), pPg->pPoly, pPg->pCoax, G4Transform3D(rm, G4ThreeVector()));

        // passive area behind the detector (placed as daughter of the original crystal)
        if (usePassive) {
          if (pPg->passThick1 > 0.) {
            zSliceGe = new G4double[2];
            zSliceGe[0] = pPg->zCenter + pPg->tubL / 2. - pPg->passThick1;
            zSliceGe[1] = pPg->zCenter + pPg->tubL / 2.;

            InnRadGe = new G4double[2];
            InnRadGe[0] = pPg->tubr;
            InnRadGe[1] = pPg->tubr;

            OutRadGe = new G4double[2];
            OutRadGe[0] = pPg->tubR;
            OutRadGe[1] = pPg->tubR;

            sprintf(sName, "geTubsB%2.2d", nGe);
            pPg->pTubs1 = new G4Polycone(G4String(sName), 0. * deg, 360. * deg, 2, zSliceGe, InnRadGe, OutRadGe);
            sprintf(sName, "gePassB%2.2d", nGe);
            pPg->pCaps1 =
                new G4IntersectionSolid(G4String(sName), pPg->pPoly, pPg->pTubs1, G4Transform3D(rm, G4ThreeVector()));
            sprintf(sName, "gePassBL%2.2d", nGe);
            pPg->pDetL1 = new G4LogicalVolume(pPg->pCaps1, matCryst, G4String(sName), 0, 0, 0);
            pPg->pDetL1->SetVisAttributes(pPg->pDetVA);
          }
          else
            pPg->pDetL1 = NULL;
        }
      }
      else {
        zSliceGe = new G4double[4];
        zSliceGe[0] = pPg->zCenter - pPg->tubL / 2.;
        zSliceGe[1] = zFace1 + pPg->thick;
        zSliceGe[2] = zFace1 + pPg->thick + 0.0001;
        zSliceGe[3] = pPg->zCenter + pPg->tubL / 2.;

        InnRadGe = new G4double[4];
        InnRadGe[0] = 0.;
        InnRadGe[1] = 0.;
        InnRadGe[2] = pPg->tubr;
        InnRadGe[3] = pPg->tubr;

        OutRadGe = new G4double[4];
        OutRadGe[0] = pPg->tubR;
        OutRadGe[1] = pPg->tubR;
        OutRadGe[2] = pPg->tubR;
        OutRadGe[3] = pPg->tubR;

        sprintf(sName, "gePcone%2.2d", nGe);
        pPg->pCoax = new G4Polycone(G4String(sName), 0. * deg, 360. * deg, 4, zSliceGe, InnRadGe, OutRadGe);
        sprintf(sName, "geCapsPolyPcone%2.2d", nGe);
        pPg->pCaps =
            new G4IntersectionSolid(G4String(sName), pPg->pPoly, pPg->pCoax, G4Transform3D(rm, G4ThreeVector()));

        if (usePassive) {
          if (pPg->passThick1 > 0.) {
            zSliceGe = new G4double[2];
            zSliceGe[0] = pPg->zCenter + pPg->tubL / 2. - pPg->passThick1;
            zSliceGe[1] = pPg->zCenter + pPg->tubL / 2.;

            InnRadGe = new G4double[2];
            InnRadGe[0] = pPg->tubr;
            InnRadGe[1] = pPg->tubr;

            OutRadGe = new G4double[2];
            OutRadGe[0] = pPg->tubR;
            OutRadGe[1] = pPg->tubR;

            // passive area behind the detector (placed later as daughter of the original crystal)
            sprintf(sName, "geTubsB%2.2d", nGe);
            pPg->pTubs1 = new G4Polycone(G4String(sName), 0. * deg, 360. * deg, 2, zSliceGe, InnRadGe, OutRadGe);
            sprintf(sName, "gePassB%2.2d", nGe);
            pPg->pCaps1 =
                new G4IntersectionSolid(G4String(sName), pPg->pPoly, pPg->pTubs1, G4Transform3D(rm, G4ThreeVector()));

            sprintf(sName, "gePassBL%2.2d", nGe);
            pPg->pDetL1 = new G4LogicalVolume(pPg->pCaps1, matCryst, G4String(sName), 0, 0, 0);
            pPg->pDetL1->SetVisAttributes(pPg->pDetVA);
          }
          else
            pPg->pDetL1 = NULL;

          if (pPg->passThick2 > 0.) {
            // passive area at the coaxial hole (placed later as daughter of the original crystal)
            zSliceGe = new G4double[4];
            zSliceGe[0] = pPg->zFace1 + pPg->thick - pPg->passThick2;
            zSliceGe[1] = pPg->zFace1 + pPg->thick;
            zSliceGe[2] = pPg->zFace1 + pPg->thick;
            zSliceGe[3] = pPg->zCenter + pPg->tubL / 2. - pPg->passThick1;

            InnRadGe = new G4double[4];
            InnRadGe[0] = 0.;
            InnRadGe[1] = 0.;
            InnRadGe[2] = pPg->tubr;
            InnRadGe[3] = pPg->tubr;

            OutRadGe = new G4double[4];
            OutRadGe[0] = pPg->tubr + pPg->passThick2;
            OutRadGe[1] = pPg->tubr + pPg->passThick2;
            OutRadGe[2] = pPg->tubr + pPg->passThick2;
            OutRadGe[3] = pPg->tubr + pPg->passThick2;

            sprintf(sName, "geTubsC%2.2d", nGe);
            pPg->pCoax2 = new G4Polycone(G4String(sName), 0. * deg, 360. * deg, 4, zSliceGe, InnRadGe, OutRadGe);
            sprintf(sName, "gePassC%2.2d", nGe);
            pPg->pDetL2 = new G4LogicalVolume(pPg->pCoax2, matCryst, G4String(sName), 0, 0, 0);
            pPg->pDetL2->SetVisAttributes(pPg->pDetVA);
          }
          else
            pPg->pDetL2 = NULL;
        }
      }
      sprintf(sName, "geDetCapsL%2.2d", nGe);
      pPg->pDetL = new G4LogicalVolume(pPg->pCaps, matCryst, G4String(sName), 0, 0, 0); // intersezione di Poly e Tubs
    }
    else {
      sprintf(sName, "geDetPolyL%2.2d", nGe);
      pPg->pDetL = new G4LogicalVolume(pPg->pPoly, matCryst, G4String(sName), 0, 0, 0); // solo i poliedri

      // passive area behind the detector (placed later as daughter of the original crystal)
      if (usePassive) {
        if (pPg->passThick1 > 0.) {
          zSliceGe = new G4double[2];
          zSliceGe[0] = pPg->zCenter + pPg->tubL / 2. - pPg->passThick1;
          zSliceGe[1] = pPg->zCenter + pPg->tubL / 2.;

          InnRadGe = new G4double[2];
          InnRadGe[0] = 0.;
          InnRadGe[1] = 0.;

          OutRadGe = new G4double[2];
          OutRadGe[0] = pPg->minR;
          OutRadGe[1] = pPg->minR;

          sprintf(sName, "geTubsB%2.2d", nGe);
          pPg->pTubs1 = new G4Polycone(G4String(sName), 0. * deg, 360. * deg, 2, zSliceGe, InnRadGe, OutRadGe);
          sprintf(sName, "gePassB%2.2d", nGe);
          pPg->pCaps1 =
              new G4IntersectionSolid(G4String(sName), pPg->pPoly, pPg->pTubs1, G4Transform3D(rm, G4ThreeVector()));
          sprintf(sName, "gePassBL%2.2d", nGe);
          pPg->pDetL1 = new G4LogicalVolume(pPg->pCaps1, matCryst, G4String(sName), 0, 0, 0);
          pPg->pDetL1->SetVisAttributes(pPg->pDetVA);
        }
        else
          pPg->pDetL1 = NULL;
      }
    }

    if (usePassive) {
      if (pPg->pDetL1) {
        pPg->pDetP1 = NULL;
        sprintf(sName, "gePassBP%3.3d", nPh);
        pPg->pDetP1 = new G4PVPlacement(0, G4ThreeVector(), pPg->pDetL1, G4String(sName), pPg->pDetL, false, 0);
      }
      if (pPg->pDetL2) {
        pPg->pDetP2 = NULL;
        sprintf(sName, "gePassCP%3.3d", nPh);
        pPg->pDetP2 = new G4PVPlacement(0, G4ThreeVector(), pPg->pDetL2, G4String(sName), pPg->pDetL, false, 0);
      }
    }

    pPg->pDetL->SetVisAttributes(pPg->pDetVA);
    pPg->pDetL->SetSensitiveDetector(m_HicariScorer);

    ngen++;

    G4cout << "\n  Total Ge volume (" << pPg->pCaps->GetName() << ")     = " << pPg->pCaps->GetCubicVolume() / cm3
           << " cm3" << G4endl;
    G4cout << "    Back dead layer volume (" << pPg->pCaps1->GetName()
           << ")     = " << pPg->pCaps1->GetCubicVolume() / cm3 << " cm3" << G4endl;
    G4cout << "    Coaxial dead layer volume (" << pPg->pCoax2->GetName()
           << ") = " << pPg->pCoax2->GetCubicVolume() / cm3 << " cm3\n"
           << G4endl;
  }

  G4cout << "Number of generated crystals is " << ngen << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Hicari::ConstructTheCapsules() {
  G4int nPg, nGe, nSid;
  G4bool movePlane;
  G4double dist1 = 0.;
  G4double dist2 = 0.;
  char sName[128];

  CpolyhPoints* pPg = NULL; // germanium
  CpolyhPoints* pPv = NULL; // vacuum
  CpolyhPoints* pPc = NULL; // capsule

  G4RotationMatrix rm;
  rm.set(0, 0, 0);

  if (!matWalls || !matBackWalls || !matHole) {
    G4cout << G4endl << "----> Missing materials, cannot build the capsules!" << G4endl;
    return;
  }

  G4cout << G4endl << "Generating the capsules ... " << G4endl;

  capsI.clear();
  capsO.clear();

  capsI.resize(nPgons);
  capsO.resize(nPgons);

  //  dist2 += dist1;

  for (nPg = 0; nPg < nPgons; nPg++) {
    pPg = &pgons[nPg];

    nGe = pPg->whichGe;

    if (pPg->isPlanar) { // planar: no capsule, guardring
                         //   dist1 = pPg->guardThick;
                         //   dist2 = dist1;

      if (makeCapsule) {
        pPv = &capsO[nPg];
        sprintf(sName, "plaPoly%2.2d", nGe);
        pPv->pPoly = new CConvexPolyhedron(G4String(sName), pPg->vertex);
        //      for( nSid=0; nSid<pPg->pPoly->GetnPlanes()-2; nSid++ )
        //	for( nSid=2; nSid<pPg->pPoly->GetnPlanes(); nSid++ )
        //	  movePlane = pPv->pPoly->MovePlane( nSid, pPg->guardThick[nSid%4] );
        pPv->whichGe = nGe;
        sprintf(sName, "plaPolyL%2.2d", nGe);
        pPv->pDetL = new G4LogicalVolume(pPv->pPoly, matCryst, G4String(sName), 0, 0, 0);
        pPv->pDetVA = new G4VisAttributes(G4Color(pPg->colx, pPg->coly, pPg->colz));
        pPv->pDetVA->SetForceWireframe(true); // KW commented
        pPv->pDetL->SetVisAttributes(pPv->pDetVA);

        pPc = &capsI[nPg];
        sprintf(sName, "plbPoly%2.2d", nGe);
        pPc->pPoly = new CConvexPolyhedron(G4String(sName), pPg->vertex);
        //	for( nSid=0; nSid<pPg->pPoly->GetnPlanes()-2; nSid++ )
        //	for( nSid=2; nSid<pPg->pPoly->GetnPlanes(); nSid++ )
        //	  movePlane = pPc->pPoly->MovePlane( nSid, pPg->guardThick[nSid%4] );
        pPc->whichGe = nGe;
        sprintf(sName, "plbPolyL%2.2d", nGe);
        pPc->pDetL = new G4LogicalVolume(pPc->pPoly, matCryst, G4String(sName), 0, 0, 0);
        pPc->pDetVA = new G4VisAttributes(G4Color(pPg->colx, pPg->coly, pPg->colz));
        pPc->pDetVA->SetForceWireframe(true); // KW commented
        pPc->pDetL->SetVisAttributes(pPc->pDetVA);

        new G4PVPlacement(0, G4ThreeVector(), pPc->pDetL, G4String(sName), pPv->pDetL, false, 0);
        new G4PVPlacement(0, G4ThreeVector(), pPg->pDetL, G4String(sName), pPc->pDetL, false, 0);
      }
      else {
        pPv = &capsO[nPg];
        sprintf(sName, "plaPoly%2.2d", nGe);
        pPv->pPoly = new CConvexPolyhedron(G4String(sName), pPg->vertex);
        //      for( nSid=0; nSid<pPg->pPoly->GetnPlanes()-2; nSid++ )
        //	for( nSid=2; nSid<pPg->pPoly->GetnPlanes(); nSid++ )
        //	  movePlane = pPv->pPoly->MovePlane( nSid, pPg->guardThick[nSid%4] );
        pPv->whichGe = nGe;
        sprintf(sName, "plaPolyL%2.2d", nGe);
        pPv->pDetL = new G4LogicalVolume(pPv->pPoly, matCryst, G4String(sName), 0, 0, 0);
        pPv->pDetVA = new G4VisAttributes(G4Color(pPg->colx, pPg->coly, pPg->colz));
        pPv->pDetVA->SetForceWireframe(true);
        pPv->pDetL->SetVisAttributes(pPv->pDetVA);
        new G4PVPlacement(0, G4ThreeVector(), pPg->pDetL, G4String(sName), pPv->pDetL, false, 0);
      }
    }
    else if (!pPg->makeCapsule) { // no capsule
      if (makeCapsule) {
        pPv = &capsO[nPg];
        sprintf(sName, "cryPoly%2.2d", nGe);
        pPv->pPoly = new CConvexPolyhedron(G4String(sName), pPg->vertex);
        sprintf(sName, "cryCoax%2.2d", nGe);
        pPv->pCoax = new G4Polycone(*(pPg->pCoax));
        sprintf(sName, "cryCaps%2.2d", nGe);
        pPv->pCaps =
            new G4IntersectionSolid(G4String(sName), pPv->pPoly, pPv->pCoax, G4Transform3D(rm, G4ThreeVector()));
        pPv->whichGe = nGe;
        sprintf(sName, "cryPolyL%2.2d", nGe);
        pPv->pDetL = new G4LogicalVolume(pPv->pCaps, matHole, G4String(sName), 0, 0, 0);
        pPv->pDetVA = new G4VisAttributes(G4Color(pPg->colx, pPg->coly, pPg->colz));
        pPv->pDetVA->SetForceWireframe(true);
        pPv->pDetVA->SetVisibility(false);
        pPv->pDetVA->SetDaughtersInvisible(false);
        pPv->pDetL->SetVisAttributes(pPv->pDetVA);

        pPc = &capsI[nPg];
        sprintf(sName, "crzPoly%2.2d", nGe);
        pPc->pPoly = new CConvexPolyhedron(G4String(sName), pPg->vertex);
        sprintf(sName, "crzCoax%2.2d", nGe);
        pPc->pCoax = new G4Polycone(*(pPg->pCoax));
        sprintf(sName, "crzCaps%2.2d", nGe);
        pPc->pCaps =
            new G4IntersectionSolid(G4String(sName), pPc->pPoly, pPc->pCoax, G4Transform3D(rm, G4ThreeVector()));
        pPc->whichGe = nGe;
        sprintf(sName, "crzPolyL%2.2d", nGe);
        pPc->pDetL = new G4LogicalVolume(pPc->pCaps, matHole, G4String(sName), 0, 0, 0);
        pPc->pDetVA = new G4VisAttributes(G4Color(pPg->colx, pPg->coly, pPg->colz));
        pPc->pDetVA->SetForceWireframe(true);
        pPc->pDetVA->SetVisibility(false);
        pPc->pDetVA->SetDaughtersInvisible(false);
        pPc->pDetL->SetVisAttributes(pPc->pDetVA);

        new G4PVPlacement(0, G4ThreeVector(), pPg->pDetL, G4String(sName), pPc->pDetL, false, 0);
        new G4PVPlacement(0, G4ThreeVector(), pPc->pDetL, G4String(sName), pPv->pDetL, false, 0);
      }
      else {
        pPv = &capsO[nPg];
        sprintf(sName, "cryPoly%2.2d", nGe);
        pPv->pPoly = new CConvexPolyhedron(G4String(sName), pPg->vertex);
        sprintf(sName, "cryCoax%2.2d", nGe);
        pPv->pCoax = new G4Polycone(*(pPg->pCoax));
        sprintf(sName, "cryCaps%2.2d", nGe);
        pPv->pCaps =
            new G4IntersectionSolid(G4String(sName), pPv->pPoly, pPv->pCoax, G4Transform3D(rm, G4ThreeVector()));
        pPv->whichGe = nGe;
        sprintf(sName, "cryPolyL%2.2d", nGe);
        pPv->pDetL = new G4LogicalVolume(pPv->pCaps, matHole, G4String(sName), 0, 0, 0);
        pPv->pDetVA = new G4VisAttributes(G4Color(pPg->colx, pPg->coly, pPg->colz));
        pPv->pDetVA->SetForceWireframe(true);
        pPv->pDetVA->SetVisibility(false);
        pPv->pDetVA->SetDaughtersInvisible(false);
        pPv->pDetL->SetVisAttributes(pPv->pDetVA);
        new G4PVPlacement(0, G4ThreeVector(), pPg->pDetL, G4String(sName), pPv->pDetL, false, 0);
      }
    } // else if( !pPg->makeCapsule )
    else {
      dist1 = pPg->capSpace;
      dist2 = dist1 + pPg->capThick;
      if (makeCapsule) {
        // G4cout << " in in " << G4endl;
        // G4cout << " in in " << G4endl;
        // G4cout << " in in " << G4endl;

        // for(int nn=0; nn<pPg->npoints; nn++ ) {
        //   G4cout << pPg->vertex[nn].x() << "\t"  << pPg->vertex[nn].y() << "\t"  << pPg->vertex[nn].z() << G4endl;
        // }

        // vacuum
        pPv = &capsI[nPg];
        pPv->whichGe = nGe;
        sprintf(sName, "vaPoly%2.2d", nGe);
        pPv->pPoly = new CConvexPolyhedron(G4String(sName), pPg->vertex);
        for (nSid = 0; nSid < pPg->pPoly->GetnPlanes(); nSid++)
          movePlane = pPv->pPoly->MovePlane(nSid, dist1);

        sprintf(sName, "vaCoax%2.2d", nGe);
        pPv->pTubs = new G4Tubs(G4String(sName), 0., pPg->tubR + dist1, pPg->tubL, 0. * deg, 360. * deg);

        sprintf(sName, "vaCaps%2.2d", nGe);
        pPv->pCaps =
            new G4IntersectionSolid(G4String(sName), pPv->pPoly, pPv->pTubs, G4Transform3D(rm, G4ThreeVector()));

        sprintf(sName, "geVacPolyL%2.2d", nGe);
        pPv->pDetL = new G4LogicalVolume(pPv->pPoly, matHole, G4String(sName), 0, 0, 0);
        pPv->pDetVA = new G4VisAttributes(G4Color(0, 0, 0));
        pPv->pDetVA->SetForceWireframe(false);
        pPv->pDetVA->SetVisibility(false);
        pPv->pDetL->SetVisAttributes(pPv->pDetVA);

        // actual capsule
        pPc = &capsO[nPg];
        pPc->whichGe = nGe;
        sprintf(sName, "caPoly%2.2d", nGe);
        pPc->pPoly = new CConvexPolyhedron(G4String(sName), pPg->vertex);
        for (nSid = 0; nSid < pPg->pPoly->GetnPlanes(); nSid++)
          movePlane = pPc->pPoly->MovePlane(nSid, dist2);

        sprintf(sName, "caCoax%2.2d", nGe);
        pPc->pTubs = new G4Tubs(G4String(sName), 0., pPg->tubR + dist2, pPg->tubL, 0. * deg, 360. * deg);

        sprintf(sName, "caCaps%2.2d", nGe);
        pPc->pCaps =
            new G4IntersectionSolid(G4String(sName), pPc->pPoly, pPc->pTubs, G4Transform3D(rm, G4ThreeVector()));

        sprintf(sName, "geCapPolyL%2.2d", nGe);
        pPc->pDetL = new G4LogicalVolume(pPc->pCaps, matWalls, G4String(sName), 0, 0, 0);
        pPc->pDetVA = new G4VisAttributes(G4Color(pPg->colx, pPg->coly, pPg->colz));
        // pPc->pDetVA  = new G4VisAttributes( G4Color(1, 0, 0.0) );
        pPc->pDetVA->SetForceWireframe(true);
        pPc->pDetVA->SetVisibility(true);
        pPc->pDetL->SetVisAttributes(pPc->pDetVA);

        new G4PVPlacement(0, G4ThreeVector(), pPg->pDetL, G4String(sName), pPv->pDetL, false, 0);
        new G4PVPlacement(0, G4ThreeVector(), pPv->pDetL, G4String(sName), pPc->pDetL, false, 0);
      }
      else {
        pPv = &capsO[nPg];
        sprintf(sName, "cryPoly%2.2d", nGe);
        pPv->pPoly = new CConvexPolyhedron(G4String(sName), pPg->vertex);
        //	for( nSid=0; nSid<pPg->pPoly->GetnPlanes(); nSid++ )
        //	  movePlane = pPv->pPoly->MovePlane( nSid, dist1 );
        pPv->whichGe = nGe;
        sprintf(sName, "cryPolyL%2.2d", nGe);
        pPv->pDetL = new G4LogicalVolume(pPv->pPoly, matHole, G4String(sName), 0, 0, 0);
        pPv->pDetVA = new G4VisAttributes(G4Color(pPg->colx, pPg->coly, pPg->colz));
        pPv->pDetVA->SetForceWireframe(true);
        pPv->pDetVA->SetVisibility(false);
        pPv->pDetVA->SetDaughtersInvisible(false);
        pPv->pDetL->SetVisAttributes(pPv->pDetVA);
        new G4PVPlacement(0, G4ThreeVector(), pPg->pDetL, G4String(sName), pPv->pDetL, false, 0);
      }
    }
  }
  G4cout << "Number of generated capsules is " << nPgons << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Hicari::ConstructTheClusters() {

  G4int nCa, nPg, nSo;
  CclusterAngles* pCa;
  CpolyhPoints* pPg;
  CeulerAngles* pEa;
  G4RotationMatrix rm;
  G4ThreeVector rotatedPos;
  G4Transform3D transf;

  G4cout << G4endl << "Building the clusters ..." << G4endl;

  for (nCa = 0; nCa < nClAng; nCa++) {
    G4cout << " Cluster #" << nCa << G4endl;
    pCa = &clust[nCa];
    for (nSo = 0; nSo < pCa->nsolids; nSo++) {
      pEa = &pCa->solids[nSo];

      rm = pEa->rotMat;

      rotatedPos = pEa->trasl;

      // germanium detectors
      for (nPg = 0; nPg < nPgons; nPg++) {
        //        if( makeCapsule && pgons[nPg].makeCapsule )
        //          pPg = &capsO[nPg];
        //        else
        //          pPg = &pgons[nPg];
        pPg = &capsO[nPg];

        if (pPg->whichGe != pEa->whichGe)
          continue;
        if (!pPg->pDetL)
          continue;

        if (pCa->pAssV) {
          transf = G4Transform3D(*(pEa->pTransf));
          pCa->pAssV->AddPlacedVolume(pPg->pDetL, transf);
        }

        printf("  Solid %4d      %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n", pPg->whichGe, pEa->ps / deg, pEa->th / deg,
               pEa->ph / deg, rotatedPos.x() / cm, rotatedPos.y() / cm, rotatedPos.z() / cm);
      }
    }
    // the walls
    /*

    #ifdef GRETA
        for( nSo=0; nSo<pCa->nsolids; nSo++ ) {
          pEa = &pCa->solids[nSo];

          rm = pEa->rotMat;

          rotatedPos = pEa->trasl;

          for(nPg=0; nPg<nWalls; nPg++) {
            pPg = &walls[nPg];

            if( pPg->whichGe != pCa->whichClus )
              continue;
    //        if( pPg->whichCrystal != nSo )
            if( pPg->whichCrystal != pEa->numPhys )
              continue;
            if( !pPg->pDetL )
              continue;
            if( pCa->pAssV ) {
              transf = G4Transform3D( *(pEa->pTransf) );
              pCa->pAssV->AddPlacedVolume( pPg->pDetL, transf );
            }
            printf( "   Wall %4d      %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n",
              pPg->whichWall, pEa->ps/deg, pEa->th/deg, pEa->ph/deg, rotatedPos.x()/cm, rotatedPos.y()/cm,
    rotatedPos.z()/cm );
          }
        }
    */
    rm.set(0, 0, 0);
    G4double psi = 0.;
    G4double the = 0.;
    G4double phi = 0.;
    rotatedPos = G4ThreeVector();

    // walls
    for (nPg = 0; nPg < nWalls; nPg++) {
      pPg = &walls[nPg];

      if (pPg->whichGe != pCa->whichClus)
        continue;
      if (!pPg->pDetL)
        continue;

      if (pCa->pAssV) {
        transf = G4Transform3D(rm, rotatedPos);
        pCa->pAssV->AddPlacedVolume(pPg->pDetL, transf);
        pCa->nwalls++;
      }

      printf("  Wall  %4d      %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n", pPg->whichWall, psi / deg, the / deg, phi / deg,
             rotatedPos.x() / cm, rotatedPos.y() / cm, rotatedPos.z() / cm);
    }
  }
}

void Hicari::ConstructTheWalls() {
  char sName[50];
  G4int ngen, nPg, nGe, nPt;

  if (!matWalls || !matBackWalls) {
    G4cout << G4endl << "----> Missing material, cannot build the walls!" << G4endl;
    return;
  }

  G4cout << G4endl << "Generating walls ... " << G4endl;

  ngen = 0;
  for (nPg = 0; nPg < nWalls; nPg++) {
    CpolyhPoints* pPg = &walls[nPg];
    nGe = pPg->whichGe;
    if (nGe < 0 || nGe >= nPgons)
      continue;
    nPt = pPg->npoints;
    if (nPt >= 6) {
      sprintf(sName, "wlPoly%2.2d", nGe);
      pPg->pPoly = new CConvexPolyhedron(G4String(sName), pPg->vertex);

      sprintf(sName, "wlDetL%2.2d", nGe);
      // Use a different material for the back walls                        //LR
      if (nPg < 18) // LR
        pPg->pDetL = new G4LogicalVolume(pPg->pPoly, matWalls, G4String(sName), 0, 0, 0);
      else                                                                                    // LR
        pPg->pDetL = new G4LogicalVolume(pPg->pPoly, matBackWalls, G4String(sName), 0, 0, 0); // LR

      // LR      pPg->pDetVA = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5) );
      // pPg->pDetVA = new G4VisAttributes( G4Colour(1, 1, 1) );
      pPg->pDetVA = new G4VisAttributes(G4Colour(1, 1, 1, 0.3)); // purple
      pPg->pDetVA->SetForceWireframe(false);
      // orig
      pPg->pDetVA->SetVisibility(true);
      // KW walls not visible
      // pPg->pDetVA->SetVisibility(false);
      pPg->pDetL->SetVisAttributes(pPg->pDetVA);

      ngen++;
    }
  }
  G4cout << "Number of generated walls is " << ngen << G4endl;
}

void Hicari::PlaceTheClusters(G4LogicalVolume* world) {
  G4int nGe, nCl, nEa, nCa, nPg, nSol, nPt, indexP;
  G4int ii, jj;

  CclusterAngles* pCa = NULL;
  CpolyhPoints* pPg = NULL;
  CeulerAngles* pEa = NULL;
  CeulerAngles* pEc = NULL;

  G4RotationMatrix rm;
  G4RotationMatrix rm1;
  G4RotationMatrix radd;
  G4RotationMatrix rmP; // PRISMA rotation
  G4ThreeVector rotatedPos;
  G4ThreeVector rotatedPos1;
  G4ThreeVector rotatedPos2;
  G4Transform3D transf;

  G4int iClTot = 0;
  G4int iClMin = -1;
  G4int iClMax = -1;

  nDets = 0;
  nWlTot = 0;
  nClus = 0;

  G4cout << G4endl << "Placing clusters ... " << G4endl;

  // G4RunManager* runManager = G4RunManager::GetRunManager();
  // DetectorConstruction* theDetector  = (DetectorConstruction*) runManager->GetUserDetectorConstruction();

  G4int* iCl = new G4int[nClAng];
  memset(iCl, 0, nClAng * sizeof(G4int));

  arrayRmin = 1.e10;
  arrayRmax = -1.e10;

  crystType.clear();
  crystType.resize(nEuler);

  planarLUT.clear();
  planarLUT.resize(nEuler);

  rmP.set(0, 0, 0);
  if (thetaPrisma != 0.)
    rmP.rotateX(thetaPrisma);

  // For the cryostats
  G4Polycone* Cryostat = new G4Polycone("Cryostat", 0., 360. * deg, 7, cryostatZplanes, cryostatRinner, cryostatRouter);
  G4LogicalVolume* logicCryostat = new G4LogicalVolume(Cryostat, matCryo, "Cryostat_log", 0, 0, 0);

  for (nEa = 0; nEa < nEuler; nEa++) {
    pEc = &euler[nEa];
    nCl = pEc->whichGe;
    if (nCl < 0)
      continue;

    nGe = pEc->numPhys;

    for (nCa = 0; nCa < nClAng; nCa++) {
      pCa = &clust[nCa];
      if (pCa->whichClus != nCl)
        continue;
      if (!pCa->pAssV)
        continue;

      G4cout << nEa << " nEuler = " << nEuler << " pCa whichClus = " << pCa->whichClus << G4endl;

      rm = pEc->rotMat;

      rotatedPos = pEc->trasl + posShift;

      if ((thetaShift * phiShift != 0.) || (thetaShift + phiShift != 0.)) {
        radd.set(0, 0, 0);
        radd.rotateY(thetaShift);
        radd.rotateZ(phiShift);
        rotatedPos = radd(rotatedPos);
        rm = radd * rm;
      }

      if (thetaPrisma != 0.) {
        rotatedPos = rmP(rotatedPos);
        rm = rmP * rm;
      }

      indexP = 1000 * nGe + maxSolids * nGe;
      G4cout << "indexP " << indexP << G4endl;

      transf = G4Transform3D(rm, rotatedPos);
      // pCa->pAssV->MakeImprint(theDetector->HallLog(), transf, indexP-1);
      pCa->pAssV->MakeImprint(world, transf, indexP - 1);

      // Place a Cryostat                                                   //LR
      if (cryostatStatus) {
        cryostatPos = cryostatPos0;
        cryostatPos.rotateZ(pEc->ps);
        cryostatPos.rotateY(pEc->th);
        cryostatPos.rotateZ(pEc->ph);
        cryostatRot = G4RotationMatrix::IDENTITY;
        cryostatRot.rotateY(cryostatPos.getTheta());
        cryostatRot.rotateZ(cryostatPos.getPhi());

        new G4PVPlacement(G4Transform3D(cryostatRot, cryostatPos), logicCryostat, "Cryostat", world, false, 0);
      }

      // Since the solids are defined centered in the origin, need to recalculate
      // the size of the equivalent shell with the crystals placed
      // For this we can neglect the additional rotation radd!
      for (nSol = 0; nSol < pCa->nsolids; nSol++) {
        pEa = &pCa->solids[nSol];

        rm1 = pEa->rotMat;

        rotatedPos1 = pEa->trasl;

        for (nPg = 0; nPg < nPgons; nPg++) {
          pPg = &pgons[nPg];
          if (pPg->whichGe != pEa->whichGe)
            continue;

          for (nPt = 0; nPt < pPg->npoints; nPt++) {
            rotatedPos2 = G4ThreeVector(pPg->vertex[nPt]);
            rotatedPos2 = rm(rm1(rotatedPos2) + rotatedPos1) + rotatedPos;
            arrayRmin = min(arrayRmin, rotatedPos2.mag());
            arrayRmax = max(arrayRmax, rotatedPos2.mag());
          }
          // should consider also the centres of the faces!!!
          rotatedPos2 = G4ThreeVector(pPg->centerFace1);
          rotatedPos2 = rm(rm1(rotatedPos2) + rotatedPos1) + rotatedPos;
          arrayRmin = min(arrayRmin, rotatedPos2.mag());
          arrayRmax = max(arrayRmax, rotatedPos2.mag());
          rotatedPos2 = G4ThreeVector(pPg->centerFace2);
          rotatedPos2 = rm(rm1(rotatedPos2) + rotatedPos1) + rotatedPos;
          arrayRmin = min(arrayRmin, rotatedPos2.mag());
          arrayRmax = max(arrayRmax, rotatedPos2.mag());
        }
      }

      nDets += pCa->nsolids;
      nWlTot += pCa->nwalls;
      nClus++;

      if (iClMin < 0 || nGe < iClMin)
        iClMin = nGe;
      if (iClMax < 0 || nGe > iClMax)
        iClMax = nGe;

      printf("%4d %4d %4d %8d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n", iClTot, nGe, nCl, indexP, pEc->ps / deg,
             pEc->th / deg, pEc->ph / deg, pEc->trasl.x() / cm, pEc->trasl.y() / cm, pEc->trasl.z() / cm);
      iCl[nCl]++;
      iClTot++;
    }
  }
  nClus = iClMax - iClMin + 1;

  // store the crystal type for pulse shape calculations
  crystType.resize((1 + iClMax) * maxSolids);
  planarLUT.resize((1 + iClMax) * maxSolids);
  for (ii = 0; ii < ((G4int)crystType.size()); ii++)
    crystType[ii] = -1;
  for (ii = 0; ii < ((G4int)planarLUT.size()); ii++)
    planarLUT[ii] = 0; // by default, coaxial

  for (nEa = 0; nEa < nEuler; nEa++) {
    pEc = &euler[nEa];
    nCl = pEc->whichGe;
    if (nCl < 0)
      continue;

    nGe = pEc->numPhys;

    for (nCa = 0; nCa < nClAng; nCa++) {
      pCa = &clust[nCa];
      if (pCa->whichClus != nCl)
        continue;
      if (!pCa->pAssV)
        continue;
      for (ii = 0; ii < ((G4int)pCa->solids.size()); ii++) {
        pEa = &pCa->solids[ii];
        crystType[nGe * maxSolids + ii] = pEa->whichGe;

        for (jj = 0; jj < nPgons; jj++) {
          pPg = &pgons[jj];
          if (pPg->whichGe != pEa->whichGe)
            continue;
          if (!pPg->isPlanar)
            continue;
          planarLUT[nGe * maxSolids + ii] = 1;
        }
        // G4cout << nGe*maxSolids+ii << " " << pEa->whichGe << G4endl;
      }
    }
  }

  iCMin = iClMin;
  iCMax = iClMax;
  iGMin = maxSolids * iClMin;
  iGMax = maxSolids * (iClMax + 1) - 1;

  // G4cout << " MaxSolids is " << maxSolids << G4endl;

  G4cout << "Number of placed clusters is " << iClTot << "  [ ";
  for (nPg = 0; nPg < nClAng; nPg++)
    G4cout << iCl[nPg] << " ";
  G4cout << "]" << G4endl;
  G4cout << "Cluster  Index ranging from " << iClMin << " to " << iClMax << G4endl;
  G4cout << "Detector Index ranging from " << maxSolids * iClMin << " to " << maxSolids * (iClMax + 1) - 1 << G4endl
         << G4endl;
  G4cout << "Number of placed walls is " << nWlTot << G4endl << G4endl;
  delete[] iCl;

  G4cout << "The equivalent shell extends from " << arrayRmin / cm << " cm to " << arrayRmax / cm << " cm" << G4endl
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// This method calculates the vertexes of the segments for poly nn
// To allow the thickness at the detector axis to be different from the geometrical (outer) slice
// the segments are decomposed into 4 parts: first they are cut into Left and Right of the reference edge.
// Each half is then split into a lower pyramide (quadrangolar base on the external face, tip at z_inner)
// and an upper tetrahedron (based on the upper pyramide-side-face, tip at Z_inner)

// Warning! these solids are not centered on zero! This avoids a further translation when placing the volumes.

G4int Hicari::CalculateSegments(G4int iPg) {
  G4int sector, slice;

  CpolyhPoints* ppg = &pgons[iPg];

  if (ppg->isPlanar)
    return 0;

  G4int npoints = ppg->npoints;
  G4int nsides = npoints / 2;
  G4int nslices = ppg->nslice;
  G4int nsegs = nsides * nslices;

  // vertices of inner and outer face of crystal
  //#ifdef G4V10
  G4Point3D* vertexF1;
  G4Point3D* vertexF2;
  vertexF1 = new G4Point3D[nsides];
  vertexF2 = new G4Point3D[nsides];
  /*
  #else
    G4Point3DVector vertexF1;
    G4Point3DVector vertexF2;
    vertexF1.resize(nsides);
    vertexF2.resize(nsides);
  #endif
  */
  // for consistency, the points taken directly from CConvexPolyhedron
  for (sector = 0; sector < nsides; sector++) {
    vertexF1[sector] = ppg->pPoly->GetPoints(sector);
    vertexF2[sector] = ppg->pPoly->GetPoints(sector + nsides);
  }
  // center (at cylinder axis) of inner and outer face of crystal
  // tubX, tubY should be zero!
  G4Point3D centerF1(ppg->tubX, ppg->tubY, ppg->zFace1);
  G4Point3D centerF2(ppg->tubX, ppg->tubY, ppg->zFace2);

  G4Plane3D xyPlane;         // a plane normal to zAxis
  G4Plane3D zzPlane;         // a plane passing through the zAxis
  G4Point3D pz, p1, pm, p2;  // the points on the lower segment-face
  G4Point3D PZ, P1, PM, P2;  // the points on the upper segment-face
  CpolyhPoints *ppsl, *ppsu; // pointers to the lower & upper decomposition of half-segment
  nsegs = 0;
  G4int isA, isB;
  for (slice = 0; slice < nslices; slice++) {
    for (sector = 0; sector < nsides; sector++, nsegs++) {
      isA = sector;
      for (int n = 0; n < 2; n++) { // loop on the two faces of the edge
        if (n == 0) {
          isB = (isA + nsides - 1) % nsides; // first towards the previous edge
          pgSegLl.push_back(CpolyhPoints());
          pgSegLu.push_back(CpolyhPoints());
          ppsl = &pgSegLl.back();
          ppsu = &pgSegLu.back();
        }
        else {
          isB = (isA + 1) % nsides; // than towards the next edge
          pgSegRl.push_back(CpolyhPoints());
          pgSegRu.push_back(CpolyhPoints());
          ppsl = &pgSegRl.back();
          ppsu = &pgSegRu.back();
        }
        ppsl->whichGe = iPg;
        ppsu->whichGe = iPg;

        xyPlane = G4Plane3D(0., 0., 1., -ppg->zSliceI[slice]);  // xy-plane at inner lower-level
        pz = XPlaneLine(xyPlane, centerF1, centerF2);           // inner lower-point
        xyPlane = G4Plane3D(0., 0., 1., -ppg->zSliceO[slice]);  // plane at outer lower-level
        p1 = XPlaneLine(xyPlane, vertexF1[isA], vertexF2[isA]); // intercept edge
        p2 = XPlaneLine(xyPlane, vertexF1[isB], vertexF2[isB]); // intercept next/previous edge
        pm = (p1 + p2) / 2;                                     // midpoint at lower-level

        xyPlane = G4Plane3D(0., 0., 1., -ppg->zSliceI[slice + 1]); // plane at inner upper-level
        PZ = XPlaneLine(xyPlane, centerF1, centerF2);              // inner upper-point
        xyPlane = G4Plane3D(0., 0., 1., -ppg->zSliceO[slice + 1]); // plane at outer upper-level
        P1 = XPlaneLine(xyPlane, vertexF1[isA], vertexF2[isA]);    // intercept edge
        P2 = XPlaneLine(xyPlane, vertexF1[isB], vertexF2[isB]);    // intercept next/previous edge
        PM = (P1 + P2) / 2;                                        // midpoint at loupper-level

        // the points of the lower part
        ppsl->vertex.resize(5);
        ppsl->npoints = 5;
        ppsl->vertex[0] = pz;
        ppsl->vertex[1] = p1;
        ppsl->vertex[2] = pm;
        ppsl->vertex[3] = PM;
        ppsl->vertex[4] = P1;

        // description of the pyramide
        ppsl->ifaces.clear();
        ppsl->nfaces = 0;
        ppsl->ifaces.push_back(4); // the quadrangular basis
        ppsl->ifaces.push_back(1);
        ppsl->ifaces.push_back(2);
        ppsl->ifaces.push_back(3);
        ppsl->ifaces.push_back(4);
        ppsl->nfaces++;
        ppsl->ifaces.push_back(3); // the 4 side triangular faces
        ppsl->ifaces.push_back(0);
        ppsl->ifaces.push_back(1);
        ppsl->ifaces.push_back(2);
        ppsl->nfaces++;
        ppsl->ifaces.push_back(3);
        ppsl->ifaces.push_back(0);
        ppsl->ifaces.push_back(2);
        ppsl->ifaces.push_back(3);
        ppsl->nfaces++;
        ppsl->ifaces.push_back(3);
        ppsl->ifaces.push_back(0);
        ppsl->ifaces.push_back(3);
        ppsl->ifaces.push_back(4);
        ppsl->nfaces++;
        ppsl->ifaces.push_back(3);
        ppsl->ifaces.push_back(0);
        ppsl->ifaces.push_back(4);
        ppsl->ifaces.push_back(1);
        ppsl->nfaces++;
        ppsl->ifaces.push_back(-1);

        // the points of the upper part
        ppsu->vertex.resize(4);
        ppsu->npoints = 4;
        ppsu->vertex[0] = pz;
        ppsu->vertex[1] = P1;
        ppsu->vertex[2] = PM;
        ppsu->vertex[3] = PZ;

        // description of the tetrahedron
        ppsu->ifaces.clear();
        ppsu->nfaces = 0;
        ppsu->ifaces.push_back(3); // the "upper" face
        ppsu->ifaces.push_back(1);
        ppsu->ifaces.push_back(2);
        ppsu->ifaces.push_back(3);
        ppsu->nfaces++;
        ppsu->ifaces.push_back(3); // the 3 "side" faces
        ppsu->ifaces.push_back(0);
        ppsu->ifaces.push_back(1);
        ppsu->ifaces.push_back(2);
        ppsu->nfaces++;
        ppsu->ifaces.push_back(3);
        ppsu->ifaces.push_back(0);
        ppsu->ifaces.push_back(2);
        ppsu->ifaces.push_back(3);
        ppsu->nfaces++;
        ppsu->ifaces.push_back(3);
        ppsu->ifaces.push_back(0);
        ppsu->ifaces.push_back(3);
        ppsu->ifaces.push_back(1);
        ppsu->nfaces++;
        ppsu->ifaces.push_back(-1);
      }
    }
  }
  return nsegs;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Hicari::ConstructSegments() {
  char sName1[50], sName2[50];
  G4int nGe;
  G4int iPg, sector, slice;

  G4cout << G4endl << "Generating segments for the ReadOut geometry... " << G4endl;

  G4VisAttributes* segVA[4];
  segVA[0] = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  segVA[1] = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
  segVA[2] = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  segVA[3] = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));

  G4VisAttributes* altVA = new G4VisAttributes(G4Colour(0.1, 0.1, 0.1));
  altVA->SetForceWireframe(true);

  G4int indexS; // index of segment in pgSeg...
  CpolyhPoints* ppgerm;
  CpolyhPoints* ppseg = NULL;
  nSeg = 0;
  G4int* iSeg = new G4int[nPgons];
  memset(iSeg, 0, nPgons * sizeof(G4int));

  for (iPg = 0; iPg < nPgons; iPg++) {
    ppgerm = &pgons[iPg];

    if (ppgerm->isPlanar)
      continue;

    indexS = tSegments[iPg];

    if (printVolumes) {
      G4cout << " Crystal type " << iPg << ": " << G4endl;
      G4cout << "   Segment volumes:" << G4endl;
      G4cout << "                     Total      Coax Passive  Back Passive" << G4endl;
      G4cout << "   Segment           [cm3]         [cm3]         [cm3]" << G4endl;
    }
    G4double segVol, backPassiveVol, coaxPassiveVol;

    for (slice = 0; slice < ppgerm->nslice; slice++) {
      for (sector = 0; sector < ppgerm->npoints / 2; sector++, indexS++) {
        nGe = 100 * iPg + 10 * slice + sector; // --> PPPSs (P=CrystalShape, S =Slice, s=sector)
        segVol = 0.;
        backPassiveVol = 0.;
        coaxPassiveVol = 0.;
        // the four parts composing the segment
        for (int ss = 0; ss < 4; ss++) {
          switch (ss) {
          case 0:
            ppseg = &pgSegLl[indexS]; // the lower segment at the left
            sprintf(sName1, "SegmLl_%5.5d", nGe);
            sprintf(sName2, "SegmLl_L_%5.5d", nGe);
            break;
          case 1:
            ppseg = &pgSegLu[indexS]; // the upper segment at the left
            sprintf(sName1, "SegmLu_%5.5d", nGe);
            sprintf(sName2, "SegmLu_L_%5.5d", nGe);
            break;
          case 2:
            ppseg = &pgSegRl[indexS]; // the lower segment at the right
            sprintf(sName1, "SegmRl_%5.5d", nGe);
            sprintf(sName2, "SegmRl_L_%5.5d", nGe);
            break;
          case 3:
            ppseg = &pgSegRu[indexS]; // the upper segment at the right
            sprintf(sName1, "SegmRu_%5.5d", nGe);
            sprintf(sName2, "SegmRu_L_%5.5d", nGe);
            break;
          }
          ppseg->pPoly = new CConvexPolyhedron(G4String(sName1), ppseg->vertex, ppseg->nfaces, ppseg->ifaces);
          //          ppseg->pDetL  = new G4LogicalVolume( ppseg->pPoly, matCryst, G4String(sName2), 0, 0, 0 );
          // LR: segment volumes were too large (sum > crystal volume). Intersect with crystal.
          G4RotationMatrix rm;
          rm.set(0, 0, 0);
          ppseg->pCaps = new G4IntersectionSolid(G4String(sName1), ppseg->pPoly, ppgerm->pCaps,
                                                 G4Transform3D(rm, G4ThreeVector()));

          // We make these to calculate the dead volume in the segment.
          sprintf(sName1, "SegmPassB_%5.5d", nGe);
          ppseg->pCaps1 = new G4IntersectionSolid(G4String(sName1), ppseg->pPoly, ppgerm->pCaps1,
                                                  G4Transform3D(rm, G4ThreeVector()));
          sprintf(sName1, "SegmPassC_%5.5d", nGe);
          ppseg->pCaps2 = new G4IntersectionSolid(G4String(sName1), ppseg->pPoly, ppgerm->pCoax2,
                                                  G4Transform3D(rm, G4ThreeVector()));

          ppseg->pDetL = new G4LogicalVolume(ppseg->pCaps, matCryst, G4String(sName2), 0, 0, 0);
          ppseg->pDetL->SetVisAttributes(segVA[(indexS) % 4]); // in this way they get also the same color
          if (drawReadOut) {
            new G4PVPlacement(0, G4ThreeVector(), ppseg->pDetL, G4String(sName2), ppgerm->pDetL, false, 0);
            ppgerm->pDetL->SetVisAttributes(altVA);
          }
          // These take a while, so only compute them if we're printing them.
          if (printVolumes) {
            segVol += ppseg->pCaps->GetCubicVolume();
            backPassiveVol += ppseg->pCaps1->GetCubicVolume();
            coaxPassiveVol += ppseg->pCaps2->GetCubicVolume();
          }
          nSeg++;
          iSeg[iPg]++;
        }
        if (printVolumes) {
          G4cout << "   " << std::setw(8) << nGe << std::fixed << std::setprecision(2) << std::setw(14) << segVol / cm3
                 << std::fixed << std::setprecision(2) << std::setw(14) << coaxPassiveVol / cm3 << std::fixed
                 << std::setprecision(2) << std::setw(14) << backPassiveVol / cm3 << G4endl;
        }
      }
    }
    G4cout << " " << iSeg[iPg] << " segments" << G4endl;
  }
  G4cout << " --> Total number of generated sub-segments (4 sub-segments per segment) is " << nSeg << " [ ";
  for (iPg = 0; iPg < nPgons; iPg++)
    G4cout << iSeg[iPg] << " ";
  G4cout << "]" << G4endl;

  G4cout << G4endl << "Checking consistency of segments ..." << G4endl;
  G4int nproblems = 0;
  for (iPg = 0; iPg < nPgons; iPg++) {
    nproblems += CheckOverlap(iPg, tSegments[iPg], nSegments[iPg]);
  }
  if (nproblems)
    G4cout << nproblems << " points with problems" << G4endl;
  else
    G4cout << "all OK" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Hicari::ShowStatus() {
  G4cout << G4endl;
  G4int prec = G4cout.precision(3);
  G4cout.setf(ios::fixed);
  G4cout << " Array composed of " << std::setw(3) << nDets << " detectors" << G4endl;
  G4cout << "       arranged in " << std::setw(3) << nClus << " clusters" << G4endl;
  G4cout << " Array composed of " << std::setw(3) << nWlTot << " walls" << G4endl;
  G4cout << " Description of detectors                read from " << solidFile << G4endl;
  G4cout << " Description of dead materials (walls)   read from " << wallsFile << G4endl;
  G4cout << " Description of clusters                 read from " << clustFile << G4endl;
  G4cout << " Euler angles for clusters               read from " << eulerFile << G4endl;

  if (readOut) {
    G4cout << " Slicing planes for read out of segments read from " << sliceFile << G4endl;
    G4cout << " Generated " << nSeg << " sub-segments " << G4endl;
  }
  if (makeCapsule)
    G4cout << " The capsules have been generated with the proper thickness and spacing." << G4endl;
  else
    G4cout << " The capsules have not been generated." << G4endl;
  if (useCylinder)
    G4cout << " The intersection with a cylinder has been considered in generating the crystals." << G4endl;
  else
    G4cout << " The intersection with a cylinder has not been considered in generating the crystals." << G4endl;
  if (usePassive)
    G4cout << " The passivated zones have been considered." << G4endl;
  else
    G4cout << " The passivated zones have not been considered." << G4endl;

  G4cout << " The detectors   material is " << matCrystName << G4endl;
  G4cout << " The walls       material is " << matWallsName << G4endl;
  if (thetaShift || phiShift)
    G4cout << " The array is rotated by theta, phi = " << thetaShift / deg << ", " << phiShift / deg << " degrees"
           << G4endl;
  if (thetaPrisma != 0.)
    G4cout << " PRISMA rotation is theta = " << thetaPrisma / deg << " degrees" << G4endl;
  if (posShift.mag2())
    G4cout << " The array is shifted by " << posShift / mm << " mm" << G4endl;

  G4cout.unsetf(ios::fixed);
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// plane vv = (a,b,c,d)  pV=(a,b,c)
// line A-->B = AB
// intercept pX = (pV cross (pA cross pB) - d * AB ) / pV dot AB
G4Point3D Hicari::XPlaneLine(const G4Plane3D& vv, const G4Point3D& pA, const G4Point3D& pB) {
  G4Point3D AB;
  G4Normal3D pV;
  G4Point3D AxB, VAB;
  G4double xp;

  pV = vv.normal();

  AB = pB - pA;

  xp = pV.dot(AB);
  if (!xp)
    return G4Point3D();

  AxB = pA.cross(pB);
  VAB = pV.cross(AxB);

  G4Point3D pX;
  pX = (VAB - vv.d() * AB) / xp;

  return pX;
}

G4int Hicari::CheckOverlap(G4int iPg, G4int start, G4int nsegs) {
  CpolyhPoints** ppsegs = new CpolyhPoints*[4 * nsegs]; // collect here the pointers to all segs of this shape

  G4int nstot = 0;
  for (G4int n = 0; n < nsegs; n++) {
    ppsegs[nstot++] = &pgSegLl[start + n];
    ppsegs[nstot++] = &pgSegLu[start + n];
    ppsegs[nstot++] = &pgSegRl[start + n];
    ppsegs[nstot++] = &pgSegRu[start + n];
  }

  CpolyhPoints* ppS1;
  CpolyhPoints* ppS2;
  G4int i1, np1, i2;
  G4Point3D pt1;
  EInside inside;

  G4int nproblems = 0;
  for (i1 = 0; i1 < nstot; i1++) {
    ppS1 = ppsegs[i1];
    for (np1 = 0; np1 < ppS1->npoints; np1++) {
      pt1 = ppS1->pPoly->GetPoints(np1);
      inside = pgons[iPg].pPoly->Inside(G4ThreeVector(pt1));
      if (inside == kOutside) {
        printf("Warning: crystal %d : point %3d of segment %3d(%d) is outside its crystal\n", iPg, np1, i1 / 4, i1 % 4);
        nproblems++;
      }
      for (i2 = 0; i2 < nstot; i2++) {
        if (i2 == i1)
          continue; // no check with itself
        ppS2 = ppsegs[i2];
        inside = ppS2->pPoly->Inside(G4ThreeVector(pt1));
        if (inside == kInside) {
          printf("Warning: crystal %d : point %3d of segmentL %3d(%d) is inside segment %3d(%d)\n", iPg, np1, i1 / 4,
                 i1 % 4, i2 / 4, i2 % 4);
          nproblems++;
        }
      }
    }
  }
  delete[] ppsegs;
  return nproblems;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int Hicari::GetSegmentNumber(G4int nGe, G4ThreeVector position) {
  if (!readOut)
    return 0;
  //  else if( planarLUT[nGe%1000] )
  //    return GetPlanSegmentNumber( nGe, position );
  else
    return GetCoaxSegmentNumber(nGe, position);

  return 0;
}

G4int Hicari::GetCoaxSegmentNumber(G4int nGe, G4ThreeVector position) {
  EInside inside;

  CeulerAngles* pEa = NULL;
  CclusterAngles* pCa = NULL;
  CpolyhPoints* ppgerm = NULL;
  CpolyhPoints* ppseg = NULL;

  G4int detNum = nGe % 1000;
  G4int cluNum = nGe / 1000;
  G4int subIndex = detNum % maxSolids;

  // G4cout << " nGe, det, clu, ind " << nGe << " " << detNum<< " " << cluNum << " " << subIndex << G4endl;

  G4int slice = 0, sector = 0, ss, nCa, nPg, indexS, whichGe;

  G4int whichClus = -100;

  for (nCa = 0; nCa < (G4int)(euler.size()); nCa++) {
    pEa = &euler[nCa];
    // G4cout << nCa << " " << pEa->numPhys << " " << pEa->whichGe << G4endl;
    if (pEa->numPhys == cluNum) {
      whichClus = pEa->whichGe;
      break;
    }
  }

  if (whichClus < 0) {
    G4cout << " Warning! Could not find any detector containing this point: " << position / cm << " cm" << G4endl;
    return -1;
  }

  for (nCa = 0; nCa < nClAng; nCa++) {
    pCa = &clust[nCa];
    if (pCa->whichClus != whichClus)
      continue;

    whichGe = pCa->solids[subIndex].whichGe;
    // G4cout << " whichGe " << whichGe << G4endl;

    for (nPg = 0; nPg < nPgons; nPg++) {
      ppgerm = &pgons[nPg];
      if (ppgerm->whichGe != whichGe)
        continue;
      // G4cout << " ppgerm->whichGe " << nGe << " " << whichGe << " " << ppgerm->whichGe << " " << position/mm <<
      // G4endl;
      indexS = tSegments[nPg];
      for (slice = 0; slice < ppgerm->nslice; slice++) {
        for (sector = 0; sector < ppgerm->npoints / 2; sector++, indexS++) {
          // the four parts composing the segment
          for (ss = 0; ss < 4; ss++) {
            switch (ss) {
            case 0:
              ppseg = &pgSegLl[indexS]; // the lower segment at the left
              break;
            case 1:
              ppseg = &pgSegLu[indexS]; // the upper segment at the left
              break;
            case 2:
              ppseg = &pgSegRl[indexS]; // the lower segment at the right
              break;
            case 3:
              ppseg = &pgSegRu[indexS]; // the upper segment at the right
              break;
            }
            inside = ppseg->pPoly->Inside(position);
            if (inside != kOutside) {
              // G4cout << "slice = " << slice << "   sector = " << sector << G4endl;
              return 10 * slice + sector;
            }
          }
        }
      }

      G4cout << " Warning! Could not find any segment containing this point: " << position / cm << " cm" << G4endl;
      return -1;
    }
  }
  G4cout << " Warning! Could not find any detector containing this point: " << position / cm << " cm" << G4endl;
  return -1;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Hicari::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree* pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("Hicari")) {
    pTree->Branch("Hicari", "THicariData", &m_Event);
  }
  pTree->SetBranchAddress("Hicari", &m_Event);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Hicari::ReadSensitive(const G4Event*) {

  m_Event->Clear();

  ///////////
  GeScorers::PS_GeDetector* Scorer = (GeScorers::PS_GeDetector*)m_HicariScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult();
  // std::cout << "in Hicari, mult"<< size << std::endl;
  for (unsigned int i = 0; i < size; i++) {
    // double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Hicari_NS::ResoEnergy);
    double Energy = Scorer->GetEnergy(i);
    // if(Energy>Hicari_NS::EnergyThreshold){
    double Time = Scorer->GetTime(i);

    int detCode = Scorer->GetIndex(i);
    // detcode = clu_index * 1000 + clu_index * maxSolids (4 in our case)
    // int cluNum = detCode / 1000;
    int detNum = detCode % 1000;
    // int subIndex = detNum % maxSolids;
    int fcluster = detNum / 4;   // same as clNum
    int fcrystalid = detNum % 4; // same as subindex
    /*
          std::cout << "\t Energy"<< Energy << std::endl;
          std::cout << "\t detCode"<< detCode << std::endl;
          std::cout << "\t cluNum"<< cluNum << std::endl;
          std::cout << "\t detNum"<< detNum << std::endl;
          std::cout << "\t subIndex"<< subIndex << std::endl;
          std::cout << "\t fcluster"<< fcluster << std::endl;
          std::cout << "\t fcrystalid"<< fcrystalid << std::endl;
    */

    G4ThreeVector HitPos = Scorer->GetPos(i);
    int segCode = GetSegmentNumber(detCode, HitPos);
    // std::cout << "segCode:" << segCode << std::endl;

    // cout << "read detcode " << detCode<< ", segcode " << segCode << ", detNum " << detNum ;
    //  Modify sector number to match GRETINA data stream
    G4int slice = segCode / 10;
    G4int sector = segCode % 10;
    if (detNum < FIRSTMB) {
      // Type B crystal (offset -1)
      if (sector > 0)
        sector--;
      else
        sector = 5;
      // Type A crystal (offset -2)
      if (detNum % 2) {
        if (sector > 0)
          sector--;
        else
          sector = 5;
      }
    }
    segCode = sector + 6 * slice;

    // cout << "detNum " << detNum << ", sector " << sector << ", slice " << slice << ", segcode " << segCode;
    // mapping MB segments
    if (detNum > FIRSTMB - 1 && detNum < FIRSTMB + 4 * MBCLUST)
      segCode = MBsegments[segCode];
    else if (detNum > FIRSTMB + 4 * MBCLUST - 1)
      segCode = CLsegments[segCode];

    // FOR 8-fold Super CLOVERS this is wrong!!!
    // else if(detNum>FIRSTMB+4*MBCLUST-1)
    //   segCode = sector + 4 * slice;
    // cout << " -> " << segCode << endl;

    // Check if the particle has interact before, if yes, add up the energies.
    unsigned int j = m_Event->Find(fcluster, fcrystalid, segCode);
    if (j != 99999) {
      m_Event->AddE(j, Energy);
    }
    else {
      m_Event->SetCluster(fcluster);
      m_Event->SetCrystal(fcrystalid); // equivalent to subindex
      m_Event->SetSegment(segCode);
      m_Event->SetEnergy(Energy);
      m_Event->SetTime(Time);
    }
    /*
        std::cout <<"mult:"<< m_Event->GetMult() << std::endl;
        //Blur the enrgies using the resolution
        for(unsigned int i =0; i<m_Event->GetMult(); i++){
          std::cout <<"before:"<< m_Event->GetEnergy(i) << std::endl;
          double tempE = RandGauss::shoot(m_Event->GetEnergy(i),Hicari_NS::ResoEnergy);
          std::cout<<"after:" << tempE << std::endl;
          m_Event->ChangeE(i, tempE);
        }
    */
  }

  /*
      //G4int Nhits = gammaCollection->entries();

      unsigned int Nhits = Scorer->GetMult();
      std::cout << "Nhits, number of interaction points in Hicari: "<< size << std::endl;

      if(Nhits>0) {

          // Packing: consolidate interaction points within segments
          // based on proximity.
          int MAX_SEGS=8;
          int MAX_INTPTS=16;
          int MBCLUST=6;
          int CLOVERS=4;
          double packingRes=0.*mm;
          double hitTolerance = 0.00001*mm;
  //#define FIRSTMB 800
  //#define MBSEGS 6
  //#define CLSEGS 4

          G4int trackID[100*MAX_INTPTS];
          G4int detNum[100*MAX_INTPTS];
          G4int segNum[100*MAX_INTPTS];
          G4double measuredEdep[100*MAX_INTPTS];
          G4double segmentEdep[100*MAX_INTPTS];
          G4double measuredX[100*MAX_INTPTS];
          G4double measuredY[100*MAX_INTPTS];
          G4double measuredZ[100*MAX_INTPTS];
          G4double X0[100*MAX_INTPTS];
          G4double Y0[100*MAX_INTPTS];
          G4double Z0[100*MAX_INTPTS];
          G4int NCons[100*MAX_INTPTS];
          G4double packingRes2 = packingRes*packingRes;


          G4int NMeasured = 0;
          G4double totalEdep = 0;

          G4bool foundfirst[(MBCLUST+CLOVERS)*MAXCRYSTALNO][MAXSEGNO];
          for(int cl=0; cl<MBCLUST+CLOVERS; cl++){
              for(int cr=0; cr<MAXCRYSTALNO; cr++){
                  for(int se=0; se<MAXSEGNO; se++){
                      foundfirst[cl*MAXCRYSTALNO+cr][se] = false;
                  }
              }
          }


          for(G4int i = 0; i < Nhits; i++){

              G4double x, y, z;
              x = Scorer->GetPos(i).getX()/mm;
              y = Scorer->GetPos(i).getY()/mm;
              z = Scorer->GetPos(i).getZ()/mm;

              G4double en  = Scorer->GetEnery(i)/keV;
              totalEdep += en;

              NCons[i] = -1;
              G4bool processed = false;

              // Initialize a new interaction point for each gamma-ray hit.
              if((*gammaCollection)[i]->GetParticleID() == "gamma"){

                  // Combine multiple gamma hits at the same position.
                  // (This is rare, but it happens.)
                  if(i > 0
                          && (x - measuredX[i-1])*(x - measuredX[i-1]) < hitTolerance*hitTolerance
                          && (y - measuredY[i-1])*(y - measuredY[i-1]) < hitTolerance*hitTolerance
                          && (z - measuredZ[i-1])*(z - measuredZ[i-1]) < hitTolerance*hitTolerance){

                      measuredEdep[NMeasured-1] += en;
                      processed = true;

                  } else {

                      trackID[NMeasured]      = (*gammaCollection)[i]->GetTrackID();
                      detNum[NMeasured]       = (*gammaCollection)[i]->GetDetNumb();
                      segNum[NMeasured]       = (*gammaCollection)[i]->GetSegNumb();
                      //G4cout << trackID[NMeasured] << "\t" << detNum[NMeasured] << "\t" << segNum[NMeasured] <<
  G4endl;

                      // This becomes the total energy deposit associated with this
                      // interaction.
                      measuredEdep[NMeasured] = en;

                      // This becomes the barycenter of all energy depositions associated
                      // with this interaction.
                      measuredX[NMeasured]    = x;
                      measuredY[NMeasured]    = y;
                      measuredZ[NMeasured]    = z;
                      trackID[NMeasured]      = (*gammaCollection)[i]->GetTrackID();

                      // Position of the initial interaction. We use position to identify
                      // the tracks produced by this interaction.
                      X0[NMeasured] = (*gammaCollection)[i]->GetPos().getX()/mm;
                      Y0[NMeasured] = (*gammaCollection)[i]->GetPos().getY()/mm;
                      Z0[NMeasured] = (*gammaCollection)[i]->GetPos().getZ()/mm;

                      //find first interaction for HiCARI
                      //cout << "interaction in det = " << detNum[NMeasured] << ", seg = "<< segNum[NMeasured] << endl;
                      if(detNum[NMeasured]>=FIRSTMB){
                          if(!foundfirst[detNum[NMeasured]-FIRSTMB][segNum[NMeasured]]){
                              //		  cout << "first interaction in det = " << detNum[NMeasured] << " - FIRSTMB "
  <<detNum[NMeasured]-FIRSTMB << ", seg = "<< segNum[NMeasured] << " at (x,y,z) = (" << x <<" ," << y <<" ," << z <<")"
  << endl; foundfirst[detNum[NMeasured]-FIRSTMB][segNum[NMeasured]] = true; G4ThreeVector gamma(x,y,z);
                              firstposition[detNum[NMeasured]-FIRSTMB][segNum[NMeasured]]+=gamma;
                              firstctr[detNum[NMeasured]-FIRSTMB][segNum[NMeasured]]++;
                          }
                      }

                      NCons[NMeasured] = 1;
                      NMeasured++;
                      processed = true;
                  }

                  // Combine secondary-particle hits with their parent interaction points.
              } else {

                  // Compare hit i with existing "measured" interaction points.
                  for(G4int j = 0; j < NMeasured; j++){

                      G4double x0 = (*gammaCollection)[i]->GetTrackOrigin().getX()/mm;
                      G4double y0 = (*gammaCollection)[i]->GetTrackOrigin().getY()/mm;
                      G4double z0 = (*gammaCollection)[i]->GetTrackOrigin().getZ()/mm;

                      // G4cout << "(*gammaCollection)["<< i <<"]->GetParentTrackID() = "
                      // 	   << (*gammaCollection)[i]->GetParentTrackID()
                      // 	   << "   trackID[j] = " << trackID[j]
                      // 	   << "   (x0 - X0[" << j << "]) = " << (x0 - X0[j])
                      // 	   << "   (y0 - Y0[" << j << "]) = " << (y0 - Y0[j])
                      // 	   << "   (z0 - Z0[" << j << "]) = " << (z0 - Z0[j])
                      // 	   << "   (*gammaCollection)[" << i << "]->GetDetNumb()"
                      // 	   << (*gammaCollection)[i]->GetDetNumb()
                      // 	   << "   detNum[" << j << "] = " << detNum[j] << G4endl;

                      if( (*gammaCollection)[i]->GetParentTrackID() == trackID[j]  // correct parent
                              && (x0 - X0[j])*(x0 - X0[j]) < hitTolerance*hitTolerance
                              && (y0 - Y0[j])*(y0 - Y0[j]) < hitTolerance*hitTolerance
                              && (z0 - Z0[j])*(z0 - Z0[j]) < hitTolerance*hitTolerance // correct interaction point
                              && (*gammaCollection)[i]->GetDetNumb() == detNum[j]){        // same crystal

                          // Energy-weighted average position (barycenter)
                          measuredX[j] = (measuredEdep[j]*measuredX[j] + en*x)/(measuredEdep[j] + en);
                          measuredY[j] = (measuredEdep[j]*measuredY[j] + en*y)/(measuredEdep[j] + en);
                          measuredZ[j] = (measuredEdep[j]*measuredZ[j] + en*z)/(measuredEdep[j] + en);
                          measuredEdep[j] += en;

                          NCons[j]++;
                          processed = true;

                      }
                  }
              }

              // If hit i is not a gamma-ray hit and cannot be consolidated with
              // an existing gamma-ray interaction point, it's a positron or an
              // electron multiple-scattering event that isn't associated with
              // energy deposition in the sensitive volume of GRETINA by a gamma
              // ray. (The gamma-ray interaction happened in dead material.)
              // We'll initialize a new interaction point and treat it as a
              // gamma-ray interaction.
              if(!processed){

                  trackID[NMeasured]      = (*gammaCollection)[i]->GetTrackID();
                  detNum[NMeasured]       = (*gammaCollection)[i]->GetDetNumb();
                  segNum[NMeasured]       = (*gammaCollection)[i]->GetSegNumb();
                  measuredEdep[NMeasured] = en;
                  measuredX[NMeasured]    = x;
                  measuredY[NMeasured]    = y;
                  measuredZ[NMeasured]    = z;

                  // This is not a gamma ray. We need to trick its siblings into
                  // treating it as the parent gamma.
                  trackID[NMeasured]      = (*gammaCollection)[i]->GetParentTrackID();
                  X0[NMeasured] = (*gammaCollection)[i]->GetTrackOrigin().getX()/mm;
                  Y0[NMeasured] = (*gammaCollection)[i]->GetTrackOrigin().getY()/mm;
                  Z0[NMeasured] = (*gammaCollection)[i]->GetTrackOrigin().getZ()/mm;

                  NCons[NMeasured] = 1;
                  NMeasured++;
                  processed = true;

              }

              if(!processed)
                  G4cout << "Warning: Could not find a home for hit " << i
                      << " of event " << event_id
                      << G4endl;

              if(NMeasured >= 100*MAX_INTPTS){
                  G4cout << "Error: too many decomposed hits. Increase hit processing array dimension."
                      << G4endl;
                  exit(EXIT_FAILURE);
              }

          }

          // Packing: Consolidate the "measured" gamma-ray interaction points
          // within a single segment that are closer than the PackingRes
          // parameter.
          G4int NGammaHits = NMeasured;
          for(G4int i = 0; i < NMeasured; i++){
              for(G4int j = i+1; j < NMeasured; j++){

                  if( ( (measuredX[i] - measuredX[j])*(measuredX[i] - measuredX[j])
                              + (measuredY[i] - measuredY[j])*(measuredY[i] - measuredY[j])
                              + (measuredZ[i] - measuredZ[j])*(measuredZ[i] - measuredZ[j])
                              < packingRes2 )                                  // proximal
                          && detNum[i] == detNum[j]                          // same crystal
                          && segNum[i] == segNum[j]                          // same segment
                          &&  (NCons[i] > 0 && NCons[j] > 0) ){  // not already consolidated

                      // Energy-weighted average
                      measuredX[i] = (measuredEdep[i]*measuredX[i] +
  measuredEdep[j]*measuredX[j])/(measuredEdep[i]+measuredEdep[j]); measuredY[i] = (measuredEdep[i]*measuredY[i] +
  measuredEdep[j]*measuredY[j])/(measuredEdep[i]+measuredEdep[j]); measuredZ[i] = (measuredEdep[i]*measuredZ[i] +
  measuredEdep[j]*measuredZ[j])/(measuredEdep[i]+measuredEdep[j]); measuredEdep[i] += measuredEdep[j];

                      NCons[j] = -1;
                      NGammaHits--;

                  }

              }

          }

          // Calculate the total energy deposited in each segment.
          for(G4int i = 0; i < NMeasured; i++) // initialize
              if(NCons[i] > 0)
                  segmentEdep[i] = measuredEdep[i];

          G4bool singleDetector = true;
          for(G4int i = 0; i < NMeasured; i++){
              for(G4int j = i+1; j < NMeasured; j++){
                  if(NCons[i] > 0 && NCons[j] > 0
                          && detNum[i] == detNum[j]    // same crystal
                          && segNum[i] == segNum[j]){ // same segment
                      segmentEdep[i] += measuredEdep[j];
                      segmentEdep[j] += measuredEdep[i];
                  }
                  if( detNum[i] != detNum[j] )
                      singleDetector = false;
              }
          }

          // Identify events in which the full emitted gamma-ray energy
          // is deposited in a single crystal
          // (only evaluated for emitted multiplicity = 1 events).
          if( eventInfo->GetNEmittedGammas() == 1 ){
              if( singleDetector &&
                      (totalEdep - eventInfo->GetEmittedGammaEnergy(0))
                      *(totalEdep - eventInfo->GetEmittedGammaEnergy(0))
                      < 0.001*keV*0.001*keV )
                  eventInfo->SetFullEnergy(1);
              else
                  eventInfo->SetFullEnergy(0);
          }


          // Write event to the output file
          if(fisInBeam){
              //write out zero degree "data"
              writeZeroDeg(timestamp,eventInfo->GetATA(), eventInfo->GetBTA(), eventInfo->GetXTA(), eventInfo->GetYTA(),
  eventInfo->GetBetaTA()); if(fisMINOS){ writeMINOS(timestamp,eventInfo->GetVertex(), eventInfo->GetBetaRE());
              }
          }

          // Fold position resolution measuredX, measuredY, measuredZ
          // WARNING: this "smears" positions across the boundaries of
          //          the active volume, which doesn't correspond to
          //          the behavior of signal decomposition.
          if(posRes > 0){
              for(int i=0; i<NMeasured; i++) {
                  measuredX[i] += CLHEP::RandGauss::shoot(0, posRes);
                  measuredY[i] += CLHEP::RandGauss::shoot(0, posRes);
                  measuredZ[i] += CLHEP::RandGauss::shoot(0, posRes);
              }
          }
          //NMeasured, detNum, segNum

          if(print){
              // Write decomposed gamma event(s) to the output file
              cout << "writing NMeasured "<<NMeasured << endl;
              for(int nn =0;nn<NMeasured;nn++){
                  cout << "detNum: " << detNum[nn] << ", segNum: " << segNum[nn] << ", Ncons " << NCons[nn];
                  //for(int jj=0;jj>NCons[nn];jj++)
                  cout << ", edep " << measuredEdep[nn] <<", segment "<<  segmentEdep[nn] << endl;
              }
          }
          writeDecomp(timestamp,
                  NMeasured,
                  detNum, segNum, NCons,
                  measuredX, measuredY, measuredZ,
                  measuredEdep, segmentEdep);

      }//loop hits

  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void Hicari::WriteCrystalAngles(G4String file) {
  G4int nGe, nPh, nCl;
  G4int iPg, iCa, ne, nSol;

  G4double theta, phi;

  CeulerAngles* peA;
  CeulerAngles* peA1;
  CclusterAngles* pcA;
  CpolyhPoints* ppG;
  char line[128];

  G4RotationMatrix rm, radd, rm1, rmP, frameRot;
  G4ThreeVector trasl, trasl1;
  G4ThreeVector rotatedPos;

  std::ofstream outFileLMD;

  outFileLMD.open(file);
  if (!outFileLMD.is_open()) {
    G4cout << " --> Could not open " << file << " output file, aborting ..." << G4endl;
    return;
  }

  G4cout << " --> Writing out crystal angles to " << file << " file..." << G4endl;
  for (ne = 0; ne < nEuler; ne++) {

    peA = &euler[ne];
    nCl = peA->whichGe;
    if (nCl < 0)
      continue;
    nPh = peA->numPhys * maxSolids;

    rm = peA->rotMat;
    trasl = peA->trasl + posShift;

    if ((thetaShift * phiShift != 0.) || (thetaShift + phiShift != 0.)) {
      radd.set(0, 0, 0);
      radd.rotateY(thetaShift);
      radd.rotateZ(phiShift);
      trasl = radd(trasl);
      rm = radd * rm;
    }

    if (thetaPrisma != 0.) {
      rmP.set(0, 0, 0);
      rmP.rotateX(thetaPrisma);
      rm = rmP * rm;
      trasl = rmP(trasl);
    }

    for (iCa = 0; iCa < nClAng; iCa++) {
      pcA = &clust[iCa];
      if (pcA->whichClus != nCl)
        continue;

      for (nSol = 0; nSol < pcA->nsolids; nSol++) {
        peA1 = &pcA->solids[nSol];
        nGe = peA1->whichGe;
        if (nGe < 0)
          continue;

        rm1 = peA1->rotMat;

        frameRot = rm * rm1;

        trasl1 = peA1->trasl;

        for (iPg = 0; iPg < nPgons; iPg++) {
          ppG = &pgons[iPg];
          if (ppG->whichGe != nGe)
            continue; // looks for the right solid
          rotatedPos = rm(trasl1) + trasl;

          theta = rotatedPos.theta();
          phi = rotatedPos.phi();
          if (phi < 0.)
            phi += 360. * deg;

          sprintf(line, "  Riv#%4.1d    Theta= %9.4f        Phi= %9.4f\n", nPh, theta / deg, phi / deg);
          // sprintf( line, " %3d %9d %9.4f %9.4f\n", 0, nPh, theta/deg, phi/deg );
          outFileLMD << line;
          G4cout << line;
          nPh++;
        }
      }
    }
  }
  outFileLMD.close();
  G4cout << " --> Crystal angles successfully written out to " << file << " file." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void Hicari::WriteCrystalPositions(std::ofstream &outFileLMD, G4double unit)
void Hicari::WriteCrystalPositions(G4String file, G4double unit) {
  G4int nGe, nPh, nCl;
  G4int iPg, iCa, ne, nSol;

  CeulerAngles* peA;
  CeulerAngles* peA1;
  CclusterAngles* pcA;
  CpolyhPoints* ppG;
  char line[128];

  G4RotationMatrix rm, radd, rm1, rmP, frameRot;
  G4ThreeVector trasl, trasl1;
  G4ThreeVector rotatedPos;

  std::ofstream outFileLMD;
  outFileLMD.open(file);
  if (!outFileLMD.is_open()) {
    G4cout << " --> Could not open " << file << " output file, aborting ..." << G4endl;
    return;
  }
  G4cout << " --> Writing out crystal positiond to " << file << " file..." << G4endl;

  for (ne = 0; ne < nEuler; ne++) {

    peA = &euler[ne];
    nCl = peA->whichGe;
    if (nCl < 0)
      continue;
    nPh = peA->numPhys * maxSolids;

    rm = peA->rotMat;
    trasl = peA->trasl + posShift;

    if ((thetaShift * phiShift != 0.) || (thetaShift + phiShift != 0.)) {
      radd.set(0, 0, 0);
      radd.rotateY(thetaShift);
      radd.rotateZ(phiShift);
      trasl = radd(trasl);
      rm = radd * rm;
    }

    if (thetaPrisma != 0.) {
      rmP.set(0, 0, 0);
      rmP.rotateX(thetaPrisma);
      rm = rmP * rm;
      trasl = rmP(trasl);
    }

    for (iCa = 0; iCa < nClAng; iCa++) {
      pcA = &clust[iCa];
      if (pcA->whichClus != nCl)
        continue;

      for (nSol = 0; nSol < pcA->nsolids; nSol++) {
        peA1 = &pcA->solids[nSol];
        nGe = peA1->whichGe;
        if (nGe < 0)
          continue;

        rm1 = peA1->rotMat;

        frameRot = rm * rm1;

        trasl1 = peA1->trasl;

        for (iPg = 0; iPg < nPgons; iPg++) {
          ppG = &pgons[iPg];
          if (ppG->whichGe != nGe)
            continue; // looks for the right solid
          rotatedPos = rm(trasl1) + trasl;

          sprintf(line, " %3d  0    %10.5f %10.5f %10.5f\n", nPh, rotatedPos.x() / unit, rotatedPos.y() / unit,
                  rotatedPos.z() / unit);
          outFileLMD << line;
          sprintf(line, "      1    %10.5f %10.5f %10.5f\n", frameRot.xx(), frameRot.xy(), frameRot.xz());
          outFileLMD << line;
          sprintf(line, "      2    %10.5f %10.5f %10.5f\n", frameRot.yx(), frameRot.yy(), frameRot.yz());
          outFileLMD << line;
          sprintf(line, "      3    %10.5f %10.5f %10.5f\n", frameRot.zx(), frameRot.zy(), frameRot.zz());
          outFileLMD << line;
          nPh++;
        }
      }
    }
  }

  outFileLMD.close();
  G4cout << " --> Crystal positions successfully written out to " << file << " file." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// void Hicari::WriteSegmentPositions(std::ofstream &outFileLMD, G4double unit)
void Hicari::WriteSegmentPositions(G4String file, G4double unit) {
  G4int nGe, nPh, nCl;
  G4int iPg, iCa, ne, nSol, sector, slice;
  G4double step = 1. * mm / ((G4double)stepFactor); // integration step for CalculateVolumeAndCenter()

  if ((segVolume.size() == 0) || stepHasChanged) {
    G4cout << "Calculating volume and center of segments (with step = " << step / mm << "mm ) ..." << G4endl;
    segVolume.resize(totSegments);
    segCenter.resize(totSegments);
    for (iPg = 0; iPg < nPgons; iPg++) {
      CalculateVolumeAndCenter(iPg, tSegments[iPg], nSegments[iPg], step);
    }
    stepHasChanged = false;
  }

  CeulerAngles* peA;
  CeulerAngles* peA1;
  CclusterAngles* pcA;
  CpolyhPoints* ppG;
  G4int indexS; // index of segment in segCenter
  char line[128];

  G4RotationMatrix rm, rm1, radd, rmP;
  G4ThreeVector trasl, trasl1;
  G4Point3D centreSeg, rotatedPos;
  G4Transform3D clusterToWorld, crystalToCluster, crystalToWorld;

  std::ofstream outFileLMD;
  outFileLMD.open(file);
  if (!outFileLMD.is_open()) {
    G4cout << " --> Could not open " << file << " output file, aborting ..." << G4endl;
    return;
  }
  G4cout << " --> Writing out segment positions to " << file << " file..." << G4endl;

  for (ne = 0; ne < nEuler; ne++) {

    peA = &euler[ne];
    nCl = peA->whichGe;
    if (nCl < 0)
      continue;
    nPh = peA->numPhys * maxSolids;

    trasl = peA->trasl + posShift;
    rm = peA->rotMat;

    if ((thetaShift * phiShift != 0.) || (thetaShift + phiShift != 0.)) {
      radd.set(0, 0, 0);
      radd.rotateY(thetaShift);
      radd.rotateZ(phiShift);
      trasl = radd(trasl);
      rm = radd * rm;
    }

    if (thetaPrisma != 0.) {
      rmP.set(0, 0, 0);
      rmP.rotateX(thetaPrisma);
      rm = rmP * rm;
      trasl = rmP(trasl);
    }

    clusterToWorld = G4Transform3D(rm, trasl);

    for (iCa = 0; iCa < nClAng; iCa++) {
      pcA = &clust[iCa];
      if (pcA->whichClus != nCl)
        continue;

      for (nSol = 0; nSol < pcA->nsolids; nSol++) {
        peA1 = &pcA->solids[nSol];
        nGe = peA1->whichGe;
        if (nGe < 0)
          continue;

        rm1 = peA1->rotMat;

        trasl1 = peA1->trasl;

        crystalToCluster = G4Transform3D(rm1, trasl1);

        crystalToWorld = clusterToWorld * crystalToCluster;

        for (iPg = 0; iPg < nPgons; iPg++) {
          ppG = &pgons[iPg];
          if (ppG->whichGe != nGe)
            continue; // looks for the right solid
          indexS = tSegments[iPg];

          for (slice = 0; slice < ppG->nslice; slice++) {
            for (sector = 0; sector < ppG->npoints / 2; sector++, indexS++) {

              centreSeg = G4Point3D(segCenter[indexS]);

              // rotatedPos  = rm( rm1( centreSeg ) + trasl1 ) + trasl; // old style!
              //  could be: rotatedPos = clusterToWorld * (crystalToCluster*centreSeg)
              rotatedPos = crystalToWorld * centreSeg;

              sprintf(line, " %3d %2d %2d %10.5f %10.5f %10.5f  %10.5f\n", nPh, slice, sector, rotatedPos.x() / unit,
                      rotatedPos.y() / unit, rotatedPos.z() / unit, segVolume[indexS] / unit * unit * unit);
              outFileLMD << line;
            }
          }

          nPh++;
        }
      }
    }
  }

  outFileLMD.close();
  G4cout << " --> Segment positions successfully written out to " << file << " file." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void Hicari::CalculateVolumeAndCenter(G4int iPg, G4int start, G4int nsegs, G4double step) {
  G4int nse, ss, npts, nn;
  G4double xMin, yMin, zMin;
  G4double xMax, yMax, zMax;
  std::vector<CpolyhPoints*> ppsegs; // collect here the pointers the four sub-segments
  ppsegs.resize(4);
  for (nse = 0; nse < nsegs; nse++) {
    // gather the composing segments
    ppsegs[0] = &pgSegLl[start + nse];
    ppsegs[1] = &pgSegLu[start + nse];
    ppsegs[2] = &pgSegRl[start + nse];
    ppsegs[3] = &pgSegRu[start + nse];

    // define a bounding box enclosing all the sub-segments
    // (an own one is more convenient here than a G4BoundingBox object)
    xMin = xMax = ppsegs[0]->pPoly->GetPoints(0).x();
    yMin = yMax = ppsegs[0]->pPoly->GetPoints(0).y();
    zMin = zMax = ppsegs[0]->pPoly->GetPoints(0).z();
    for (ss = 0; ss < 4; ss++) {
      G4Point3D pts;
      npts = ppsegs[ss]->pPoly->GetnPoints();
      for (nn = 0; nn < npts; nn++) {
        pts = ppsegs[ss]->pPoly->GetPoints(nn);
        xMin = min(xMin, pts.x());
        xMax = max(xMax, pts.x());
        yMin = min(yMin, pts.y());
        yMax = max(yMax, pts.y());
        zMin = min(zMin, pts.z());
        zMax = max(zMax, pts.z());
      }
    }
    // some rounding of limits
    xMin = floor(xMin);
    yMin = floor(yMin);
    zMin = floor(zMin);
    xMax = ceil(xMax);
    yMax = ceil(yMax);
    zMax = ceil(zMax);

    G4double rr2, rc2 = 0., Rc2 = 0., zch = 0.;
    G4bool nowuseCylinder = this->useCylinder && pgons[iPg].cylinderMakesSense;
    if (nowuseCylinder) {
      rc2 = pow(pgons[iPg].tubr, 2.);
      Rc2 = pow(pgons[iPg].tubR, 2.);
      zch = pgons[iPg].zFace1 + pgons[iPg].thick;
    }

    G4Point3D theCenter(0., 0., 0.);
    G4Point3D myPoint;
    G4double theVolume;
    G4int nseen = 0;
    G4int ntot = 0;
    G4int ssold = -1;
    G4double xx, yy, zz;
    G4double dd1 = step;
    G4double dd2 = dd1 / 2.;
    EInside inside;
    G4ThreeVector pnt;

    for (zz = zMin + dd2; zz < zMax; zz += dd1) {
      for (yy = yMin + dd2; yy < yMax; yy += dd1) {
        for (xx = xMin + dd2; xx < xMax; xx += dd1) {
          myPoint.set(xx, yy, zz);
          ntot++;
          // Intersection with the cylinder: we check directly
          // (however, we neglect the passive areas)
          if (nowuseCylinder) {
            rr2 = myPoint.x() * myPoint.x() + myPoint.y() * myPoint.y();
            if (rr2 > Rc2)
              continue; // outside cyl
            if ((myPoint.z() > zch) && (rr2 < rc2))
              continue; // coax hole
          }

          // check if inside one of the sub-segments composing this segment
          pnt = G4ThreeVector(myPoint);
          if (ssold >= 0) { // start with the old one
            inside = ppsegs[ssold]->pPoly->Inside(pnt);
            if (inside == kInside) {
              nseen++;
              theCenter += myPoint;
              continue; // no need to check the others
            }
          }
          for (ss = 0; ss < 4; ss++) { // look into the (other) segments
            if (ss == ssold)
              continue;
            inside = ppsegs[ss]->pPoly->Inside(pnt);
            if (inside == kInside) {
              nseen++;
              theCenter += myPoint;
              ssold = ss;
              break; // no need to check the others
            }
          }
          ssold = -1; // outside all of them
        }
      }
    }

    theCenter = theCenter / nseen;
    theVolume = nseen * pow(dd1, 3.);
    printf(" Crystal# %2d  Segment# %3d : %6d shot %6d seen --> Volume = %6.3f cm3   Center = (%7.3f %7.3f %7.3f) cm\n",
           iPg, nse, ntot, nseen, theVolume / cm3, theCenter.x() / cm, theCenter.y() / cm, theCenter.z() / cm);
    segVolume[start + nse] = theVolume;
    segCenter[start + nse] = theCenter;
  }
  ppsegs.clear();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
//
void Hicari::InitializeScorers() {

  bool already_exist = false;
  m_HicariScorer = CheckScorer("HicariScorer", already_exist);

  if (already_exist)
    return;

  // Otherwise the scorer is initialised
  // vector<int> level({0, 1});
  int level = 1;
  /*
    m_HicariScorer->RegisterPrimitive(
            new CalorimeterScorers::PS_Calorimeter("Cristal",level, 0));
    G4SDManager::GetSDMpointer()->AddNewDetector(m_HicariScorer) ;
  */

  m_HicariScorer->RegisterPrimitive(new GeScorers::PS_GeDetector("Cristal", level, 1));
  G4SDManager::GetSDMpointer()->AddNewDetector(m_HicariScorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Hicari::Construct() { return (NPS::VDetector*)new Hicari(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_nps_Hicari {
 public:
  proxy_nps_Hicari() {
    NPS::DetectorFactory::getInstance()->AddToken("Hicari", "Hicari");
    NPS::DetectorFactory::getInstance()->AddDetector("Hicari", Hicari::Construct);
  }
};

proxy_nps_Hicari p_nps_Hicari;
}
