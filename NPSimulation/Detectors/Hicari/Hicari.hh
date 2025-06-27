#ifndef Hicari_h
#define Hicari_h 1
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Goigoux  contact address: thomas.goigoux@cea.fr                        *
 *                                                                           *
 * Creation Date  : july 2019                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Hicari simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ header
#include <string>
#include <vector>
using namespace std;

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "Hicari_Helper.hh"
#include "CConvexPolyhedron.hh"
#include "THicariData.h"
#include "NPInputParser.h"

#define MAX_SEGS 8              /* max. number of segments to take in events */
#define MAX_INTPTS (2*MAX_SEGS) /* max. number of interaction points */
#define FIRSTMB 800
#define MBSEGS 6
#define MBCLUST 6
#define CLOVERS 4
#define CLSEGS 4

class Hicari : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Hicari() ;
    virtual ~Hicari() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Spherical
    void AddDetector(double X,double Y, double Z, double ThetaX, double ThetaY, double ThetaZ);  

    G4int InitializeMaterials();
    void BuildClover(int i_clo, G4LogicalVolume* world);
    void BuildSideCatcher();
    void BuildBackCatcher();
    void BuildSideShield();
    void BuildCollimator();
  
    ////////////////////////////////////////////////////
    //////  Inherite from NPS::VDetector class /////////
    ////////////////////////////////////////////////////
  public:
    // Read stream at Configfile to pick-up parameters of detector (Position,...)
    // Called in DetecorConstruction::ReadDetextorConfiguration Method
    void ReadConfiguration(NPL::InputParser) ;


    // Construct detector and inialise sensitive part.
    // Called After DetecorConstruction::AddDetector Method
    void ConstructDetector(G4LogicalVolume* world) ;

    // Add Detector branch to the EventTree.
    // Called After DetecorConstruction::AddDetector Method
    void InitializeRootOutput() ;

    // Read sensitive part and fill the Root tree.
    // Called at in the EventAction::EndOfEventAvtion
    void ReadSensitive(const G4Event* event) ;

  public:   // Scorer
    //   Initialize all Scorer used by the MUST2Array
    void InitializeScorers() ;

    //   Associated Scorer
    G4MultiFunctionalDetector* m_HicariScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    THicariData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
 /*
  private: // Geometry
    // Detector Coordinate (of the front of the Germanium)
    vector<double>  m_X; 
    vector<double>  m_Y; 
    vector<double>  m_Z; 
    // Detector orientation
    vector<double>  m_ThetaX; //rotation angles to X, Y, Z axis
    vector<double>  m_ThetaY;
    vector<double>  m_ThetaZ;
*/

  /////////////////////////////////////////////////
  /// Files from which actual geometry is read 
  ///////////////////////////////////////////////// 
  private:
    G4String                    iniPath;     //> directory where the files are located
    G4String                    eulerFile;   //> angles and positions to place the clusters into space
    G4String                    solidFile;   //> shape of the crystals
    G4String                    sliceFile;   //> segmentation
    G4String                    wallsFile;   //> cryostats
    G4String                    clustFile;   //> arrangement of the crystals and cryostats within a cluster
    int MBsegments[6] = {2, 3, 5, 4, 1, 0};
    int CLsegments[4] = {2, 1, 0, 3};
    
  private:
    G4String                    directoryName;  //> for the command line  
    
  /////////////////////////////////////
  /// materials (pointers and names)
  ///////////////////////////////////// 
  private:
    G4Material                  *matCryst;       //> crystals
    G4Material                  *matWalls;       //> encapsulation
    G4Material                  *matBackWalls;   //> behind the crystals
    G4Material                  *matHole;        //> vacuum within the cryostat
    G4Material                  *matCryo;        //> cryostats

  private:
    G4String                    matCrystName;     //> crystals
    G4String                    matWallsName;     //> cryostats and encapsulation
    G4String                    matBackWallsName; //> behind the crystals
    G4String                    matHoleName;      //> vacuum within the cryostat
    G4String                    matCryoName;      //> cryostats

  ///////////////////////////////////////////////////////////////////////
  /// structures needed to store geometry data during the construction  
  ///////////////////////////////////////////////////////////////////////
  private:
    std::vector<CeulerAngles>   euler;        //> angles and positions to place the clusters into space
    std::vector<CpolyhPoints>   pgons;        //> shape of the crystals
    std::vector<CpolyhPoints>   walls;        //> cryostats
    std::vector<CclusterAngles> clust;        //> arrangement of the crystals and cryostats within a cluster
    std::vector<CpolyhPoints>   capsO;        //> encapsulation (outer size)
    std::vector<CpolyhPoints>   capsI;        //> encapsulation (inner size)

  private:
    G4int                       nEuler;       //> number of clusters composing the array
    G4int                       nPgons;       //> number of different crystal shapes within the array
    G4int                       nClAng;       //> number of crystals composing a cluster
    G4int                       nWalls;       //> number of cryostat parts within a cluster
    G4int                       maxPgons;     //> maximum index of crystal shapes
    G4int                       nDets;
    G4int                       nClus;
    G4int                       iCMin;
    G4int                       iCMax;
    G4int                       iGMin;
    G4int                       iGMax;
    G4int                       maxSec;
    G4int                       maxSli;

  private:
    std::vector<G4int>          crystType;    //> lookup table detector number --> crystal shape
    std::vector<G4int>          planarLUT;    //> lookup table detector number --> planar or not

  private:
    G4int                       nWlTot;       //> total number of cryostat parts within the array
    G4int                       maxSolids;    //> maximum number of solids within a cluster

  /////////////////////////////////////////////////
  /// structures needed to build the segmentation
  ////////////////////////////////////////////////
  private:
    std::vector<CpolyhPoints>   pgSegLl;      //> segments on lower Left  side of edges
    std::vector<CpolyhPoints>   pgSegLu;      //> segments on upper Left  side of edges
    std::vector<CpolyhPoints>   pgSegRl;      //> segments on lower Right side of edges
    std::vector<CpolyhPoints>   pgSegRu;      //> segments on upper Right side of edges

  private:
    std::vector<G4int>          nSegments;
    std::vector<G4int>          tSegments;
    G4int                       totSegments;
    G4int                       nSeg;
    
  private:
    std::vector<G4double>       segVolume;    //> the volume of the (composite) segment
    std::vector<G4Point3D>      segCenter;    //> the center of mass of the (composite) segment

  private:
    G4int                       stepFactor;    //> integration step for the calculation of segment volume
    G4bool                      stepHasChanged;//> true: integration step was changed and segment volumes
                                               //>       should be recomputed
   
  ////////////////////////////////////////////
  /// size of the equivalent germanium shell
  ////////////////////////////////////////////
  private:
    G4double                    arrayRmin;     //> inner radius
    G4double                    arrayRmax;     //> outer radius
    
  ///////////////////////////////////////////
  /// rotation applied to the whole array
  //////////////////////////////////////////
  private:
    G4double                    thetaShift;   //> theta
    G4double                    phiShift;     //> phi
    G4double                    thetaPrisma;  //> thetaPrisma

  ///////////////////////////////////////////
  /// traslation applied to the whole array
  //////////////////////////////////////////
  private:
    G4ThreeVector               posShift;   

  /////////////////
  /// some flags 
  //////////////// 
  private:
    G4bool                      usePassive;     //> true: passive areas of the crystals will be generated
    G4bool                      drawReadOut;    //> true: segments will be visualized
    G4bool                      useAncillary;   //> true: ancillary detectors will be constructed
    G4bool                      makeCapsule;    //> true: encapsulation will be generated
    G4bool                      useCylinder;    //> true: the intersection with the cylinder will be considered
    G4bool                      cryostatStatus; //> true: include cryostats behind clusters //LR
    G4bool                      readOut;        //> true: a segmentation have been defined  
    G4bool                      printVolumes;   //> true: report crystal, segment, and passive volumes

  ///////////////////////////////////////////                               //LR
  /// Cryostats             
  //////////////////////////////////////////
  public:
    void SetCryostats(){cryostatStatus = true;};                            //LR

  private:
    G4ThreeVector               cryostatPos0;                               //LR
    G4ThreeVector               cryostatPos;                                //LR
    G4RotationMatrix            cryostatRot;                                //LR
  //G4double                    cryostatRadius;                             //LR
  //G4double                    cryostatLength;                             //LR
    G4double                    cryostatZplanes[7];                         //LR
    G4double                    cryostatRinner[7];                          //LR
    G4double                    cryostatRouter[7];                          //LR

  //////////////////////////////
  ///////// Methods  ///////////
  //////////////////////////////
  private:
    void     InitData();

  //////////////////////////
  /// read the input files
  /////////////////////////
  private:
    void      ReadEulerFile();
    void      ReadSolidFile();
    void      ReadSliceFile();
    void      ReadWallsFile();
    void      ReadClustFile();
  
  //////////////////////////////////////////////////////////////
  /// look for the materials starting from the material names
  /////////////////////////////////////////////////////////////
  private:
    G4int     FindMaterials();

  /////////////////////////////////////////////////////////
  /// Construct the various elements composing the array
  /////////////////////////////////////////////////////////  
  private:
    void      ConstructGeCrystals  ();    
    void      ConstructTheCapsules ();    
    void      ConstructTheClusters ();    
    void      ConstructTheWalls    ();
    
  ////////////////////////////////
  /// placement of the elements  
  ////////////////////////////////
  private:
    void      PlaceTheClusters     (G4LogicalVolume*);
    void      ShowStatus ();

  //////////////////////////////////
  /// Construction of the segments
  //////////////////////////////////
  
  private:
    //////////////////////////////
    /// Construct the segments
    //////////////////////////////
    void      ConstructSegments        ();
    /////////////////////////////////////////////////////////////////////////////////
    /// Calculate the vertexes of the segments starting from the original polyhedra
    /////////////////////////////////////////////////////////////////////////////////
    G4int     CalculateSegments        ( G4int );
    ///////////////////////////////////
    /// Checks for possible overlaps
    //////////////////////////////////
    G4int     CheckOverlap             ( G4int, G4int, G4int );
    /////////////////////////////////////////////////////////////////////////////////////////
    /// Calculates volume and center of the segments (each of them composed of more parts!)
    //////////////////////////////////////////////////////////////////////////////////////////
    void      CalculateVolumeAndCenter ( G4int, G4int, G4int, G4double );
    ////////////////////////////////////////////////////////////////////////////////////////
    /// Calculates the intersection between a plane and a line passing through two points
    ///////////////////////////////////////////////////////////////////////////////////////
    G4Point3D XPlaneLine               ( const G4Plane3D &vv, const G4Point3D &pA,  const G4Point3D &pB );
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Calculates segment number (corresponding to a given detector and position relative to the crystal)
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    G4int     GetCoaxSegmentNumber         ( G4int, G4ThreeVector );      

  //////////////////////////////////
  /// writes out the information
  /////////////////////////////////
  private:
    void WritePositions               ( std::ofstream &outFileLMD, G4double=1.*mm );
    //void WriteSegmentPositions        ( std::ofstream &outFileLMD, G4double=1.*mm );
    void WriteSegmentPositions        ( G4String, G4double=1*mm );
    //void WriteCrystalPositions        ( std::ofstream &outFileLMD, G4double=1.*mm );
    void WriteCrystalPositions        ( G4String, G4double=1*mm );
    void WriteCrystalTransformations  ( std::ofstream &outFileLMD, G4double=1.*mm );
    void WriteHeader                  ( std::ofstream &outFileLMD, G4double=1.*mm );
    void WriteSegmentAngles           ( G4String, G4int = 0 );
    void WriteCrystalAngles           ( G4String );

  //////////////////////////////////////////////
  /// public interface to the private method!
  ////////////////////////////////////////////////
  public:
    G4int     GetSegmentNumber         ( G4int, G4ThreeVector );    

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
