#ifndef STARK_h
#define STARK_h 1
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Serge Franchoo  contact address: franchoo@ipno.in2p3.fr  *
 *                                                                           *
 * Creation Date  : February 2018                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describes Tina simulation                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include <string>
#include <vector>

#include "G4AssemblyVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4ThreeVector.hh"
#include "NPInputParser.h"
#include "NPSVDetector.hh"
#include "TSTARKData.h"
#include "TSTARKRaw.h"

////////////////////////////////////////////////////////////////////////////////
// namespace for STARK
//  for some constants
////////////////////////////////////////////////////////////////////////////////
namespace STARKNS {
  extern G4double EnergyThreshold;
  ////////////////////////////////////////////////////////////
  // Resolution
  ////////////////////////////////////////////////////////////
  extern G4double X6_TRes;
  extern G4double X6_ERes;
  extern G4double BB10_TRes;
  extern G4double BB10_ERes;
  extern G4double QQQ5_TRes;
  extern G4double QQQ5_ERes;
  extern G4double CsI_TRes;
  extern G4double CsI_ERes;
  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // Geometry
  ////////////////////////////////////////////////////////////
  //
  // X6
  ////////////////////////////////////////////////////////////
  extern const G4double X6_PCBX;
  extern const G4double X6_PCBY;
  extern const G4double X6_PCBZ;
  extern const G4double X6_PCBSub1X;
  extern const G4double X6_PCBSub1Y;
  extern const G4double X6_PCBSub1Z;
  extern const G4double X6_PCBSub1XOffset;
  extern const G4double X6_PCBSub1YOffset;
  extern const G4double X6_PCBSub1ZOffset;
  extern const G4double X6_PCBSub2X;
  extern const G4double X6_PCBSub2Y;
  extern const G4double X6_PCBSub2Z;
  extern const G4double X6_PCBSub2XOffset;
  extern const G4double X6_PCBSub2YOffset;
  extern const G4double X6_PCBSub2ZOffset;

  extern const G4double X6_SiX;
  extern const G4double X6_SiY;
  extern const G4double X6_SiZ;
  extern const G4double X6_SiXOffset;
  extern const G4double X6_SiYOffset;
  extern const G4double X6_SiZOffset;
  extern const G4double X6_SiActiveX;
  extern const G4double X6_SiActiveY;
  extern const G4double X6_SiActiveZ;

  extern const G4int X6_NFrontStrips;
  extern const G4int X6_NBackStrips;
  ////////////////////////////////////////////////////////////
  //
  // BB10
  ////////////////////////////////////////////////////////////
  extern const G4double BB10_PCBX;
  extern const G4double BB10_PCBY;
  extern const G4double BB10_PCBZ;
  extern const G4double BB10_PCBSub1X;
  extern const G4double BB10_PCBSub1Y;
  extern const G4double BB10_PCBSub1Z;
  extern const G4double BB10_PCBSub1XOffset;
  extern const G4double BB10_PCBSub1YOffset;
  extern const G4double BB10_PCBSub1ZOffset;
  extern const G4double BB10_PCBSub2X;
  extern const G4double BB10_PCBSub2Y;
  extern const G4double BB10_PCBSub2Z;
  extern const G4double BB10_PCBSub2XOffset;
  extern const G4double BB10_PCBSub2YOffset;
  extern const G4double BB10_PCBSub2ZOffset;

  extern const G4double BB10_SiX;
  extern const G4double BB10_SiY;
  extern const G4double BB10_SiZ;
  extern const G4double BB10_SiXOffset;
  extern const G4double BB10_SiYOffset;
  extern const G4double BB10_SiZOffset;
  extern const G4double BB10_SiActiveX;
  extern const G4double BB10_SiActiveY;
  extern const G4double BB10_SiActiveZ;

  extern const G4int BB10_NFrontStrips;
  extern const G4int BB10_NBackStrips;
  ////////////////////////////////////////////////////////////
  //
  // QQQ5
  ////////////////////////////////////////////////////////////
  extern const G4double QQQ5_PCBOutR;
  extern const G4double QQQ5_PCBInR;
  extern const G4double QQQ5_PCBPhi0;
  extern const G4double QQQ5_PCBPhi1;
  extern const G4double QQQ5_PCBT;
  extern const G4double QQQ5_PCBCutX;
  extern const G4double QQQ5_PCBCutY;
  extern const G4double QQQ5_PCBCutZ;
  extern const G4double QQQ5_PCBCutXOffset;
  extern const G4double QQQ5_PCBCutYOffset;
  extern const G4double QQQ5_PCBCutZOffset;

  ////////////////////////////////////////////////////////////
  //
  // QQQ Wafer
  ////////////////////////////////////////////////////////////
  extern const G4double QQQ5_SiOutR;
  extern const G4double QQQ5_SiInR;
  extern const G4double QQQ5_SiT;
  extern const G4double QQQ5_SiPhi0;
  extern const G4double QQQ5_SiPhi1;
  extern const G4double QQQ5_SiActiveOutR;
  extern const G4double QQQ5_SiActiveInR;
  extern const G4double QQQ5_SiCut1X;
  extern const G4double QQQ5_SiCut1Y;
  extern const G4double QQQ5_SiCut1Z;
  extern const G4double QQQ5_SiCut1XOffset;
  extern const G4double QQQ5_SiCut1YOffset;
  extern const G4double QQQ5_SiCut1ZOffset;
  extern const G4double QQQ5_SiCut2X;
  extern const G4double QQQ5_SiCut2Y;
  extern const G4double QQQ5_SiCut2Z;
  extern const G4double QQQ5_SiCut2XOffset;
  extern const G4double QQQ5_SiCut2YOffset;
  extern const G4double QQQ5_SiCut2ZOffset;

  extern const G4int QQQ5_NRStrip;
  extern const G4int QQQ5_NAStrip;

  ////////////////////////////////////////////////////////////
  //
  // Connector
  ////////////////////////////////////////////////////////////
  extern const G4double Conn_X;
  extern const G4double Conn_Y;
  extern const G4double Conn_Z;

  ////////////////////////////////////////////////////////////
  //
  // CsI
  ////////////////////////////////////////////////////////////
  extern const G4double CsI_X;
  extern const G4double CsI_Y;
  extern const G4double CsI_Z;
  extern const G4double CsI_X6_XOffset;
  extern const G4double CsI_X6_YOffset1;
  extern const G4double CsI_X6_YOffset2;
  extern const G4double CsI_X6_ZOffset;
  extern const G4double CsI_BB10_XOffset;
  extern const G4double CsI_BB10_YOffset1;
  extern const G4double CsI_BB10_YOffset2;
  extern const G4double CsI_BB10_ZOffset;

  ////////////////////////////////////////////////////////////
  //
  // ANASEN CsI for QQQ3
  ////////////////////////////////////////////////////////////
  extern const G4double CsI_QQQ5_XOffset;
  extern const G4double CsI_QQQ5_YOffset;
  extern const G4double CsI_QQQ5_ZOffset;
  extern const G4double m_ANASENQQQ3CsIhypotenuse;
  extern const G4double m_ANASENQQQ3CsIThickness;
  extern const G4double m_ANASENQQQ3CsIWidthBot;
  extern const G4double m_ANASENQQQ3CsIWidthTop;
  extern const G4double m_ANASENQQQ3CsITotHyp;
  extern const G4double m_ANASENQQQ3CsITotHeight;
  extern const G4double m_ANASENQQQ3CsIHeight;
  ////////////////////////////////////////////////////////////

} // namespace STARKNS

////////////////////////////////////////////////////////////////////////////////

class STARK : public NPS::VDetector {
 public:
  STARK();
  virtual ~STARK();

  void AddDetector(string Type, G4ThreeVector POS, bool Flip, bool Rev, double Beta, int CsI, int Group, string mvName);
  void AddDetector(string Type, G4ThreeVector POS, G4ThreeVector Rot, int CsI, int Group, string mvName);
  G4AssemblyVolume* BuildX6Detector(int nCsI = 0);
  G4AssemblyVolume* BuildBB10Detector(int nCsI = 0);
  G4AssemblyVolume* BuildQQQ5Detector(int nCsI = 0);
  G4AssemblyVolume* BuildTarget();
  G4LogicalVolume* BuildSquareCsI();
  G4LogicalVolume* BuildANASENQQQ3CsI();

  // Reaction Region
  G4Region* m_ReactionRegion = nullptr;
  void SetReactionRegion(G4LogicalVolume* world);
  void SetMotherVolume(G4LogicalVolume* motherVolume);

  // Inherited from NPS::VDetector class /////////////
 public:
  // Read stream at Configfile to pick-up parameters of detector
  // called in DetectorConstruction::ReadDetectorConfiguration
  void ReadConfiguration(NPL::InputParser);

  // Construct detector and initialise sensitive part
  // called after DetectorConstruction::AddDetector
  void ConstructDetector(G4LogicalVolume* world);

  // Add detector branch to the EventTree
  // called after DetectorConstruction::AddDetector
  void InitializeRootOutput();

  // Read sensitive part and fill the Root tree
  // called in EventAction::EndOfEventAvtion
  void ReadSensitive(const G4Event* event);

  // Initialise all scorers used by the detector
  void InitializeScorers();
  G4MultiFunctionalDetector* m_X6Det;
  G4MultiFunctionalDetector* m_BB10Det;
  G4MultiFunctionalDetector* m_QQQ5Det;
  G4MultiFunctionalDetector* m_CsIDet;
  ////////////////////////////////////////////////////

 private:
  G4int HCID_X6;
  G4int HCID_BB10;
  G4int HCID_QQQ5;
  G4int HCID_CsI;

  G4AssemblyVolume* m_X6;
  G4AssemblyVolume* m_BB10;
  G4AssemblyVolume* m_QQQ5;
  G4AssemblyVolume* m_X6_wCsI[4];
  G4AssemblyVolume* m_BB10_wCsI[4];
  G4AssemblyVolume* m_QQQ5_wCsI[4];
  G4AssemblyVolume* m_Target;
  G4LogicalVolume* m_logicTarget;
  G4LogicalVolume* m_logicSquareCsI;
  G4LogicalVolume* m_logicANASENQQQ3CsI;

  // Event class to store data
  TSTARKData* m_Event;
  TSTARKRaw* m_Raw;

  // Type & Geometry
  vector<string> m_Type;
  vector<G4ThreeVector> m_Pos;
  vector<G4ThreeVector> m_Rot; // Rotation of the detector
  vector<bool>
      m_AutoRotateFacingBeamAxis; // If true, the detector will be rotated to face the beam axis. (do not use m_Rot)
  vector<bool> m_Flip;            // Which surface facing to the beam line.
  vector<bool> m_Rev;             // Reverse along the beam axis. (connector direction)
  vector<G4double> m_Beta;        // Rotation along the beam axis. (only for QQQ5)
  vector<G4int> m_Group;          // Detector group for dE-E analysis
  vector<int> m_CsI;              // Number of CsI layers
  vector<string> m_MVName;        // Name of the mother volume

  // Target
  bool m_useTarget;
  string m_TargetMaterial;
  double m_Pressure;
  double m_Temperature;
  double m_Radius;
  double m_Z;

  // Visualisation
  G4VisAttributes *m_VisX6, *m_VisX6PCB;
  G4VisAttributes *m_VisBB10, *m_VisBB10PCB;
  G4VisAttributes *m_VisQQQ5, *m_VisQQQ5PCB;
  G4VisAttributes *m_VisConn, *m_VisTarget;
  G4VisAttributes* m_VisCsI;

 public:
  // Dynamic loading of the library
  static NPS::VDetector* Construct();
};
#endif
