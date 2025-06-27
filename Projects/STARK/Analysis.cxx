/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: XAUTHORX  contact address: XMAILX                        *
 *                                                                           *
 * Creation Date  : XMONTHX XYEARX                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  TiNA analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
#include <numeric>

using namespace std;
#include "Analysis.h"
#include "NPFunction.h"
#include "NPAnalysisFactory.h"
#include "NPDetectorManager.h"
#include "NPOptionManager.h"

#include "TString.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis():
  BeamReacE(0),
  nearHitIdx(-1),
  dE(0), E(0), Ex(0), QValue(0), ThetaLab(0), ThetaCM(0),
  X(0), Y(0), Z(0),
  myInit(NULL), myReac(NULL)
{}

////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  // initialize input and output branches
  cout << "!!!!!!!!!!!!!!! Initializing analysis !!!!!!!!!!!!!!!!!!" << endl;
  InitOutputBranch();
  InitInputBranch();
    
  // get STARK objects
  std::vector<std::string> detList = m_DetectorManager->GetDetectorList();
  for (auto it = detList.begin() ; it != detList.end() ; it++) {
    if      ((*it) == "STARK") 
      STARK  = static_cast<TSTARKPhysics*> (m_DetectorManager->GetDetector("STARK"));}
  
  // get reaction information
  myReaction.ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  
  // initialize various parameters
  Rand = TRandom3();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  // Reinitiate calculated variable
  ReInitValue();

  if (myReac) {
    ReacVertex = myReac->GetVertexPosition();
    BeamReacE = myInit->GetIncidentInitialKineticEnergy(); 
    myReaction.SetBeamEnergy(BeamReacE);}

  if (!STARK) return;
  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////// LOOP for STARK ////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  // dE-E analysis
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  // Simple criterion
  // 1. E  = the total E of all hit
  // 2. dE = the E of the hit at BB10
  ////////////////////////////////////////////////////////////////////////////

  int nhit = STARK->nhit;
  if (nhit == 0) return; // no hits

  ////////////////////////////////////////////////////////////////////////////
  for (int i = 0 ; i < nhit ; i++) {
    if (STARK->type[i] == 1) { // hit at BB10
      nearHitIdx = i; break; }
  }
  if (nearHitIdx < 0) return;
  
  TVector3 HitPosition = STARK->hPosArr[nearHitIdx];
  TVector3 HitDirection = HitPosition; // assumption = ReacVertex == (0,0,0)
  X = HitPosition.X();
  Y = HitPosition.Y();
  Z = HitPosition.Z();
  
  //==================================================//
  // Part 1 : Impact Angle
  ThetaLab = HitDirection.Theta(); // assumption = BeamDirection == (0,0,1)
        
  //==================================================//
  // Part 2 : Impact Energy
  E = 0;
  for (int i = 0 ; i < nhit ; i++) { // total energy of all hits
    E += STARK->sumE[i]; }
                
  dE = STARK->sumE[nearHitIdx];
  //==================================================//
  
  //==================================================//
  // Part 3 : Excitation Energy Calculation
  Ex = myReaction.ReconstructRelativistic( E, ThetaLab );
  //==================================================//
        
  //==================================================//
  // Part 4 : Theta CM Calculation
  ThetaCM  = myReaction.EnergyLabToThetaCM( E, ThetaLab)/deg;
  ThetaLab = ThetaLab/deg;
  //==================================================//
        
  //==================================================//
  // Part 5 : QValue Calculation
  QValue = myReaction.GetQValue();
  QValue = QValue - Ex;
  //==================================================//
  return;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("ReacVertex","TVector3",&ReacVertex);
  RootOutput::getInstance()->GetTree()->Branch("BeamReacE",&BeamReacE,"BeamReacE/D");

  RootOutput::getInstance()->GetTree()->Branch("nearHitIdx",&nearHitIdx,"nearHitIdx/I");
    
  RootOutput::getInstance()->GetTree()->Branch("dE",&dE,"dE/D");
  RootOutput::getInstance()->GetTree()->Branch("E",&E,"E/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab,"ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");

  RootOutput::getInstance()->GetTree()->Branch("X",&X,"X/D");
  RootOutput::getInstance()->GetTree()->Branch("Y",&Y,"Y/D");
  RootOutput::getInstance()->GetTree()->Branch("Z",&Z,"Z/D");
  RootOutput::getInstance()->GetTree()->Branch("QValue",&QValue,"QValue/D");
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
    RootInput::getInstance()->GetChain()->SetBranchStatus("InitialConditions",true);
    RootInput::getInstance()->GetChain()->SetBranchAddress("InitialConditions",&myInit);
    myInit = new TInitialConditions();
    if (RootInput::getInstance()->GetChain()->GetBranch("ReactionConditions")) {
      RootInput::getInstance()->GetChain()->SetBranchStatus("ReactionConditions",true);
      RootInput::getInstance()->GetChain()->SetBranchAddress("ReactionConditions",&myReac);
      myReac = new TReactionConditions();}
    else
      myReac = NULL;}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  ReacVertex.SetXYZ(0,0,0);
  BeamReacE = 0;

  nearHitIdx = -1;
  dE =       TMath::QuietNaN();
  E  =       TMath::QuietNaN();
  Ex =       TMath::QuietNaN();
  QValue =   TMath::QuietNaN();
  ThetaLab = TMath::QuietNaN();
  ThetaCM =  TMath::QuietNaN();
  X =        TMath::QuietNaN();
  Y =        TMath::QuietNaN();
  Z =        TMath::QuietNaN();}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
    return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_analysis{
public:
    proxy_analysis(){
        NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
    }
};

proxy_analysis p_analysis;
}

