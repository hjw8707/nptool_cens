/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elia Pilotto, Omar Nasr                                  *
 * contact address: pilottoelia@gmail.com, omar.nasr@etu.unicaen.fr          *
 *                                                                           *
 * Creation Date  : September 2021                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Use to kill the beam track and replace it with the reaction product       *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

//C++ libraries
#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>
#include <Randomize.hh>
#include <fstream>
//G4 libraries
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4EmCalculator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4IonTable.hh"
#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
//nptool libraries
#include "NPFunction.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "NPSFunction.hh"
//other
#include "GladFieldPropagation.hh"
#include "RootOutput.h"
#include "TLorentzVector.h"

////////////////////////////////////////////////////////////////////////////////
NPS::GladFieldPropagation::GladFieldPropagation(G4String modelName, G4Region* envelope)
  : G4VFastSimulationModel(modelName, envelope) {

  m_Map = NULL;
  m_Initialized = false;
}

////////////////////////////////////////////////////////////////////////////////
NPS::GladFieldPropagation::GladFieldPropagation(G4String modelName)
  : G4VFastSimulationModel(modelName) {}

////////////////////////////////////////////////////////////////////////////////
NPS::GladFieldPropagation::~GladFieldPropagation() {}

////////////////////////////////////////////////////////////////////////////////
void NPS::GladFieldPropagation::AttachReactionConditions() {
  // Reasssigned the branch address
  /*
  if (RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
    RootOutput::getInstance()->GetTree()->SetBranchAddress(
		  "ReactionConditions", &m_ReactionConditions);*/
}

////////////////////////////////////////////////////////////////////////////////
G4bool NPS::GladFieldPropagation::IsApplicable(const G4ParticleDefinition& particleType) {
  if (particleType.GetPDGCharge() == 0) return false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
G4bool NPS::GladFieldPropagation::ModelTrigger(const G4FastTrack& fastTrack) {
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void NPS::GladFieldPropagation::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) {
  if(!m_Initialized){
    m_Map = new GladFieldMap();

    //Needed for the RungeKutta method
    m_Map->SetCurrent(m_Current);
    m_Map->SetLimit( 10.*m / m_StepSize );
    m_Map->SetPropagationTimeInterval(0.8);
    m_Map->SetGladEntrance(m_GladEntrance.X(), m_GladEntrance.Y(), m_GladEntrance.Z());
    m_Map->SetGladTiltAngle(m_GladTiltAngle);
    m_Map->Set_MWPC3_Position(m_MW3_POS.X(), m_MW3_POS.Y(), m_MW3_POS.Z());
    m_Map->SetCentralTheta(m_MW3_CentralTheta);
    m_Map->LoadMap(m_FieldMap);
    m_Initialized = true;
  }
  
  RungeKuttaPropagation(fastTrack, fastStep);
  
  return;
}

/////////////////////////////////////////////////////////////////////////
void NPS::GladFieldPropagation::RungeKuttaPropagation (const G4FastTrack& fastTrack, G4FastStep& fastStep){

  static bool inside = false;//true if previous step is inside the volume
  static vector<TVector3> trajectory;
  static int count;//keeps track of the step reached inside trajectory

  // Get the track info
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();

  static G4VSolid* solid;//Stores the current solid 
  solid = PrimaryTrack->GetTouchable()->GetSolid();
  
  G4ThreeVector localDir = fastTrack.GetPrimaryTrackLocalDirection();
  G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();
  G4ThreeVector localMomentum = fastTrack.GetPrimaryTrackLocalMomentum();

  double speed = PrimaryTrack->GetVelocity();
  double charge = PrimaryTrack->GetParticleDefinition()->GetPDGCharge() / coulomb;
  double ConF_p = 1 / (joule * c_light / (m/s)) ; // MeV/c to kg*m/s (SI units)

  //Initially inside is false -> calculate trajectory
  if (!inside){
    count = 2;//skip first two positions as they are the same as the current position
    double Brho = localMomentum.mag() * ConF_p / charge;
    TVector3 pos (localPosition.x(), localPosition.y(), localPosition.z());
    TVector3 dir (localDir.x(), localDir.y(), localDir.z());

    trajectory.clear();
    trajectory = m_Map->Propagate(Brho, pos, dir,true);

    inside = true;
  }

  G4ThreeVector newPosition (trajectory[count].x(), trajectory[count].y(), trajectory[count].z());
  //benchmark
  //G4ThreeVector newDir = (newPosition - localPosition).unit();
  //G4ThreeVector newMomentum = newDir * localMomentum.mag(); 

  //Check if newPosition is not inside
  if (solid->Inside(newPosition) != kInside){
    inside = false;
    G4ThreeVector toOut = solid->DistanceToOut(localPosition, localDir) * localDir;
    newPosition = localPosition + toOut;
    //benchmark
    //newDir = (newPosition - localPosition).unit();
    //newMomentum = newDir * localMomentum.mag();

    //benchmark
    //if(single_particle) PrintData(m_StepSize, newPosition, newMomentum, outRKsp);
    //else PrintData(m_StepSize, newPosition, newMomentum, outRK);
    //counter++;
  }

  //benchmark
  G4ThreeVector newDir = (newPosition - localPosition).unit();
  G4ThreeVector newMomentum = newDir * localMomentum.mag(); 


  double time  = PrimaryTrack->GetGlobalTime()+(newPosition - localPosition).mag()/speed;

  fastStep.ProposePrimaryTrackFinalPosition( newPosition );
  fastStep.SetPrimaryTrackFinalMomentum ( newMomentum );
  fastStep.ProposePrimaryTrackFinalTime( time );

  count++;
  return;
}


