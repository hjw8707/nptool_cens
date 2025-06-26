/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : november 2009                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold CsITl  Physics                                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

//   NPL
#include "TCsITlPhysics.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "NPDetectorFactory.h"
#include "NPCalibrationManager.h"
#include "NPOptionManager.h"

//   STL
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <stdlib.h>
#include <cmath>
using namespace std;

//   ROOT
#include "TChain.h"
#include "TVector3.h"

ClassImp(TCsITlPhysics)
///////////////////////////////////////////////////////////////////////////
TCsITlPhysics::TCsITlPhysics()
{
    EventData = new TCsITlData ;
    EventPhysics = this ;
}

///////////////////////////////////////////////////////////////////////////
TCsITlPhysics::~TCsITlPhysics()
{}

///////////////////////////////////////////////////////////////////////////
void TCsITlPhysics::Clear(){
  nhit = 0; }


///////////////////////////////////////////////////////////////////////////
void TCsITlPhysics::ReadConfiguration(NPL::InputParser parser){
}


///////////////////////////////////////////////////////////////////////////
void TCsITlPhysics::InitializeRootInputRaw()
{
    TChain* inputChain = RootInput::getInstance()->GetChain()     ;
    inputChain->SetBranchStatus ( "CsITl"       , true )        ;
    inputChain->SetBranchAddress( "CsITl"       , &EventData )  ;
}
///////////////////////////////////////////////////////////////////////////
void TCsITlPhysics::InitializeRootInputPhysics()
{
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus ( "CsITl", true );
    inputChain->SetBranchAddress( "CsITl", &EventPhysics );
}
///////////////////////////////////////////////////////////////////////////
void TCsITlPhysics::InitializeRootOutput()
{
    TTree* outputTree = RootOutput::getInstance()->GetTree()            ;
    outputTree->Branch( "CsITl" , "TCsITlPhysics" , &EventPhysics ) ;
    outputTree->Branch("nhit" ,&nhit,"nhit/I");
    outputTree->Branch("detN" ,detN, "detN[nhit]/I");
    outputTree->Branch("E",E,"E[nhit]/D");
    outputTree->Branch("T",T,"T[nhit]/D");
}

///////////////////////////////////////////////////////////////////////////
void TCsITlPhysics::BuildPhysicalEvent()
{
    BuildSimplePhysicalEvent()   ;
}

///////////////////////////////////////////////////////////////////////////
void TCsITlPhysics::BuildSimplePhysicalEvent()
{
  Clear();
  
  ////////////////////////////////////////////////////////////
  // Loop for All
  nhit = EventData->GetMult();
  for (Int_t i = 0 ; i < EventData->GetMult() ; i++) {
    detN[i] = EventData->GetDetN(i);
    E   [i] = EventData->GetE(i);
    T   [i] = EventData->GetT(i);
  }    
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TCsITlPhysics::Construct(){
    return (NPL::VDetector*) new TCsITlPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
    class proxy_csi{
    public:
        proxy_csi(){
	  NPL::DetectorFactory::getInstance()->AddToken("CsITl","CsITl");
            NPL::DetectorFactory::getInstance()->AddDetector("CsITl",TCsITlPhysics::Construct);
        }
    };
  
    proxy_csi p;
}

