#ifndef Analysis_h 
#define Analysis_h
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

#include "NPVAnalysis.h"
#include "NPEnergyLoss.h"
#include "NPReaction.h"
#include "NPBeam.h"
#include "RootOutput.h"
#include "RootInput.h"
#include "TSTARKPhysics.h"
#include "TInitialConditions.h"
#include "TReactionConditions.h"
#include <TRandom3.h>
#include <TVector3.h>
#include <TMath.h>
#include <vector>

class Analysis: public NPL::VAnalysis{
  public:
    Analysis();
    ~Analysis();

  public: 
    void Init();
    void TreatEvent();
    void End();

    void InitOutputBranch();
    void InitInputBranch();
    void ReInitValue();
    static NPL::VAnalysis* Construct();

  private:
    ////////////////////////////////////////////////////////////
    // Vertex
    TVector3 ReacVertex;
    // Beam energy at Reaction
    double BeamReacE;
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // for STARK
    int nearHitIdx;
    double dE, E, Ex, QValue, ThetaLab, ThetaCM;
    
    NPL::Reaction myReaction;

    TRandom3 Rand ;
    double X, Y, Z ;
    // Branches and detectors
    TSTARKPhysics* STARK;

    ////////////////////////////////////////////////////////////
    
    TInitialConditions* myInit ;
    TReactionConditions* myReac ;

};
#endif
