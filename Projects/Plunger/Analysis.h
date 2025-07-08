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
 *  This class describe  Plunger analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include "NPVAnalysis.h"
#include "TASGARDData.h"
#include "TASGARDPhysics.h"
#include "TInitialConditions.h"
#include "TPlungerData.h"
#include "TPlungerPhysics.h"
#include "TReactionConditions.h"
class Analysis : public NPL::VAnalysis {
   public:
    Analysis();
    ~Analysis();

   public:
    void Init();
    void TreatEvent();
    void End();

    void InitializeRootInput();
    void InitializeRootOutput();

    static NPL::VAnalysis* Construct();

   private:
    TInitialConditions* InitialConditions;
    TReactionConditions* ReactionConditions;
    TPlungerData* PlungerData;
    TASGARDData* ASGARDData;

   private:
    double flagVelocity;
    double flagKineticEnergy;

    TPlungerPhysics* PlungerPhysics;
    TASGARDPhysics* ASGARDPhysics;
};
#endif
