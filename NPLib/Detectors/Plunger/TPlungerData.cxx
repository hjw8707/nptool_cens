/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: hjw8707@gmail.com                        *
 *                                                                           *
 * Creation Date  : 7월 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Plunger Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "TPlungerData.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

ClassImp(TPlungerData)

    //////////////////////////////////////////////////////////////////////
    TPlungerData::TPlungerData() {}

//////////////////////////////////////////////////////////////////////
TPlungerData::~TPlungerData() {}

//////////////////////////////////////////////////////////////////////
void TPlungerData::Clear() {
    // Particle flying to the stopper
    fPlunger_ParticleName.clear();
    fPlunger_Velocity.clear();
    fPlunger_KineticEnergy.clear();
    fPlunger_Position.clear();
}

//////////////////////////////////////////////////////////////////////
void TPlungerData::Dump() const {
    // This method is very useful for debuging and worth the dev.
    cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TPlungerData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

    // Particle flying to the stopper
    size_t mysize = fPlunger_ParticleName.size();
    cout << "Plunger_ParticleName_Mult: " << mysize << endl;

    for (size_t i = 0; i < mysize; i++) {
        cout << "ParticleName: " << fPlunger_ParticleName[i] << " Velocity: " << fPlunger_Velocity[i]
             << " KineticEnergy: " << fPlunger_KineticEnergy[i] << " Position: " << fPlunger_Position[i].X() << " "
             << fPlunger_Position[i].Y() << " " << fPlunger_Position[i].Z() << endl;
    }
}
