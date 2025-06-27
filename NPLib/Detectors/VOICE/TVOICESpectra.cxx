/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: shbae  contact address: shbae2703@gmail.com                        *
 *                                                                           *
 * Creation Date  : March 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold VOICE Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// class header
#include "TVOICESpectra.h"

// STL
#include <iostream>
#include <string>
using namespace std;

// NPTool header
#include "NPOptionManager.h"

////////////////////////////////////////////////////////////////////////////////
TVOICESpectra::TVOICESpectra() : fNumberOfDetectors(0) { SetName("VOICE"); }

////////////////////////////////////////////////////////////////////////////////
TVOICESpectra::TVOICESpectra(unsigned int NumberOfDetectors) {
    if (NPOptionManager::getInstance()->GetVerboseLevel() > 0)
        cout << "************************************************" << endl
             << "TVOICESpectra : Initalizing control spectra for " << NumberOfDetectors << " Detectors" << endl
             << "************************************************" << endl;
    SetName("VOICE");
    fNumberOfDetectors = NumberOfDetectors;

    InitRawSpectra();
    InitPreTreatedSpectra();
    InitPhysicsSpectra();
}

////////////////////////////////////////////////////////////////////////////////
TVOICESpectra::~TVOICESpectra() {}

////////////////////////////////////////////////////////////////////////////////
void TVOICESpectra::InitRawSpectra() {
    static string name;
    for (unsigned int i = 0; i < fNumberOfDetectors; i++) {  // loop on number of detectors
        // Energy
        name = "VOICE" + NPL::itoa(i + 1) + "_ENERGY_RAW";
        AddHisto1D(name, name, 4096, 0, 16384, "VOICE/RAW");
        // Time
        name = "VOICE" + NPL::itoa(i + 1) + "_TIME_RAW";
        AddHisto1D(name, name, 4096, 0, 16384, "VOICE/RAW");
    }  // end loop on number of detectors
}

////////////////////////////////////////////////////////////////////////////////
void TVOICESpectra::InitPreTreatedSpectra() {
    static string name;
    for (unsigned int i = 0; i < fNumberOfDetectors; i++) {  // loop on number of detectors
        // Energy
        name = "VOICE" + NPL::itoa(i + 1) + "_ENERGY_CAL";
        AddHisto1D(name, name, 500, 0, 25, "VOICE/CAL");
        // Time
        name = "VOICE" + NPL::itoa(i + 1) + "_TIME_CAL";
        AddHisto1D(name, name, 500, 0, 25, "VOICE/CAL");

    }  // end loop on number of detectors
}

////////////////////////////////////////////////////////////////////////////////
void TVOICESpectra::InitPhysicsSpectra() {
    static string name;
    // Kinematic Plot
    name = "VOICE_ENERGY_TIME";
    AddHisto2D(name, name, 500, 0, 500, 500, 0, 50, "VOICE/PHY");
}

////////////////////////////////////////////////////////////////////////////////
void TVOICESpectra::FillRawSpectra(TVOICEData* RawData) {
    static string name;
    static string family;

    // Energy
    unsigned int sizeE = RawData->GetMult();
    for (unsigned int i = 0; i < sizeE; i++) {
        name = "VOICE" + NPL::itoa(i) + "_ENERGY_RAW";
        family = "VOICE/RAW";

        FillSpectra(family, name, RawData->GetEnergy(i));
    }

    // Time
    unsigned int sizeT = RawData->GetMult();
    for (unsigned int i = 0; i < sizeT; i++) {
        name = "VOICE" + NPL::itoa(i) + "_TIME_RAW";
        family = "VOICE/RAW";

        FillSpectra(family, name, RawData->GetTime(i));
    }
}

////////////////////////////////////////////////////////////////////////////////
void TVOICESpectra::FillPreTreatedSpectra(TVOICEData* PreTreatedData) {
    static string name;
    static string family;

    // Energy
    unsigned int sizeE = PreTreatedData->GetMult();
    for (unsigned int i = 0; i < sizeE; i++) {
        name = "VOICE" + NPL::itoa(i) + "_ENERGY_CAL";
        family = "VOICE/CAL";

        FillSpectra(family, name, PreTreatedData->GetEnergy(i));
    }

    // Time
    unsigned int sizeT = PreTreatedData->GetMult();
    for (unsigned int i = 0; i < sizeT; i++) {
        name = "VOICE" + NPL::itoa(i) + "_TIME_CAL";
        family = "VOICE/CAL";

        FillSpectra(family, name, PreTreatedData->GetTime(i));
    }
}

////////////////////////////////////////////////////////////////////////////////
void TVOICESpectra::FillPhysicsSpectra(TVOICEPhysics* Physics) {
    static string name;
    static string family;
    family = "VOICE/PHY";

    // Energy vs time
    // unsigned int sizeE = Physics->Energy.size();
    // for (unsigned int i = 0; i < sizeE; i++) {
    //     name = "VOICE_ENERGY_TIME";
    //     FillSpectra(family, name, Physics->Energy[i], Physics->Time[i]);
    // }
}
