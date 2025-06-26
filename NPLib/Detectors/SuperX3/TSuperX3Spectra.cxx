/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 ****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville         address: deserevi@ipno.in2p3.fr  *
 *                                                                           *
 * Creation Date  : November 2015                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for SuperX3                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// class header
#include "TSuperX3Spectra.h"

// C++ headers
#include <iostream>
#include <string>
using namespace std;

// NPTool headers
#include "NPOptionManager.h"

////////////////////////////////////////////////////////////////////////////////
TSuperX3Spectra::TSuperX3Spectra()
    : fNumberOfDetectors(0), fNumberOfStripsFront(16), fNumberOfStripsBack(16), fNumberOfCounters(10) {
  SetName("SuperX3");
}

////////////////////////////////////////////////////////////////////////////////
TSuperX3Spectra::TSuperX3Spectra(unsigned int NumberOfDetectors) {
  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}

////////////////////////////////////////////////////////////////////////////////
TSuperX3Spectra::~TSuperX3Spectra() {}

////////////////////////////////////////////////////////////////////////////////
void TSuperX3Spectra::InitRawSpectra() {}

////////////////////////////////////////////////////////////////////////////////
void TSuperX3Spectra::InitPreTreatedSpectra() {}

////////////////////////////////////////////////////////////////////////////////
void TSuperX3Spectra::InitPhysicsSpectra() {}

////////////////////////////////////////////////////////////////////////////////
void TSuperX3Spectra::FillRawSpectra(TSuperX3Data* RawData) {}

////////////////////////////////////////////////////////////////////////////////
void TSuperX3Spectra::FillPreTreatedSpectra(TSuperX3Data* PreTreatedData) {}

////////////////////////////////////////////////////////////////////////////////
void TSuperX3Spectra::FillPhysicsSpectra(TSuperX3Physics* Physics) {}

