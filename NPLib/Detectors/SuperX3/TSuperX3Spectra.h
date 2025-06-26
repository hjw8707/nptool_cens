#ifndef TSuperX3SPECTRA_H
#define TSuperX3SPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 ****************************************************************************/

/*****************************************************************************
 * Original Author: n. de Sereville        address: deserevi@ipno.in2p3.fr   *
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

// NPLib headers
#include "NPVSpectra.h"
#include "TSuperX3Data.h"
#include "TSuperX3Physics.h"

// ForwardDeclaration
class TSuperX3Physics;

class TSuperX3Spectra : public VSpectra {
 public:
  // constructor and destructor
  TSuperX3Spectra();
  TSuperX3Spectra(unsigned int NumberOfDetectors);
  ~TSuperX3Spectra();

 private:
  // Initialization methods
  void InitRawSpectra();
  void InitPreTreatedSpectra();
  void InitPhysicsSpectra();

 public:
  // Filling methods
  void FillRawSpectra(TSuperX3Data*);
  void FillPreTreatedSpectra(TSuperX3Data*);
  void FillPhysicsSpectra(TSuperX3Physics*);

 private: // Information on SuperX3
  unsigned int fNumberOfDetectors;
  unsigned int fNumberOfStripsFront;
  unsigned int fNumberOfStripsBack;
  Int_t fNumberOfCounters;
};

#endif
