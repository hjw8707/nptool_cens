#ifndef TMUSETTSPECTRA_H
#define TMUSETTSPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: matta@lpccaen.in2p3.fr *
 *                                                                           *
 * Creation Date  : February 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for MUSETT                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TMUSETTData.h"
#include "TMUSETTPhysics.h"

// ForwardDeclaration
class TMUSETTPhysics ;

class TMUSETTSpectra:public VSpectra{
  public:
    // constructor and destructor
    TMUSETTSpectra();
    TMUSETTSpectra(std::map<int,int> DetectorIndex);
    ~TMUSETTSpectra();

  private:
    // Initialization methods
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  public:
    // Filling methods
    void FillRawSpectra(TMUSETTData*);
    void FillPreTreatedSpectra(TMUSETTData*);
    void FillPhysicsSpectra(TMUSETTPhysics*);

  private: // Information on MUSETT
    std::map<int,int> fDetectorToIndex;
    std::map<int,int> fIndexToDetector;
    unsigned int fNumberOfDetector;
    unsigned int fStripX;
    unsigned int fStripY;
    unsigned int fStripSecondLayer;
};

#endif
