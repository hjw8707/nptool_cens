#ifndef TISSSPECTRA_H
#define TISSSPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 ****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 * Author: M. Labiche                     address: marc.labiche@stfc.ac.uk   *
 *                                                                           *
 * Creation Date  : July 2019                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for Iss                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TIssData.h"
#include "TIssPhysics.h"

// ForwardDeclaration
class TIssPhysics ;

class TIssSpectra:public VSpectra {
  public:
    // constructor and destructor
    TIssSpectra();
    TIssSpectra(unsigned int NumberOfDetector);
    ~TIssSpectra();

  private:
    // Initialization methods
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  public:
    // Filling methods
    void FillRawSpectra(TIssData*);
    void FillPreTreatedSpectra(TIssData*);
    void FillPhysicsSpectra(TIssPhysics*);

  private: // Information on Iss
    unsigned int fNumberOfDetector;
    unsigned int fStripFront;
    unsigned int fStripBack;
//    unsigned int fPad;
};

#endif
