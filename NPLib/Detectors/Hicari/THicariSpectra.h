#ifndef THICARISPECTRA_H
#define THICARISPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : march 2011                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for Hicari                 *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    + first version (not complete yet)                                     *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "THicariData.h"
#include "THicariPhysics.h"

// ForwardDeclaration
class THicariPhysics ;

class THicariSpectra:public VSpectra {
  public:
    // constructor and destructor
    THicariSpectra();
    THicariSpectra(unsigned int NumberOfClover);
    ~THicariSpectra();

  private:
    // Initialization methods
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  public:
    // Filling methods
    void FillRawSpectra(THicariData*);
    void FillPreTreatedSpectra(THicariData*);
    void FillPhysicsSpectra(THicariPhysics*);

  private: // Information on Hicari
    unsigned int fNumberOfClover ;
    unsigned int fnbinsRaw;
    unsigned int fbinRawMin;
    unsigned int fbinRawMax;
    unsigned int fnbinsCal;
    unsigned int fbinCalMin;
    unsigned int fbinCalMax;
    unsigned int fNumberOfSegments;
    unsigned int fNumberOfCores;
};

#endif
