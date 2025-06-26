#ifndef TMonsterSPECTRA_H
#define TMonsterSPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Emile Cantacuzène  contact address: emile.cantacuzene@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : August 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Monster Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TMonsterData.h"
#include "TMonsterPhysics.h"

// Forward Declaration
class TMonsterPhysics;


class TMonsterSpectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TMonsterSpectra();
    TMonsterSpectra(unsigned int NumberOfDetectors);
    ~TMonsterSpectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TMonsterData*);
    void FillPreTreatedSpectra(TMonsterData*);
    void FillPhysicsSpectra(TMonsterPhysics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters 
  private:
    unsigned int fNumberOfDetectors;
};

#endif
