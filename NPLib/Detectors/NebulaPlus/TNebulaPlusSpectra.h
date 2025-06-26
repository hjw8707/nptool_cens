#ifndef TNebulaPlusSPECTRA_H
#define TNebulaPlusSPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : December 2019                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold NebulaPlus Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TNebulaPlusData.h"
#include "TNebulaPlusPhysics.h"

// Forward Declaration
class TNebulaPlusPhysics;


class TNebulaPlusSpectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TNebulaPlusSpectra();
    TNebulaPlusSpectra(unsigned int NumberOfDetectors);
    ~TNebulaPlusSpectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TNebulaPlusData*);
    void FillPreTreatedSpectra(TNebulaPlusData*);
    void FillPhysicsSpectra(TNebulaPlusPhysics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters 
  private:
    unsigned int fNumberOfDetectors;
};

#endif
