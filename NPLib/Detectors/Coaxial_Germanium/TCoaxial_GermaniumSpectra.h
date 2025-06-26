#ifndef TCoaxial_GermaniumSPECTRA_H
#define TCoaxial_GermaniumSPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Léo Plagnol  contact address: leo.plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : January 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Coaxial_Germanium Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TCoaxial_GermaniumData.h"
#include "TCoaxial_GermaniumPhysics.h"

// Forward Declaration
class TCoaxial_GermaniumPhysics;


class TCoaxial_GermaniumSpectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TCoaxial_GermaniumSpectra();
    TCoaxial_GermaniumSpectra(unsigned int NumberOfDetectors);
    ~TCoaxial_GermaniumSpectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TCoaxial_GermaniumData*);
    void FillPreTreatedSpectra(TCoaxial_GermaniumData*);
    void FillPhysicsSpectra(TCoaxial_GermaniumPhysics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters 
  private:
    unsigned int fNumberOfDetectors;
};

#endif
