#ifndef TNeuLANDSPECTRA_H
#define TNeuLANDSPECTRA_H
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
 *  This class hold NeuLAND Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TNeuLANDData.h"
#include "TNeuLANDPhysics.h"

// Forward Declaration
class TNeuLANDPhysics;


class TNeuLANDSpectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TNeuLANDSpectra();
    TNeuLANDSpectra(unsigned int NumberOfDetectors);
    ~TNeuLANDSpectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TNeuLANDData*);
    void FillPreTreatedSpectra(TNeuLANDData*);
    void FillPhysicsSpectra(TNeuLANDPhysics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters 
  private:
    unsigned int fNumberOfDetectors;
};

#endif
