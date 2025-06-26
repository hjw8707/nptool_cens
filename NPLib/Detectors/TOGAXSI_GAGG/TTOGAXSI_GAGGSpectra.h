#ifndef TTOGAXSI_GAGGSPECTRA_H
#define TTOGAXSI_GAGGSPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Pohl  contact address: thomas.pohl@riken.jp                        *
 *                                                                           *
 * Creation Date  : February 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TOGAXSI_GAGG Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TTOGAXSI_GAGGData.h"
#include "TTOGAXSI_GAGGPhysics.h"

// Forward Declaration
class TTOGAXSI_GAGGPhysics;


class TTOGAXSI_GAGGSpectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TTOGAXSI_GAGGSpectra();
//    TTOGAXSI_GAGGSpectra(unsigned int NumberOfRecoilDetectors,unsigned int NumberOfClusterDetectors);
    TTOGAXSI_GAGGSpectra(unsigned int NumberOfDetectors);
    ~TTOGAXSI_GAGGSpectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TTOGAXSI_GAGGData*);
    void FillPreTreatedSpectra(TTOGAXSI_GAGGData*);
    void FillPhysicsSpectra(TTOGAXSI_GAGGPhysics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters 
  private:
//    unsigned int fNumberOfRecoilDetectors;
//    unsigned int fNumberOfClusterDetectors;
    unsigned int fNumberOfDetectors;
};

#endif
