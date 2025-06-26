#ifndef TEdinburghDSSDSPECTRA_H
#define TEdinburghDSSDSPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Louis Heitz  contact address: louis.heitz@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : mars 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold EdinburghDSSD Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TEdinburghDSSDData.h"
#include "TEdinburghDSSDPhysics.h"

// Forward Declaration
class TEdinburghDSSDPhysics;


class TEdinburghDSSDSpectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TEdinburghDSSDSpectra();
    TEdinburghDSSDSpectra(std::map<int,int> DetectorIndex);
    ~TEdinburghDSSDSpectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TEdinburghDSSDData*);
    void FillPreTreatedSpectra(TEdinburghDSSDData*);
    void FillPhysicsSpectra(TEdinburghDSSDPhysics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters
  private:
    std::map<int,int> fDetectorToIndex;
    std::map<int,int> fIndexToDetector;
    unsigned int fNumberOfDetector;
    unsigned int fStripX;
    unsigned int fStripY;
    unsigned int fStripSecondLayer;
};

#endif
