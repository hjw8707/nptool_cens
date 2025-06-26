#ifndef TPlastic_BEDOSPECTRA_H
#define TPlastic_BEDOSPECTRA_H
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
 *  This class hold Plastic_BEDO Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TPlastic_BEDOData.h"
#include "TPlastic_BEDOPhysics.h"

// Forward Declaration
class TPlastic_BEDOPhysics;


class TPlastic_BEDOSpectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TPlastic_BEDOSpectra();
    TPlastic_BEDOSpectra(unsigned int NumberOfDetectors);
    ~TPlastic_BEDOSpectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TPlastic_BEDOData*);
    void FillPreTreatedSpectra(TPlastic_BEDOData*);
    void FillPhysicsSpectra(TPlastic_BEDOPhysics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters 
  private:
    unsigned int fNumberOfDetectors;
};

#endif
