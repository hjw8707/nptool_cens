#ifndef TquadranMSQ25SPECTRA_H
#define TquadranMSQ25SPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: leo.plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : January 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold quadranMSQ25 Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TquadranMSQ25Data.h"
#include "TquadranMSQ25Physics.h"

// Forward Declaration
class TquadranMSQ25Physics;


class TquadranMSQ25Spectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TquadranMSQ25Spectra();
    TquadranMSQ25Spectra(unsigned int NumberOfDetectors);
    ~TquadranMSQ25Spectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TquadranMSQ25Data*);
    void FillPreTreatedSpectra(TquadranMSQ25Data*);
    void FillPhysicsSpectra(TquadranMSQ25Physics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters 
  private:
    unsigned int fNumberOfDetectors;
};

#endif
