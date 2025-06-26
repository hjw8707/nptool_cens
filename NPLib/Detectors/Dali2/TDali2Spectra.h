#ifndef TDali2SPECTRA_H
#define TDali2SPECTRA_H
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elidiano Tronchin  contact address: elidiano.tronchin@studenti.unipd.it                        *
 *                                                                           *
 * Creation Date  : septembre 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Dali Spectra                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// NPLib headers
#include "NPVSpectra.h"
#include "TDali2Data.h"
#include "TDali2Physics.h"

// Forward Declaration
class TDali2Physics;


class TDali2Spectra : public VSpectra {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TDali2Spectra();
    TDali2Spectra(unsigned int NumberOfDetectors);
    ~TDali2Spectra();

  //////////////////////////////////////////////////////////////
  // Initialization methods
  private:
    void InitRawSpectra();
    void InitPreTreatedSpectra();
    void InitPhysicsSpectra();

  //////////////////////////////////////////////////////////////
  // Filling methods
  public:
    void FillRawSpectra(TDali2Data*);
    void FillPreTreatedSpectra(TDali2Data*);
    void FillPhysicsSpectra(TDali2Physics*);

  //////////////////////////////////////////////////////////////
  // Detector parameters 
  private:
    unsigned int fNumberOfDetectors;
};

#endif
